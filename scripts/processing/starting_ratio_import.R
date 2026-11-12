
# Figure 5a: starting-ratio experiment — data processing and summaries
# Run this script before Figure5a.R in the same R session.
# This figure contains descriptive estimates only (no hypothesis tests).
# Technical wells are averaged within dates; SE/CI95 is across dates.

suppressPackageStartupMessages(library(tidyverse))

# Extract the starting-ratio ZIP, then edit this path to its CSV directory.
# Run from the repository root if using the default relative directory.
DATA_DIR <- "/Users/calum/Documents/DocumentsMacBook/LabStuff/Elle/PK_Ba_comp"

TIME_MAX <- 24
TIME_SHIFT_HOURS <- 1       # Instrument times start at 1 h; display 0 h at that point.
UNCERTAINTY <- "SE"        # "SE" or "CI95", computed across biological dates

PSEUDOMONAS_START <- c(50, 60, 75, 90, 99, 99.9, 99.99, 100)
BA_LABELS <- c("0%", "0.01%", "0.1%", "1%", "10%", "25%", "40%", "50%")
if (!UNCERTAINTY %in% c("SE", "CI95")) stop("UNCERTAINTY must be SE or CI95.")
if (!dir.exists(DATA_DIR)) {
  stop("Starting-ratio CSV directory not found: ", DATA_DIR,
       "\nExtract the raw-data ZIP and set DATA_DIR to its location.")
}
if (TIME_MAX <= TIME_SHIFT_HOURS) stop("TIME_MAX must exceed TIME_SHIFT_HOURS.")

read_measurement <- function(path, channel) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    tidyr::pivot_longer(-Time, names_to = "Well_ID", values_to = channel) %>%
    mutate(Date = stringr::str_extract(basename(path), "\\d{4}(?=\\.csv$)"))
}
read_metadata <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(Date = stringr::str_extract(basename(path), "\\d{4}(?=\\.csv$)"))
}

files <- list.files(DATA_DIR, pattern = "\\.csv$", full.names = TRUE,
                    ignore.case = TRUE)
find_files <- function(prefix) {
  files[stringr::str_detect(basename(files),
                            stringr::regex(paste0("^", prefix), ignore_case = TRUE))]
}
file_groups <- list(
  OD = find_files("OD"), GFP = find_files("GFP"),
  RFP = find_files("RFP"), Metadata = find_files("Metadata")
)
if (any(lengths(file_groups) == 0L)) {
  stop("Missing one or more OD, GFP, RFP, Metadata CSV file groups in ", DATA_DIR)
}

measurements <- purrr::reduce(list(
  purrr::map_dfr(file_groups$OD, read_measurement, channel = "OD"),
  purrr::map_dfr(file_groups$GFP, read_measurement, channel = "GFP"),
  purrr::map_dfr(file_groups$RFP, read_measurement, channel = "RFP")
), full_join, by = c("Time", "Well_ID", "Date"))

metadata <- purrr::map_dfr(file_groups$Metadata, read_metadata)
raw <- measurements %>% left_join(metadata, by = c("Well_ID", "Date"))
required <- c("Pseudomonas", "Condition", "Pseudo_conc", "GFP", "RFP")
if (!all(required %in% names(raw))) {
  stop("Missing required columns: ", paste(setdiff(required, names(raw)), collapse = ", "))
}

# Pfrac is the starting percentage of Pseudomonas in the inoculum.
fig5a_data <- raw %>%
  drop_na(GFP, RFP) %>%
  mutate(
    Time = as.numeric(Time),
    Pseudomonas = trimws(as.character(Pseudomonas)),
    Condition = trimws(as.character(Condition)),
    Pseudo_conc = readr::parse_number(as.character(Pseudo_conc)),
    Date = as.character(Date),
    Pfrac = case_when(
      Condition == "Media Control" & Pseudomonas == "Absent" ~ NA_real_,
      Condition == "Control" & Pseudomonas == "Absent" ~ 0,
      Condition == "Control" & Pseudomonas == "Present" ~ 100,
      Condition == "Treatment" & Pseudomonas == "Present" ~ Pseudo_conc,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(Pfrac))

# Media background is the same-date mean of media-control wells.
media <- raw %>%
  mutate(Pseudomonas = trimws(as.character(Pseudomonas)),
         Condition = trimws(as.character(Condition)), Date = as.character(Date)) %>%
  filter(Condition == "Media Control", Pseudomonas == "Absent") %>%
  group_by(Date) %>%
  summarise(Media_GFP = mean(GFP, na.rm = TRUE),
            Media_RFP = mean(RFP, na.rm = TRUE), .groups = "drop")

fig5a_data <- fig5a_data %>%
  left_join(media, by = "Date") %>%
  mutate(
    Media_GFP = coalesce(Media_GFP, 0),
    Media_RFP = coalesce(Media_RFP, 0),
    GFP_bg = pmax(GFP - Media_GFP, 0),
    RFP_bg = pmax(RFP - Media_RFP, 0)
  )

# Estimate fluorescence cross-talk separately for each experimental date.
# Pfrac = 0 is Ba-only; Pfrac = 100 is Pseudomonas-only.
slopes_rg <- fig5a_data %>%
  filter(Pfrac == 0) %>% group_by(Date) %>%
  summarise(beta_RG = {
    z <- pick(GFP_bg, RFP_bg)
    if (n() >= 3 && any(z$RFP_bg > 0, na.rm = TRUE)) {
      as.numeric(coef(lm(GFP_bg ~ 0 + RFP_bg, data = z))[1])
    } else 0
  }, .groups = "drop")
slopes_gr <- fig5a_data %>%
  filter(Pfrac == 100) %>% group_by(Date) %>%
  summarise(beta_GR = {
    z <- pick(RFP_bg, GFP_bg)
    if (n() >= 3 && any(z$GFP_bg > 0, na.rm = TRUE)) {
      as.numeric(coef(lm(RFP_bg ~ 0 + GFP_bg, data = z))[1])
    } else 0
  }, .groups = "drop")

slopes <- full_join(slopes_rg, slopes_gr, by = "Date") %>%
  mutate(beta_RG = coalesce(beta_RG, 0),
         beta_GR = coalesce(beta_GR, 0))
fig5a_data <- fig5a_data %>%
  left_join(slopes, by = "Date") %>%
  mutate(Corrected_GFP = pmax(GFP_bg - beta_RG * RFP_bg, 0),
         Corrected_RFP = pmax(RFP_bg - beta_GR * GFP_bg, 0))

# Normalize each channel to its same-day monoculture maximum.
# This figure displays GFP only; the RFP correction is retained for the audit.
controls <- fig5a_data %>% group_by(Date) %>%
  summarise(
    GFP_max = max(ifelse(Pfrac == 100, Corrected_GFP, NA_real_), na.rm = TRUE),
    RFP_max = max(ifelse(Pfrac == 0, Corrected_RFP, NA_real_), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    GFP_max = if_else(is.finite(GFP_max) & GFP_max > 0, GFP_max, NA_real_),
    RFP_max = if_else(is.finite(RFP_max) & RFP_max > 0, RFP_max, NA_real_)
  )
fig5a_data <- fig5a_data %>%
  left_join(controls, by = "Date") %>%
  mutate(Norm_GFP = if_else(!is.na(GFP_max), Corrected_GFP / GFP_max, NA_real_),
         Norm_RFP = if_else(!is.na(RFP_max), Corrected_RFP / RFP_max, NA_real_))

raw_time <- range(fig5a_data$Time[is.finite(fig5a_data$Time)])
if (any(!is.finite(raw_time))) stop("No usable measurement times.")
if (!isTRUE(all.equal(raw_time[1], TIME_SHIFT_HOURS))) {
  warning("The earliest time is ", raw_time[1],
          " h but TIME_SHIFT_HOURS is ", TIME_SHIFT_HOURS, " h.")
}

label_ba <- function(x) case_when(
  near(x, 50) ~ "50%", near(x, 60) ~ "40%", near(x, 75) ~ "25%",
  near(x, 90) ~ "10%", near(x, 99) ~ "1%", near(x, 99.9) ~ "0.1%",
  near(x, 99.99) ~ "0.01%", near(x, 100) ~ "0%",
  TRUE ~ NA_character_
)

# Each date is a biological replicate; technical wells are averaged first.
fig5a_day <- fig5a_data %>%
  filter(is.finite(Norm_GFP), is.finite(Time), Time <= TIME_MAX) %>%
  mutate(
    Time_plot = Time - TIME_SHIFT_HOURS,
    Pfrac_clean = round(Pfrac, 4),
    Pfrac_clean = case_when(
      near(Pfrac_clean, 50) ~ 50, near(Pfrac_clean, 60) ~ 60,
      near(Pfrac_clean, 75) ~ 75, near(Pfrac_clean, 90) ~ 90,
      near(Pfrac_clean, 99) ~ 99, near(Pfrac_clean, 99.9) ~ 99.9,
      near(Pfrac_clean, 99.99) ~ 99.99, near(Pfrac_clean, 100) ~ 100,
      TRUE ~ Pfrac_clean
    )
  ) %>%
  filter(Pfrac_clean %in% PSEUDOMONAS_START, Time_plot >= 0) %>%
  mutate(
    Agricolaris_start = 100 - Pfrac_clean,
    Curve_group = factor(label_ba(Pfrac_clean), levels = BA_LABELS),
    Curve_role = case_when(
      Pfrac_clean < 90 ~ "high_agricolaris",
      Pfrac_clean %in% c(90, 99, 99.9) ~ "collapse",
      Pfrac_clean >= 99.99 ~ "stable", TRUE ~ "other"
    )
  ) %>%
  group_by(Date, Pfrac_clean, Agricolaris_start,
           Curve_group, Curve_role, Time_plot) %>%
  summarise(day_mean_Pk = mean(Norm_GFP), .groups = "drop")

if (!nrow(fig5a_day) || anyNA(fig5a_day$Curve_group)) {
  stop("No valid plotting curves or unrecognised starting ratios.")
}

fig5a_mean <- fig5a_day %>%
  group_by(Pfrac_clean, Agricolaris_start, Curve_group, Curve_role, Time_plot) %>%
  summarise(n_dates = sum(is.finite(day_mean_Pk)),
            mean_Pk = mean(day_mean_Pk, na.rm = TRUE),
            sd_Pk = sd(day_mean_Pk, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    se_Pk = if_else(n_dates >= 2, sd_Pk / sqrt(n_dates), NA_real_),
    ci95_Pk = if_else(n_dates >= 2, qt(0.975, n_dates - 1) * se_Pk, NA_real_),
    err = if (UNCERTAINTY == "SE") se_Pk else ci95_Pk,
    ymin = mean_Pk - err,
    ymax = mean_Pk + err
  )
error_data <- filter(fig5a_mean, is.finite(err), is.finite(ymin), is.finite(ymax))

cat("\nFigure 5a: biological dates per curve\n")
cat("Uncertainty:", UNCERTAINTY, "\n")

# Objects for Figure5a.R: fig5a_mean, error_data, BA_LABELS,
# TIME_MAX, TIME_SHIFT_HOURS and UNCERTAINTY.
cat("Figure 5a summaries are ready. Run Figure5a.R.\n")
