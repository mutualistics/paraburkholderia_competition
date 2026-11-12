
# Figure 5a: starting-ratio experiment — data processing and summaries
# Run before Figure5a.R and FigureS_starting_ratios.R in the same R session.
# Both figures contain descriptive estimates only (no hypothesis tests).
# Technical wells are averaged within dates; SE/CI95 is across dates.

suppressPackageStartupMessages(library(tidyverse))

# Extract the starting-ratio ZIP, then edit this path to its CSV directory.
# Run from the repository root if using the default relative directory.
DATA_DIR <- "path/to/starting_ratio.zip/data"

TIME_MAX <- 24
TIME_SHIFT_HOURS <- 1
UNCERTAINTY <- "SE"        # "SE" or "CI95", computed across biological dates

PSEUDOMONAS_START <- c(50, 60, 75, 90, 99, 99.9, 99.99, 100)
BA_LABELS <- c("0%", "0.01%", "0.1%", "1%", "10%", "25%", "40%", "50%")

if (!UNCERTAINTY %in% c("SE", "CI95")) {
  stop("UNCERTAINTY must be SE or CI95.")
}

if (!dir.exists(DATA_DIR)) {
  stop(
    "Starting-ratio CSV directory not found: ", DATA_DIR,
    "\nExtract the raw-data ZIP and set DATA_DIR to its location."
  )
}

if (TIME_MAX <= TIME_SHIFT_HOURS) {
  stop("TIME_MAX must exceed TIME_SHIFT_HOURS.")
}

# Import fluorescence data

read_measurement <- function(path, channel) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    tidyr::pivot_longer(
      -Time,
      names_to = "Well_ID",
      values_to = channel
    ) %>%
    mutate(
      Date = stringr::str_extract(
        basename(path), "\\d{4}(?=\\.csv$)"
      )
    )
}

read_metadata <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(
      Date = stringr::str_extract(
        basename(path), "\\d{4}(?=\\.csv$)"
      )
    )
}

files <- list.files(
  DATA_DIR,
  pattern = "\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)

find_files <- function(prefix) {
  files[stringr::str_detect(
    basename(files),
    stringr::regex(paste0("^", prefix), ignore_case = TRUE)
  )]
}

file_groups <- list(
  OD = find_files("OD"),
  GFP = find_files("GFP"),
  RFP = find_files("RFP"),
  Metadata = find_files("Metadata")
)

if (any(lengths(file_groups) == 0L)) {
  stop(
    "Missing one or more OD, GFP, RFP, Metadata CSV groups in ",
    DATA_DIR
  )
}

measurements <- purrr::reduce(
  list(
    purrr::map_dfr(file_groups$OD, read_measurement, channel = "OD"),
    purrr::map_dfr(file_groups$GFP, read_measurement, channel = "GFP"),
    purrr::map_dfr(file_groups$RFP, read_measurement, channel = "RFP")
  ),
  full_join,
  by = c("Time", "Well_ID", "Date")
)

metadata <- purrr::map_dfr(file_groups$Metadata, read_metadata)

raw <- measurements %>%
  left_join(metadata, by = c("Well_ID", "Date"))

required <- c(
  "Pseudomonas", "Condition", "Pseudo_conc", "GFP", "RFP"
)

if (!all(required %in% names(raw))) {
  stop(
    "Missing required columns: ",
    paste(setdiff(required, names(raw)), collapse = ", ")
  )
}

# Starting percentage of Pseudomonas

fig5a_data <- raw %>%
  drop_na(GFP, RFP) %>%
  mutate(
    Time = as.numeric(Time),
    Pseudomonas = trimws(as.character(Pseudomonas)),
    Condition = trimws(as.character(Condition)),
    Pseudo_conc = readr::parse_number(as.character(Pseudo_conc)),
    Date = as.character(Date),
    Pfrac = case_when(
      Condition == "Media Control" &
        Pseudomonas == "Absent" ~ NA_real_,
      Condition == "Control" &
        Pseudomonas == "Absent" ~ 0,
      Condition == "Control" &
        Pseudomonas == "Present" ~ 100,
      Condition == "Treatment" &
        Pseudomonas == "Present" ~ Pseudo_conc,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(Pfrac))

# Same-date media background

media <- raw %>%
  mutate(
    Pseudomonas = trimws(as.character(Pseudomonas)),
    Condition = trimws(as.character(Condition)),
    Date = as.character(Date)
  ) %>%
  filter(
    Condition == "Media Control",
    Pseudomonas == "Absent"
  ) %>%
  group_by(Date) %>%
  summarise(
    Media_GFP = mean(GFP, na.rm = TRUE),
    Media_RFP = mean(RFP, na.rm = TRUE),
    .groups = "drop"
  )

fig5a_data <- fig5a_data %>%
  left_join(media, by = "Date") %>%
  mutate(
    Media_GFP = coalesce(Media_GFP, 0),
    Media_RFP = coalesce(Media_RFP, 0),
    GFP_bg = pmax(GFP - Media_GFP, 0),
    RFP_bg = pmax(RFP - Media_RFP, 0)
  )

# Fluorescence cross-talk, estimated separately for each date.
# Pfrac = 0 is Ba-only; Pfrac = 100 is Pseudomonas-only.

slopes_rg <- fig5a_data %>%
  filter(Pfrac == 0) %>%
  group_by(Date) %>%
  summarise(
    beta_RG = {
      z <- pick(GFP_bg, RFP_bg)
      if (n() >= 3 && any(z$RFP_bg > 0, na.rm = TRUE)) {
        as.numeric(coef(lm(GFP_bg ~ 0 + RFP_bg, data = z))[1])
      } else 0
    },
    .groups = "drop"
  )

slopes_gr <- fig5a_data %>%
  filter(Pfrac == 100) %>%
  group_by(Date) %>%
  summarise(
    beta_GR = {
      z <- pick(RFP_bg, GFP_bg)
      if (n() >= 3 && any(z$GFP_bg > 0, na.rm = TRUE)) {
        as.numeric(coef(lm(RFP_bg ~ 0 + GFP_bg, data = z))[1])
      } else 0
    },
    .groups = "drop"
  )

slopes <- full_join(slopes_rg, slopes_gr, by = "Date") %>%
  mutate(
    beta_RG = coalesce(beta_RG, 0),
    beta_GR = coalesce(beta_GR, 0)
  )

fig5a_data <- fig5a_data %>%
  left_join(slopes, by = "Date") %>%
  mutate(
    Corrected_GFP = pmax(GFP_bg - beta_RG * RFP_bg, 0),
    Corrected_RFP = pmax(RFP_bg - beta_GR * GFP_bg, 0)
  )

# Same-day monoculture maximum normalization

controls <- fig5a_data %>%
  group_by(Date) %>%
  summarise(
    GFP_max = max(
      ifelse(Pfrac == 100, Corrected_GFP, NA_real_),
      na.rm = TRUE
    ),
    RFP_max = max(
      ifelse(Pfrac == 0, Corrected_RFP, NA_real_),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    GFP_max = if_else(
      is.finite(GFP_max) & GFP_max > 0,
      GFP_max, NA_real_
    ),
    RFP_max = if_else(
      is.finite(RFP_max) & RFP_max > 0,
      RFP_max, NA_real_
    )
  )

fig5a_data <- fig5a_data %>%
  left_join(controls, by = "Date") %>%
  mutate(
    Norm_GFP = if_else(
      !is.na(GFP_max),
      Corrected_GFP / GFP_max,
      NA_real_
    ),
    Norm_RFP = if_else(
      !is.na(RFP_max),
      Corrected_RFP / RFP_max,
      NA_real_
    )
  )

raw_time <- range(fig5a_data$Time[is.finite(fig5a_data$Time)])

if (any(!is.finite(raw_time))) {
  stop("No usable measurement times.")
}

if (!isTRUE(all.equal(raw_time[1], TIME_SHIFT_HOURS))) {
  warning(
    "The earliest time is ", raw_time[1],
    " h but TIME_SHIFT_HOURS is ", TIME_SHIFT_HOURS, " h."
  )
}

# Main Figure 5a: selected eight ratios

label_ba <- function(x) case_when(
  near(x, 50) ~ "50%",
  near(x, 60) ~ "40%",
  near(x, 75) ~ "25%",
  near(x, 90) ~ "10%",
  near(x, 99) ~ "1%",
  near(x, 99.9) ~ "0.1%",
  near(x, 99.99) ~ "0.01%",
  near(x, 100) ~ "0%",
  TRUE ~ NA_character_
)

fig5a_day <- fig5a_data %>%
  filter(
    is.finite(Norm_GFP),
    is.finite(Time),
    Time <= TIME_MAX
  ) %>%
  mutate(
    Time_plot = Time - TIME_SHIFT_HOURS,
    Pfrac_clean = round(Pfrac, 4),
    Pfrac_clean = case_when(
      near(Pfrac_clean, 50) ~ 50,
      near(Pfrac_clean, 60) ~ 60,
      near(Pfrac_clean, 75) ~ 75,
      near(Pfrac_clean, 90) ~ 90,
      near(Pfrac_clean, 99) ~ 99,
      near(Pfrac_clean, 99.9) ~ 99.9,
      near(Pfrac_clean, 99.99) ~ 99.99,
      near(Pfrac_clean, 100) ~ 100,
      TRUE ~ Pfrac_clean
    )
  ) %>%
  filter(
    Pfrac_clean %in% PSEUDOMONAS_START,
    Time_plot >= 0
  ) %>%
  mutate(
    Agricolaris_start = 100 - Pfrac_clean,
    Curve_group = factor(
      label_ba(Pfrac_clean),
      levels = BA_LABELS
    ),
    Curve_role = case_when(
      Pfrac_clean < 90 ~ "high_agricolaris",
      Pfrac_clean %in% c(90, 99, 99.9) ~ "collapse",
      Pfrac_clean >= 99.99 ~ "stable",
      TRUE ~ "other"
    )
  ) %>%
  group_by(
    Date, Pfrac_clean, Agricolaris_start,
    Curve_group, Curve_role, Time_plot
  ) %>%
  summarise(
    day_mean_Pk = mean(Norm_GFP),
    .groups = "drop"
  )

if (!nrow(fig5a_day) || anyNA(fig5a_day$Curve_group)) {
  stop("No valid Figure 5a plotting curves.")
}

fig5a_mean <- fig5a_day %>%
  group_by(
    Pfrac_clean, Agricolaris_start,
    Curve_group, Curve_role, Time_plot
  ) %>%
  summarise(
    n_dates = sum(is.finite(day_mean_Pk)),
    mean_Pk = mean(day_mean_Pk, na.rm = TRUE),
    sd_Pk = sd(day_mean_Pk, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    se_Pk = if_else(
      n_dates >= 2,
      sd_Pk / sqrt(n_dates),
      NA_real_
    ),
    ci95_Pk = if_else(
      n_dates >= 2,
      qt(0.975, n_dates - 1) * se_Pk,
      NA_real_
    ),
    err = if (UNCERTAINTY == "SE") se_Pk else ci95_Pk,
    ymin = mean_Pk - err,
    ymax = mean_Pk + err
  )

error_data <- filter(
  fig5a_mean,
  is.finite(err),
  is.finite(ymin),
  is.finite(ymax)
)

cat("\nFigure 5a: biological dates per selected curve\n")

print(
  fig5a_day %>%
    group_by(Curve_group) %>%
    summarise(n_dates = n_distinct(Date), .groups = "drop"),
  n = Inf
)

cat("Uncertainty:", UNCERTAINTY, "\n")

# Supplementary figure: both fluorescent species across all 15 starting ratios.
# The summaries are independent of the eight-ratio main-figure summaries.
# Each experimental date, not each technical well, is the replicate.

# Do not plot fluorescence from a species absent from the inoculum.

fig5a_all_day <- fig5a_data %>%
  filter(
    is.finite(Pfrac),
    is.finite(Time),
    !is.na(Date),
    Time <= TIME_MAX
  ) %>%
  mutate(
    Time_plot = Time - TIME_SHIFT_HOURS,
    Pfrac_clean = signif(Pfrac, digits = 10)
  ) %>%
  filter(Time_plot >= 0) %>%
  pivot_longer(
    cols = c(Norm_GFP, Norm_RFP),
    names_to = "Channel",
    names_prefix = "Norm_",
    values_to = "Normalized_fluorescence"
  ) %>%
  filter(
    (Channel == "GFP" & Pfrac_clean > 0) |
      (Channel == "RFP" & Pfrac_clean < 100),
    is.finite(Normalized_fluorescence)
  ) %>%
  mutate(
    Species = if_else(
      Channel == "GFP",
      "Pseudomonas",
      "P. agricolaris"
    )
  ) %>%
  group_by(Date, Pfrac_clean, Time_plot, Species) %>%
  summarise(
    day_mean = mean(Normalized_fluorescence),
    .groups = "drop"
  )

# Format starting P. agricolaris percentages for the facet headings.

format_ratio_pct <- function(x) {
  x <- pmax(0, round(as.numeric(x), digits = 10))
  lbl <- sprintf("%.10f", x)
  lbl <- sub("0+$", "", lbl)
  lbl <- sub("\\.$", "", lbl)
  ifelse(x == 0, "0", lbl)
}

fig5a_ratio_key <- fig5a_all_day %>%
  distinct(Pfrac_clean) %>%
  mutate(
    Agricolaris_start = pmax(0, 100 - Pfrac_clean),
    Ratio_label = paste0(
      format_ratio_pct(Agricolaris_start), "%"
    ),
    In_Fig5a = Pfrac_clean %in% PSEUDOMONAS_START
  ) %>%
  arrange(Agricolaris_start)

if (
  nrow(fig5a_ratio_key) != 15L ||
  n_distinct(fig5a_ratio_key$Ratio_label) != 15L
) {
  print(fig5a_ratio_key, n = Inf)
  stop("Expected 15 distinct starting ratios in the supplementary data.")
}

if (sum(fig5a_ratio_key$In_Fig5a) != 8L) {
  stop("Expected eight of the 15 ratios to be included in Figure 5a.")
}

fig5a_all_day <- fig5a_all_day %>%
  left_join(fig5a_ratio_key, by = "Pfrac_clean") %>%
  mutate(
    Ratio_label = factor(
      Ratio_label,
      levels = fig5a_ratio_key$Ratio_label
    ),
    Species = factor(
      Species,
      levels = c("P. agricolaris", "Pseudomonas")
    )
  )

fig5a_all_mean <- fig5a_all_day %>%
  group_by(
    Pfrac_clean, Agricolaris_start, Ratio_label,
    Species, In_Fig5a, Time_plot
  ) %>%
  summarise(
    n_dates = n(),
    mean_fluorescence = mean(day_mean),
    sd_fluorescence = sd(day_mean),
    .groups = "drop"
  ) %>%
  mutate(
    se_fluorescence = if_else(
      n_dates >= 2L,
      sd_fluorescence / sqrt(n_dates),
      NA_real_
    ),
    ci95_fluorescence = if_else(
      n_dates >= 2L,
      qt(0.975, df = n_dates - 1L) * se_fluorescence,
      NA_real_
    ),
    err = if (UNCERTAINTY == "SE") {
      se_fluorescence
    } else {
      ci95_fluorescence
    },
    ymin = pmax(0, mean_fluorescence - err),
    ymax = mean_fluorescence + err
  )

# Check that both channels are present unless a species was absent.

ratio_species <- fig5a_all_day %>%
  distinct(Pfrac_clean, Species) %>%
  count(Pfrac_clean, name = "n_species") %>%
  left_join(fig5a_ratio_key, by = "Pfrac_clean")

if (any(
  ratio_species$n_species != if_else(
    ratio_species$Pfrac_clean %in% c(0, 100),
    1L,
    2L
  )
)) {
  stop("Unexpected missing channel in a starting-ratio condition.")
}

cat("\nSupplementary: biological dates per ratio and species:\n")

print(
  fig5a_all_day %>%
    group_by(Ratio_label, Species) %>%
    summarise(
      n_dates = n_distinct(Date),
      .groups = "drop"
    ),
  n = Inf
)

cat(
  "\nAll 15 two-channel summaries ready. ",
  "Run FigureS_starting_ratios.R.\n",
  sep = ""
)
