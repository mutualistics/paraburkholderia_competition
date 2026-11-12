
# Figure 5b: supernatant experiment — processing and statistics
# Run this script before Figure5b.R in the same R session.
# Fixed-date models, planned contrasts, and adjusted p-values live here.

suppressPackageStartupMessages({
  library(tidyverse)
  library(MESS)
  library(emmeans)
})

# Extract the supernatant ZIP, then edit this path to its CSV directory.
# Run from the repository root if using the default relative directory.
DATA_DIR <- "replace/with/your/path/to/supernatant.zip"

CONCENTRATION <- "High"            # 90% supernatant
EXCLUDED_DATES <- c("0511", "1901", "1302")

if (!dir.exists(DATA_DIR)) {
  stop("Supernatant CSV directory not found: ", DATA_DIR,
       "\nExtract the raw-data ZIP and set DATA_DIR to its location.")
}

read_measurement <- function(path, channel) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    tidyr::pivot_longer(-Time, names_to = "Well_ID", values_to = channel) %>%
    mutate(Date = stringr::str_extract(basename(path), "\\d{4}(?=\\.csv$)"))
}
read_metadata <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(Date = stringr::str_extract(basename(path), "\\d{4}(?=\\.csv$)"))
}

files <- list.files(DATA_DIR, pattern = "\\.csv$", full.names = TRUE)
file_groups <- list(
  OD = grep("OD", files, value = TRUE),
  GFP = grep("GFP", files, value = TRUE),
  RFP = grep("RFP", files, value = TRUE),
  Metadata = grep("Metadata", files, value = TRUE)
)
if (any(lengths(file_groups) == 0L)) {
  stop("Missing one or more OD, GFP, RFP, Metadata CSV file groups in ", DATA_DIR)
}

measurements <- purrr::reduce(list(
  purrr::map_dfr(file_groups$OD, read_measurement, channel = "OD"),
  purrr::map_dfr(file_groups$GFP, read_measurement, channel = "GFP"),
  purrr::map_dfr(file_groups$RFP, read_measurement, channel = "RFP")
), dplyr::full_join, by = c("Time", "Well_ID", "Date"))

metadata <- purrr::map_dfr(file_groups$Metadata, read_metadata)
raw <- measurements %>%
  left_join(metadata, by = c("Well_ID", "Date"))

required <- c("Pseudomonas", "Competitor", "Supernatant", "Concentration",
              "Replicate", "Time", "Date", "GFP", "RFP")
if (!all(required %in% names(raw))) {
  stop("Missing required columns: ", paste(setdiff(required, names(raw)), collapse = ", "),
       "\nCheck that DATA_DIR points to the supernatant experiment, not Figure 5a.")
}
raw <- tidyr::drop_na(raw, GFP, RFP, Supernatant)

# Background: GFP measured in same-day wells without Pseudomonas.
background <- raw %>%
  filter(Pseudomonas == "Absent") %>%
  group_by(Date) %>%
  summarise(Background_GFP = mean(GFP, na.rm = TRUE), .groups = "drop")

corrected <- raw %>%
  left_join(background, by = "Date") %>%
  mutate(Corrected_GFP = GFP - Background_GFP)

# Integrate GFP fluorescence per technical well.
auc_wells <- corrected %>%
  group_by(Pseudomonas, Competitor, Supernatant, Concentration, Replicate, Date) %>%
  arrange(Time, .by_group = TRUE) %>%
  summarise(GFP_AUC = MESS::auc(Time, Corrected_GFP), .groups = "drop")

# Normalize to the same-date Pseudomonas monoculture in fresh YG.
daily_control <- auc_wells %>%
  filter(Supernatant == "Fresh YG", Pseudomonas == "Present") %>%
  group_by(Date) %>%
  summarise(Daily_Control_AUC = mean(GFP_AUC, na.rm = TRUE), .groups = "drop")

auc_normalized <- auc_wells %>%
  left_join(daily_control, by = "Date") %>%
  mutate(Normalized_AUC = GFP_AUC / Daily_Control_AUC) %>%
  filter(Supernatant != "Fresh YG") %>%
  mutate(
    Concentration = factor(Concentration, levels = c("Low", "Medium", "High")),
    Supernatant = factor(Supernatant),
    Competitor = factor(Competitor)
  )

# Original assay exclusions and 90% supernatant selection.
# Average log values across technical wells only after excluding nonpositive AUCs.
plot_data_reps1 <- auc_normalized %>%
  filter(!(Competitor == "Pseudomonas" & Supernatant == "Competition"),
         Pseudomonas != "Absent", !(Date %in% EXCLUDED_DATES),
         Concentration == CONCENTRATION) %>%
  mutate(
    Supernatant = forcats::fct_relevel(Supernatant, "Monoculture", "Competition"),
    Normalized_AUC_pos = if_else(Normalized_AUC > 0, Normalized_AUC, NA_real_),
    log10_Fitness = log10(Normalized_AUC_pos)
  ) %>%
  drop_na(log10_Fitness)

# The original response is log10 of the date-mean normalized AUC.
# Mean_log10_Fitness is retained separately as a descriptive well-level summary.
plot_data1 <- plot_data_reps1 %>%
  group_by(Date, Pseudomonas, Competitor, Supernatant, Concentration) %>%
  summarise(
    Mean_Normalized_AUC = mean(Normalized_AUC, na.rm = TRUE),
    Mean_log10_Fitness = mean(log10_Fitness, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Mean_Normalized_AUC_pos = if_else(Mean_Normalized_AUC > 0,
                                      Mean_Normalized_AUC, NA_real_),
    log10_Fitness = log10(Mean_Normalized_AUC_pos)
  ) %>%
  drop_na(log10_Fitness)

cat("\nFigure 5b: biological dates by treatment\n")
print(count(plot_data1, Competitor, Supernatant), n = Inf)

model_dat <- plot_data1 %>%
  mutate(
    Competitor_chr = as.character(Competitor),
    Species = case_when(
      stringr::str_detect(Competitor_chr,
                          stringr::regex("^Ba|agricolaris|Agricolaris", ignore_case = TRUE)) ~
        "P. agricolaris",
      stringr::str_detect(Competitor_chr,
                          stringr::regex("^Bh|hayleyella|Hayleyella", ignore_case = TRUE)) ~
        "P. hayleyella",
      stringr::str_detect(Competitor_chr,
                          stringr::regex("Pseudomonas", ignore_case = TRUE)) ~ "Pseudomonas",
      TRUE ~ Competitor_chr
    ),
    Species = factor(Species,
                     levels = c("P. agricolaris", "P. hayleyella", "Pseudomonas")),
    Supernatant = factor(Supernatant, levels = c("Monoculture", "Competition")),
    Date = factor(Date),
    Competitor = factor(Competitor)
  )

cat("\nSpecies mapping\n")
print(distinct(model_dat, Competitor, Species), n = Inf)

# Factorial model: effects of symbiont species and supernatant source.
main_dat <- model_dat %>%
  filter(Species != "Pseudomonas") %>%
  drop_na(log10_Fitness, Species, Supernatant, Date)
m_main_date_fixed <- lm(log10_Fitness ~ Species * Supernatant + Date,
                        data = main_dat)
emm_main <- emmeans(m_main_date_fixed, ~ Species * Supernatant)
within_species_factorial <- contrast(
  emmeans(m_main_date_fixed, ~ Supernatant | Species),
  method = "pairwise", adjust = "none"
)

cat("\nFactorial model\n")
print(anova(m_main_date_fixed))
print(drop1(m_main_date_fixed, test = "F"))
print(emm_main)
print(within_species_factorial)

# Treatment model includes the Pseudomonas monoculture baseline.
treatment_levels <- c(
  "Pseudomonas_mono", "Agricolaris_mono", "Agricolaris_coculture",
  "Hayleyella_mono", "Hayleyella_coculture"
)
treatment_dat <- model_dat %>%
  mutate(
    Treatment = case_when(
      Species == "P. agricolaris" & Supernatant == "Monoculture" ~ "Agricolaris_mono",
      Species == "P. agricolaris" & Supernatant == "Competition" ~ "Agricolaris_coculture",
      Species == "P. hayleyella" & Supernatant == "Monoculture" ~ "Hayleyella_mono",
      Species == "P. hayleyella" & Supernatant == "Competition" ~ "Hayleyella_coculture",
      Species == "Pseudomonas" ~ "Pseudomonas_mono",
      TRUE ~ NA_character_
    ),
    Treatment = factor(Treatment, levels = treatment_levels),
    Date = factor(Date)
  ) %>%
  drop_na(Treatment, log10_Fitness)

m_treatment <- lm(log10_Fitness ~ Treatment + Date, data = treatment_dat)
emm_treatment <- emmeans(m_treatment, ~ Treatment)
cat("\nTreatment model\n")
print(anova(m_treatment))
print(drop1(m_treatment, test = "F"))
print(emm_treatment)

# Planned comparisons. Holm adjustment is applied separately to:
# (1) four comparisons against Pseudomonas monoculture;
# (2) two within-symbiont comparisons (coculture vs monoculture).
actual_levels <- as.character(as.data.frame(emm_treatment)$Treatment)
if (length(actual_levels) != 5L ||
    !setequal(actual_levels, treatment_levels)) {
  stop("Expected five treatments in the marginal-mean grid.")
}

weights_for <- function(positive, negative) {
  w <- setNames(rep(0, length(treatment_levels)), treatment_levels)
  w[positive] <- 1
  w[negative] <- -1
  unname(w[actual_levels])
}

baseline_methods <- list(
  "Ba monoculture - Pk monoculture" =
    weights_for("Agricolaris_mono", "Pseudomonas_mono"),
  "Ba-Pk coculture - Pk monoculture" =
    weights_for("Agricolaris_coculture", "Pseudomonas_mono"),
  "Bh monoculture - Pk monoculture" =
    weights_for("Hayleyella_mono", "Pseudomonas_mono"),
  "Bh-Pk coculture - Pk monoculture" =
    weights_for("Hayleyella_coculture", "Pseudomonas_mono")
)
within_methods <- list(
  "Ba-Pk coculture - Ba monoculture" =
    weights_for("Agricolaris_coculture", "Agricolaris_mono"),
  "Bh-Pk coculture - Bh monoculture" =
    weights_for("Hayleyella_coculture", "Hayleyella_mono")
)

planned_tests <- function(methods) {
  x <- as.data.frame(summary(
    contrast(emm_treatment, method = methods, adjust = "none"),
    infer = c(TRUE, TRUE), adjust = "none"
  ))
  if (nrow(x) != length(methods) || any(!is.finite(x$p.value))) {
    stop("Incomplete planned contrasts.")
  }
  x$p_holm <- p.adjust(x$p.value, method = "holm")
  x
}

baseline_df <- planned_tests(baseline_methods)
within_df <- planned_tests(within_methods)
cat("\nPlanned contrasts vs Pseudomonas monoculture (Holm across 4)\n")
print(baseline_df, row.names = FALSE)
cat("\nWithin-species contrasts (Holm across 2)\n")
print(within_df, row.names = FALSE)

get_p <- function(tbl, label) {
  p <- tbl$p_holm[tbl$contrast == label]
  if (length(p) != 1L) stop("Missing planned contrast: ", label)
  p
}
star_label <- function(p) ifelse(
  p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))
)

fig5b_bracket_labels <- c(
  baseline_Ba_coculture = star_label(
    get_p(baseline_df, "Ba-Pk coculture - Pk monoculture")
  ),
  within_Ba = star_label(
    get_p(within_df, "Ba-Pk coculture - Ba monoculture")
  ),
  within_Bh = star_label(
    get_p(within_df, "Bh-Pk coculture - Bh monoculture")
  )
)

# Provide model EMMs and confidence intervals for the plotting script.
fig5b_plot_emm <- as.data.frame(emm_treatment)
ci_lo <- intersect(c("lower.CL", "asymp.LCL", "LCL"), names(fig5b_plot_emm))[1]
ci_hi <- intersect(c("upper.CL", "asymp.UCL", "UCL"), names(fig5b_plot_emm))[1]
if (is.na(ci_lo) || is.na(ci_hi)) {
  stop("Marginal-mean confidence intervals were not returned.")
}
fig5b_plot_emm$lower_ci <- fig5b_plot_emm[[ci_lo]]
fig5b_plot_emm$upper_ci <- fig5b_plot_emm[[ci_hi]]

cat("\nFigure 5b model and planned contrasts are ready. Run Figure5b.R.\n")
