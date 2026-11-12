# Competition plate-reader data import
# Run from the project directory containing ExperimentData/.
# Correct fluorescence, calculate monoculture-normalized AUCs, and
# prepare the date-level log10(RFP/GFP) response for Figures 2 and 3.
# All technically valid experimental dates are included.
# Technical replicates are averaged within each experimental date.
# Outputs: analysis_data_rep_all_valid.rds and analysis_data_day_all_valid.rds.

suppressPackageStartupMessages({
  library(tidyverse)
  library(MESS)
})

# Paths / inputs
# Run from the project root, with the raw CSV files in ExperimentData/.

# Directory containing the extracted raw experimental CSV files.
data_dir <- "/path/to/extracted/competition_data.zip"


if (!dir.exists(data_dir)) stop("Data folder not found: ", data_dir)

OD_FLOOR <- 1e-4
SAVE_DATE_COUNTS_XLSX <- FALSE

# AUC integration windows
# Ss is truncated at 36 h because of late atypical dynamics.
# All other competitors are integrated to 48 h.
AUC_T_END_DEFAULT <- 48
AUC_T_END_SS <- 36

# Correct Ss×wt control dates
# Keep ONLY these for GFP_Strain == "Ss" & RFP_Strain == "wt"
ss_wt_dates <- c("2002","1502","3001","1802","1002","1210","0802","1202","2909")

# Strains
rfp_strains <- c(
  "Pg", "Pu", "Pf", "Pu_wt", "Pf_wt", "Ps", "Ps_wt",
  "Bh22", "Bb862", "Bh171", "Bh11", "Ba31", "Ba161", "Ba159",
  "Bb859", "Ba80", "Bh21", "Bb433", "Bb822"
)

gfp_strains <- c("Pk", "Ss", "Sl", "Bp", "At", "Ck")

# Order for plotting / downstream factor consistency
gfp_levels <- c(
  "None", "wt", "Pk", "Sl", "Ss", "Bp", "Ck", "At",
  "Pkwt1", "Pkwt2", "Bh11", "Ba31", "Ss_wt"
)

# I/O helpers
read_measurement_data <- function(file_path, measurement_type) {
  readr::read_csv(file_path, show_col_types = FALSE) %>%
    tidyr::pivot_longer(
      cols = -Time,
      names_to = "Well_ID",
      values_to = measurement_type
    ) %>%
    dplyr::mutate(Date = gsub(".*(\\d{4}).csv$", "\\1", basename(file_path)))
}

read_metadata <- function(file_path) {
  readr::read_csv(file_path, show_col_types = FALSE) %>%
    dplyr::mutate(Date = gsub(".*(\\d{4}).csv$", "\\1", basename(file_path)))
}

# Standard curves + correction
create_standard_curves <- function(data) {
  gfp_models <- data %>%
    dplyr::filter(Condition %in% c("Alone", "Media Control")) %>%
    dplyr::group_by(RFP_Strain, Date) %>%
    dplyr::do(model = stats::lm(GFP ~ RFP, data = .)) %>%
    dplyr::filter(!is.null(model))
  
  rfp_models <- data %>%
    dplyr::filter(Condition %in% c("Alone", "Media Control")) %>%
    dplyr::group_by(GFP_Strain, Date) %>%
    dplyr::do(model = stats::lm(RFP ~ GFP, data = .)) %>%
    dplyr::filter(!is.null(model))
  
  list(
    gfp_models = gfp_models,
    rfp_models = rfp_models
  )
}

# Predict once per strain/date, rather than once per timepoint.
# The fitted lm objects and missing-model handling are unchanged.
estimate_autofluorescence <- function(data, gfp_models, rfp_models, media_means) {
  
  n <- nrow(data)
  estimate_gfp <- rep(NA_real_, n)
  estimate_rfp <- rep(NA_real_, n)
  
  rfp_strain <- as.character(data$RFP_Strain)
  gfp_strain <- as.character(data$GFP_Strain)
  date <- as.character(data$Date)
  
  # Media background for unlabelled controls.
  media_row <- match(date, as.character(media_means$Date))
  gfp_control <- rfp_strain %in% c("wt", "None")
  rfp_control <- gfp_strain %in% c("wt", "None")
  
  estimate_gfp[gfp_control] <-
    media_means$Media_GFP_Mean[media_row[gfp_control]]
  estimate_rfp[rfp_control] <-
    media_means$Media_RFP_Mean[media_row[rfp_control]]
  
  # Keep all observations in their original order. Split only the indices,
  # so each regression predicts its entire strain/date group in one call.
  key <- function(strain, day) paste(strain, day, sep = "\r")
  
  gfp_model_keys <- key(
    as.character(gfp_models$RFP_Strain),
    as.character(gfp_models$Date)
  )
  rfp_model_keys <- key(
    as.character(rfp_models$GFP_Strain),
    as.character(rfp_models$Date)
  )
  
  if (anyDuplicated(gfp_model_keys) || anyDuplicated(rfp_model_keys)) {
    stop("Duplicate strain/date models in standard_curves.")
  }
  
  gfp_groups <- split(
    which(!gfp_control),
    key(rfp_strain[!gfp_control], date[!gfp_control])
  )
  rfp_groups <- split(
    which(!rfp_control),
    key(gfp_strain[!rfp_control], date[!rfp_control])
  )
  
  for (group in names(gfp_groups)) {
    model_row <- match(group, gfp_model_keys)
    if (is.na(model_row)) next
    rows <- gfp_groups[[group]]
    estimate_gfp[rows] <- as.numeric(stats::predict(
      gfp_models$model[[model_row]],
      newdata = data.frame(RFP = data$RFP[rows])
    ))
  }
  
  for (group in names(rfp_groups)) {
    model_row <- match(group, rfp_model_keys)
    if (is.na(model_row)) next
    rows <- rfp_groups[[group]]
    estimate_rfp[rows] <- as.numeric(stats::predict(
      rfp_models$model[[model_row]],
      newdata = data.frame(GFP = data$GFP[rows])
    ))
  }
  
  data$Estimated_GFP_Autofluorescence <- estimate_gfp
  data$Estimated_RFP_Autofluorescence <- estimate_rfp
  data
}

correct_fluorescence <- function(data) {
  data %>%
    dplyr::mutate(
      Corrected_GFP = pmax(GFP - dplyr::coalesce(Estimated_GFP_Autofluorescence, 0), 0),
      Corrected_RFP = pmax(RFP - dplyr::coalesce(Estimated_RFP_Autofluorescence, 0), 0)
    ) %>%
    dplyr::select(
      -tidyselect::any_of(c(
        "Estimated_GFP_Autofluorescence",
        "Estimated_RFP_Autofluorescence"
      ))
    )
}

# Return NA for curves with fewer than two usable timepoints.
# For valid curves, use the same MESS::auc calculation as elsewhere.
auc_if_valid <- function(time, fluorescence) {
  usable <- is.finite(time) & is.finite(fluorescence)
  if (sum(usable) < 2L ||
      length(unique(time[usable])) < 2L) {
    return(NA_real_)
  }
  MESS::auc(time, fluorescence)
}

# AUC processing
#   - Ss uses Time <= 36 h
#   - all other competitors use Time <= 48 h
#   - RFP-side Ss values are normalized against a 36 h
#     same-day RFP-alone denominator
process_rfp_auc <- function(data,
                            rfp_strain,
                            T_END_DEFAULT = AUC_T_END_DEFAULT,
                            T_END_SS = AUC_T_END_SS) {
  
  data_processed <- data %>%
    dplyr::filter(RFP_Strain == rfp_strain) %>%
    dplyr::mutate(
      AUC_T_END = dplyr::if_else(
        as.character(GFP_Strain) == "Ss",
        T_END_SS,
        T_END_DEFAULT
      )
    ) %>%
    dplyr::filter(Time <= AUC_T_END) %>%
    dplyr::group_by(GFP_Strain, Date, Replicate, Time, AUC_T_END) %>%
    dplyr::summarise(
      Average_RFP = mean(Corrected_RFP, na.rm = TRUE),
      .groups = "drop"
    )
  
  auc_data <- data_processed %>%
    dplyr::group_by(GFP_Strain, Date, Replicate, AUC_T_END) %>%
    dplyr::summarise(
      AUC = auc_if_valid(Time, Average_RFP),
      .groups = "drop"
    )
  
  day_alone_auc <- purrr::map_dfr(
    c(T_END_DEFAULT, T_END_SS),
    function(t_end) {
      data %>%
        dplyr::filter(
          RFP_Strain == rfp_strain,
          GFP_Strain == "None",
          Time <= t_end
        ) %>%
        dplyr::group_by(Date, Replicate, Time) %>%
        dplyr::summarise(
          Average_RFP = mean(Corrected_RFP, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::group_by(Date, Replicate) %>%
        dplyr::summarise(
          AUC = auc_if_valid(Time, Average_RFP),
          .groups = "drop"
        ) %>%
        dplyr::group_by(Date) %>%
        dplyr::summarise(
          Day_Mean_AUC = mean(AUC, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::mutate(AUC_T_END = t_end)
    }
  )
  
  auc_data %>%
    dplyr::left_join(day_alone_auc, by = c("Date", "AUC_T_END")) %>%
    dplyr::mutate(
      Normalized_AUC = AUC / Day_Mean_AUC,
      RFP_Strain = rfp_strain,
      Species = stringr::str_sub(rfp_strain, 1, 2)
    ) %>%
    dplyr::select(
      -tidyselect::any_of(c("AUC", "Day_Mean_AUC", "AUC_T_END"))
    )
}

process_gfp_auc <- function(data,
                            gfp_strain,
                            T_END_DEFAULT = AUC_T_END_DEFAULT,
                            T_END_SS = AUC_T_END_SS) {
  
  T_END <- if (as.character(gfp_strain) == "Ss") {
    T_END_SS
  } else {
    T_END_DEFAULT
  }
  
  data_processed <- data %>%
    dplyr::filter(GFP_Strain == gfp_strain, Time <= T_END) %>%
    dplyr::group_by(Species, RFP_Strain, Date, Replicate, Time) %>%
    dplyr::summarise(
      Average_GFP = mean(Corrected_GFP, na.rm = TRUE),
      .groups = "drop"
    )
  
  auc_data <- data_processed %>%
    dplyr::group_by(Species, RFP_Strain, Date, Replicate) %>%
    dplyr::summarise(
      AUC = auc_if_valid(Time, Average_GFP),
      .groups = "drop"
    )
  
  day_alone_auc <- auc_data %>%
    dplyr::filter(RFP_Strain == "None") %>%
    dplyr::group_by(Date) %>%
    dplyr::summarise(
      Day_Mean_AUC = mean(AUC, na.rm = TRUE),
      .groups = "drop"
    )
  
  auc_data %>%
    dplyr::left_join(day_alone_auc, by = "Date") %>%
    dplyr::mutate(Normalized_AUC = AUC / Day_Mean_AUC) %>%
    dplyr::select(-tidyselect::any_of(c("AUC", "Day_Mean_AUC"))) %>%
    dplyr::filter(RFP_Strain != "None") %>%
    dplyr::mutate(GFP_Strain = gfp_strain)
}

# Load + merge plate reads and metadata
# GFP + RFP only here; OD is imported later
all_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

gfp_files <- grep("GFP", all_files, value = TRUE)
rfp_files <- grep("RFP", all_files, value = TRUE)
od_files <- grep("OD", all_files, value = TRUE)
metadata_files <- grep("Metadata", all_files, value = TRUE)

if (!length(gfp_files) || !length(rfp_files) ||
    !length(od_files) || !length(metadata_files)) {
  stop("Expected GFP, RFP, OD and Metadata CSV files in ", data_dir)
}

gfp_data <- purrr::map_df(
  gfp_files,
  read_measurement_data,
  measurement_type = "GFP"
)

rfp_data <- purrr::map_df(
  rfp_files,
  read_measurement_data,
  measurement_type = "RFP"
)

measurement_data <- purrr::reduce(
  list(gfp_data, rfp_data),
  dplyr::full_join,
  by = c("Time", "Well_ID", "Date")
)

rm(gfp_data, rfp_data)
gc()

metadata <- purrr::map_df(metadata_files, read_metadata)

full_data <- measurement_data %>%
  dplyr::left_join(metadata, by = c("Well_ID", "Date"))

rm(metadata, measurement_data)
gc()

# Clean timepoint-level data
full_data <- full_data %>%
  tidyr::drop_na(GFP, RFP, Condition)

# Fit standard curves, estimate autofluorescence, correct signals
standard_curves <- create_standard_curves(full_data)

media_control_means <- full_data %>%
  dplyr::filter(Condition == "Media Control") %>%
  dplyr::group_by(Date) %>%
  dplyr::summarise(
    Media_RFP_Mean = mean(RFP, na.rm = TRUE),
    Media_GFP_Mean = mean(GFP, na.rm = TRUE),
    .groups = "drop"
  )

full_data_corrected <- full_data %>%
  estimate_autofluorescence(
    gfp_models = standard_curves$gfp_models,
    rfp_models = standard_curves$rfp_models,
    media_means = media_control_means
  ) %>%
  correct_fluorescence()

rm(full_data, media_control_means)
gc()

# Add species label + enforce corrected Ss×wt control dates
full_data_corrected <- full_data_corrected %>%
  dplyr::mutate(Species = stringr::str_sub(RFP_Strain, 1, 2))

full_data_corrected <- full_data_corrected %>%
  dplyr::filter(
    !(GFP_Strain == "Ss" & RFP_Strain == "wt" & !(Date %in% ss_wt_dates))
  )

# Import OD separately, then join onto full_data_corrected
od_data <- purrr::map_df(
  od_files, read_measurement_data, measurement_type = "OD"
)

full_data_corrected <- full_data_corrected %>%
  dplyr::left_join(od_data, by = c("Time", "Well_ID", "Date")) %>%
  dplyr::mutate(
    OD = dplyr::if_else(is.na(OD), NA_real_, pmax(OD, OD_FLOOR))
  )

rm(od_data)
gc()

# AUC, RFP side
# Ss <= 36 h; all other competitors <= 48 h
rfp_results_list <- lapply(
  rfp_strains,
  function(strain) process_rfp_auc(full_data_corrected, strain)
)

merged_data <- dplyr::bind_rows(rfp_results_list)

merged_data$GFP_Strain <- factor(
  merged_data$GFP_Strain,
  levels = gfp_levels
)

# AUC, GFP side
# Ss <= 36 h; all other competitors <= 48 h
gfp_results_list <- lapply(
  gfp_strains,
  function(strain) process_gfp_auc(full_data_corrected, strain)
)

combined_gfp_data <- dplyr::bind_rows(gfp_results_list)

# Merge AUC datasets
merged_data <- merged_data %>%
  dplyr::rename(Normalized_AUC_RFP = Normalized_AUC)

combined_gfp_data <- combined_gfp_data %>%
  dplyr::rename(Normalized_AUC_GFP = Normalized_AUC)

combined_auc_data <- merged_data %>%
  dplyr::inner_join(
    combined_gfp_data,
    by = c("Species", "RFP_Strain", "Date", "Replicate", "GFP_Strain")
  )

# Add back the "wt-only" rows that fail the inner join
missing_wt_rows <- dplyr::bind_rows(
  combined_gfp_data %>%
    dplyr::filter(RFP_Strain == "wt") %>%
    dplyr::mutate(Normalized_AUC_RFP = NA_real_),
  merged_data %>%
    dplyr::filter(GFP_Strain == "wt") %>%
    dplyr::mutate(Normalized_AUC_GFP = NA_real_)
)

combined_auc_data <- dplyr::bind_rows(combined_auc_data, missing_wt_rows)

combined_auc_data <- combined_auc_data %>%
  dplyr::mutate(
    Reduction_Status = ifelse(
      Species %in% c("Bh", "Bb"),
      "Reduced",
      "Non-Reduced"
    )
  )

rm(combined_gfp_data, merged_data, missing_wt_rows)
gc()

# Exclude documented experimental failures from the analysis.
# Each exclusion applies to one RFP strain, GFP competitor and date.
# Valid additional experimental dates are retained.
# The 0805 date is excluded separately below.
auc_exclusions <- tibble::tribble(
  ~GFP_Strain, ~RFP_Strain, ~Date,  ~Reason,
  "Ck",        "Bb859",    "1606", "Abnormal GFP",
  "Pk",        "Ba161",    "1705", "Abnormal GFP",
  "Bp",        "Ba31",     "0206", "Abnormal GFP",
  "Bp",        "Bh11",     "0206", "Abnormal GFP",
  "Bp",        "Bh22",     "0610", "Abnormal GFP",
  "Pk",        "Bh22",     "1705", "Abnormal GFP",
  "Ck",        "Ba31",     "3009", "Abnormal GFP",
  "Ck",        "Bb859",    "2109", "Abnormal GFP",
  "Ck",        "Bh21",     "3009", "Abnormal GFP",
  "At",        "Bh21",     "1707", "Abnormal RFP",
  "Ss",        "Bh22",     "1705", "Abnormal RFP",
  "Ck",        "Bh21",     "1507", "Contamination"
)

exclusion_keys <- c("GFP_Strain", "RFP_Strain", "Date")
if (anyDuplicated(auc_exclusions[exclusion_keys])) {
  stop("Duplicate entries in auc_exclusions.")
}

cat("\nExcluded experimental runs by reason:\n")
print(dplyr::count(auc_exclusions, Reason))

filtered_data <- combined_auc_data %>%
  dplyr::filter(
    Date != "0805",
    !GFP_Strain %in% c("None", "Pkwt1", "Pkwt2", "Bh11", "Ba31")
  ) %>%
  dplyr::anti_join(auc_exclusions, by = c("GFP_Strain", "RFP_Strain", "Date"))

filtered_data <- filtered_data %>%
  dplyr::mutate(
    RFP_Group = dplyr::case_when(
      stringr::str_detect(RFP_Strain, "^Ba") ~ "Ba",
      stringr::str_detect(RFP_Strain, "^Bb") ~ "Bb",
      stringr::str_detect(RFP_Strain, "^Bh") ~ "Bh",
      RFP_Strain %in% c("Pf", "Pu", "Ps") ~ "Other",
      RFP_Strain %in% c("wt", "None") ~ RFP_Strain,
      TRUE ~ "Other"
    ),
    Species = dplyr::if_else(
      RFP_Strain %in% c("wt", "None"),
      NA_character_,
      as.character(Species)
    )
  )

# Bridge to modeling scripts
# Builds:
#   - analysis_data_rep = replicate-level log-ratio data
#   - analysis_data_day = day-mean data, technical replicates averaged
gfp_model_levels <- c("Sl", "Ss", "At", "Bp", "Ck", "Pk")
rfp_species6_levels <- c("Pf", "Pu", "Ps", "Ba", "Bb", "Bh")

analysis_data_rep <- filtered_data %>%
  dplyr::mutate(
    RFP_Species6 = dplyr::case_when(
      stringr::str_detect(RFP_Strain, "^Ba") ~ "Ba",
      stringr::str_detect(RFP_Strain, "^Bb") ~ "Bb",
      stringr::str_detect(RFP_Strain, "^Bh") ~ "Bh",
      RFP_Strain == "Pf" ~ "Pf",
      RFP_Strain == "Pu" ~ "Pu",
      RFP_Strain == "Ps" ~ "Ps",
      TRUE ~ NA_character_
    ),
    Ratio_RFP_GFP = Normalized_AUC_RFP / Normalized_AUC_GFP,
    log_ratio = log10(Ratio_RFP_GFP)
  ) %>%
  dplyr::filter(
    !is.na(RFP_Species6),
    GFP_Strain %in% gfp_model_levels,
    is.finite(Normalized_AUC_RFP),
    is.finite(Normalized_AUC_GFP),
    Normalized_AUC_RFP > 0,
    Normalized_AUC_GFP > 0,
    is.finite(log_ratio)
  ) %>%
  dplyr::mutate(
    GFP_Strain = factor(GFP_Strain, levels = gfp_model_levels),
    RFP_Species6 = factor(RFP_Species6, levels = rfp_species6_levels),
    RFP_Strain = factor(RFP_Strain),
    Date = factor(Date)
  )

analysis_data_day <- analysis_data_rep %>%
  dplyr::group_by(Date, RFP_Strain, RFP_Species6, GFP_Strain) %>%
  dplyr::summarise(
    sd_log_ratio = sd(log_ratio, na.rm = TRUE),
    n_tech = dplyr::n(),
    log_ratio = mean(log_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    GFP_Strain = factor(GFP_Strain, levels = gfp_model_levels),
    RFP_Species6 = factor(RFP_Species6, levels = rfp_species6_levels),
    RFP_Strain = factor(RFP_Strain),
    Date = factor(Date)
  )

# Check the analysis input.
if (nrow(analysis_data_day) != 282L ||
    dplyr::n_distinct(analysis_data_day$RFP_Strain) != 15L) {
  stop("Expected 282 date-level observations from 15 RFP strains.")
}

if (any(is.na(analysis_data_day$sd_log_ratio) &
        analysis_data_day$n_tech > 1L)) {
  stop("Missing SD for a date with multiple usable technical replicates.")
}

if (anyDuplicated(analysis_data_day[c(
  "Date", "RFP_Strain", "GFP_Strain"
)])) {
  stop("Duplicate experimental dates in analysis_data_day.")
}

replicate_counts <- analysis_data_day %>%
  dplyr::count(RFP_Strain, GFP_Strain, name = "Num_Dates")

if (nrow(replicate_counts) != 90L ||
    any(replicate_counts$Num_Dates < 3L |
        replicate_counts$Num_Dates > 5L)) {
  stop("Expected three to five dates for each of 90 strain/competitor pairs.")
}

if (nrow(dplyr::semi_join(
  analysis_data_day,
  auc_exclusions,
  by = exclusion_keys
)) > 0L) {
  stop("Excluded experimental runs remain in analysis_data_day.")
}

count_dates <- replicate_counts %>%
  tidyr::pivot_wider(
    names_from = GFP_Strain,
    values_from = Num_Dates,
    values_fill = 0
  )

cat("\nDate counts by RFP strain and competitor:\n")
print(count_dates, n = Inf)

if (SAVE_DATE_COUNTS_XLSX) {
  if (!requireNamespace("writexl", quietly = TRUE)) {
    stop("Install 'writexl' to export replicate counts.")
  }
  writexl::write_xlsx(count_dates, "ReplicateNumber_all_valid.xlsx")
}

print(summary(analysis_data_rep$log_ratio))
print(table(analysis_data_rep$RFP_Species6))
print(table(analysis_data_rep$GFP_Strain))
print(head(analysis_data_day))

saveRDS(analysis_data_rep, "analysis_data_rep_all_valid.rds")
saveRDS(analysis_data_day, "analysis_data_day_all_valid.rds")


cat("\nSaved model inputs:\n")
cat("  analysis_data_rep_all_valid.rds\n")
cat("  analysis_data_day_all_valid.rds\n")
