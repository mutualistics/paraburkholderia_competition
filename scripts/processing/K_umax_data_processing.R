
# K and maximum early growth-rate calculations
# Run after the final 1_import_comp_data.R in the same R session.
# Uses 282 accepted competition dates from analysis_data_day_all_valid.rds.
# The existing trait estimators, QC, and WT normalization are unchanged.
# No files are written. K: 0-36 h; mu_max: 0-24 h.

suppressPackageStartupMessages({
  library(tidyverse)
  library(growthcurver)
  library(zoo)
})

THR <- list(
  loose_sig        = 0.50,
  loose_snr        = 2.0,
  strict_sig       = 0.30,
  strict_snr       = 3.0,
  loess_span       = 0.20,
  neg_slope_thr    = -0.001,
  neg_run_len      = 3,
  K_MIN            = 50,
  fluor_floor      = 1,
  umax_rate_min    = 1e-6,
  ratio_floor      = 1e-12,
  umax_window_points   = 5,
  umax_baseline_t_max  = 2,
  umax_noise_mult      = 2,
  umax_min_above_frac  = 0.60
)

stopifnot(
  THR$K_MIN > 0,
  THR$fluor_floor > 0,
  THR$umax_rate_min > 0,
  THR$ratio_floor > 0
)

T_END_K <- 36
T_END_UMAX <- 24

# Helpers

print_tbl <- function(x, n = Inf) {
  x <- as.data.frame(x)
  if (is.finite(n) && nrow(x) > n) {
    x <- x[seq_len(n), , drop = FALSE]
  }
  print(x, row.names = FALSE)
  invisible(x)
}

get_num1 <- function(x, name, default = NA_real_) {
  if (is.null(x)) return(default)
  val <- tryCatch(x[[name]], error = function(e) NULL)
  if (is.null(val) || length(val) == 0) default else as.numeric(val[1])
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  median(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x)
}

safe_mad <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(0)
  
  out <- stats::mad(
    x,
    center = median(x),
    constant = 1.4826,
    na.rm = TRUE
  )
  
  if (!is.finite(out)) 0 else out
}

fit_one_gc_dbg <- function(t, y) {
  tryCatch({
    t <- as.numeric(t)
    y <- as.numeric(y)
    
    ok <- is.finite(t) & is.finite(y)
    t <- t[ok]
    y <- y[ok]
    
    ord <- order(t)
    t <- t[ord]
    y <- y[ord]
    
    t <- t - min(t, na.rm = TRUE)
    y <- pmax(y, THR$fluor_floor)
    
    y_max <- max(y, na.rm = TRUE)
    y_s <- y / y_max
    
    out <- growthcurver::SummarizeGrowth(
      data_t = t,
      data_n = y_s
    )
    
    list(
      ok = TRUE,
      err = NA_character_,
      vals = out$vals,
      y_max = y_max
    )
  }, error = function(e) {
    list(
      ok = FALSE,
      err = conditionMessage(e),
      vals = NULL,
      y_max = NA_real_
    )
  })
}

get_umax_signal_threshold <- function(t,
                                      y,
                                      baseline_t_max = THR$umax_baseline_t_max,
                                      noise_mult = THR$umax_noise_mult,
                                      eps = THR$fluor_floor) {
  
  base_idx <- is.finite(t) &
    is.finite(y) &
    t <= baseline_t_max
  
  if (sum(base_idx) < 3) {
    base_idx <- seq_len(min(3, length(y)))
  }
  
  baseline_signal <- median(y[base_idx], na.rm = TRUE)
  baseline_noise <- safe_mad(y[base_idx])
  
  if (!is.finite(baseline_noise) || baseline_noise == 0) {
    baseline_noise <- median(abs(diff(y[base_idx])), na.rm = TRUE)
  }
  
  if (!is.finite(baseline_noise)) {
    baseline_noise <- 0
  }
  
  signal_threshold <- max(
    eps,
    baseline_signal + noise_mult * baseline_noise,
    na.rm = TRUE
  )
  
  tibble::tibble(
    baseline_signal = baseline_signal,
    baseline_noise = baseline_noise,
    signal_threshold = signal_threshold
  )
}

calc_umax_rolling_points_0_24 <- function(t,
                                          y,
                                          window_points = THR$umax_window_points,
                                          baseline_t_max = THR$umax_baseline_t_max,
                                          noise_mult = THR$umax_noise_mult,
                                          min_above_frac = THR$umax_min_above_frac,
                                          eps = THR$fluor_floor) {
  
  empty_candidate_tbl <- function() {
    tibble::tibble(
      slope = numeric(),
      t_mid = numeric(),
      above_frac = numeric(),
      window_min_t = numeric(),
      window_max_t = numeric()
    )
  }
  
  t <- as.numeric(t)
  y <- pmax(as.numeric(y), eps)
  
  ok <- is.finite(t) &
    is.finite(y) &
    t >= 0 &
    t <= T_END_UMAX
  
  t <- t[ok]
  y <- y[ok]
  
  if (length(t) < window_points) {
    return(
      tibble::tibble(
        umax = NA_real_,
        t_umax = NA_real_,
        n_points = length(t),
        n_candidates = 0L,
        n_positive_slopes = NA_integer_,
        max_signal = ifelse(length(y) > 0, max(y, na.rm = TRUE), NA_real_),
        baseline_signal = NA_real_,
        baseline_noise = NA_real_,
        signal_threshold = NA_real_
      )
    )
  }
  
  o <- order(t)
  t <- t[o]
  y <- y[o]
  log_y <- log(y)
  
  threshold_tbl <- get_umax_signal_threshold(
    t = t,
    y = y,
    baseline_t_max = baseline_t_max,
    noise_mult = noise_mult,
    eps = eps
  )
  
  signal_threshold <- threshold_tbl$signal_threshold[1]
  starts <- seq_len(length(t) - window_points + 1)
  
  candidate_tbl <- purrr::map_dfr(
    starts,
    function(i) {
      
      idx <- i:(i + window_points - 1)
      
      t_win <- t[idx]
      y_win <- y[idx]
      log_win <- log_y[idx]
      
      above_frac <- mean(y_win >= signal_threshold, na.rm = TRUE)
      
      if (above_frac < min_above_frac) {
        return(empty_candidate_tbl())
      }
      
      slope <- tryCatch(
        stats::coef(stats::lm(log_win ~ t_win))[2],
        error = function(e) NA_real_
      )
      
      tibble::tibble(
        slope = as.numeric(slope),
        t_mid = median(t_win),
        above_frac = above_frac,
        window_min_t = min(t_win),
        window_max_t = max(t_win)
      )
    }
  ) %>%
    dplyr::filter(
      is.finite(slope),
      is.finite(t_mid)
    )
  
  if (nrow(candidate_tbl) == 0) {
    return(
      tibble::tibble(
        umax = NA_real_,
        t_umax = NA_real_,
        n_points = length(t),
        n_candidates = 0L,
        n_positive_slopes = NA_integer_,
        max_signal = max(y, na.rm = TRUE),
        baseline_signal = threshold_tbl$baseline_signal[1],
        baseline_noise = threshold_tbl$baseline_noise[1],
        signal_threshold = signal_threshold
      )
    )
  }
  
  best <- candidate_tbl %>%
    dplyr::slice_max(slope, n = 1, with_ties = FALSE)
  
  tibble::tibble(
    umax = best$slope,
    t_umax = best$t_mid,
    n_points = length(t),
    n_candidates = nrow(candidate_tbl),
    n_positive_slopes = sum(candidate_tbl$slope > 0, na.rm = TRUE),
    max_signal = max(y, na.rm = TRUE),
    baseline_signal = threshold_tbl$baseline_signal[1],
    baseline_noise = threshold_tbl$baseline_noise[1],
    signal_threshold = signal_threshold
  )
}

still_increasing_end <- function(t, y, tail_frac = 0.15, slope_thr = 0.002) {
  t <- as.numeric(t)
  y <- pmax(as.numeric(y), THR$fluor_floor)
  
  ok <- is.finite(t) & is.finite(y)
  t <- t[ok]
  y <- y[ok]
  
  if (length(t) < 8) return(FALSE)
  
  o <- order(t)
  t <- t[o]
  y <- y[o]
  
  n <- length(t)
  k <- max(3L, floor(tail_frac * n))
  i <- seq.int(n - k + 1L, n)
  
  s <- tryCatch(
    stats::coef(stats::lm(log(y[i]) ~ t[i]))[2],
    error = function(e) NA_real_
  )
  
  is.finite(s) && s > slope_thr
}

K_visual_top <- function(t, y, q = 0.98, span = THR$loess_span) {
  t <- as.numeric(t)
  y <- pmax(as.numeric(y), THR$fluor_floor)
  
  ok <- is.finite(t) & is.finite(y)
  t <- t[ok]
  y <- y[ok]
  
  if (length(t) < 6) return(NA_real_)
  
  o <- order(t)
  t <- t[o]
  y <- y[o]
  
  ys <- try(
    stats::predict(stats::loess(y ~ t, span = span)),
    silent = TRUE
  )
  
  ys <- if (inherits(ys, "try-error") || any(!is.finite(ys))) {
    y
  } else {
    ys
  }
  
  as.numeric(stats::quantile(ys, probs = q, na.rm = TRUE))
}

K_star <- function(t,
                   y,
                   span = THR$loess_span,
                   neg_run_len = THR$neg_run_len,
                   neg_slope_thr = THR$neg_slope_thr) {
  t <- as.numeric(t)
  y <- as.numeric(y)
  
  ok <- is.finite(t) & is.finite(y)
  t <- t[ok]
  y <- pmax(y[ok], THR$fluor_floor)
  
  if (length(t) < 4) return(NA_real_)
  
  o <- order(t)
  t <- t[o]
  y <- y[o]
  
  ys <- try(
    stats::predict(stats::loess(y ~ t, span = span)),
    silent = TRUE
  )
  
  if (inherits(ys, "try-error") || any(!is.finite(ys))) {
    ys <- y
  }
  
  dy <- c(NA, diff(ys)) / c(NA, diff(t))
  pk <- which.max(ys)
  n <- length(ys)
  
  if (pk < n) {
    post <- dy[(pk + 1):n]
    rr <- rle(post < neg_slope_thr)
    
    if (any(rr$values & rr$lengths >= neg_run_len)) {
      first_run <- which(rr$values & rr$lengths >= neg_run_len)[1]
      cut_idx <- pk + sum(rr$lengths[seq_len(first_run - 1)]) + 1
      return(max(y[seq_len(max(1, cut_idx - 1))], na.rm = TRUE))
    }
  }
  
  tail_n <- max(3, ceiling(0.10 * n))
  mean(tail(y, tail_n), na.rm = TRUE)
}

apply_exclusions <- function(df, excl) {
  df <- df %>%
    dplyr::mutate(
      dplyr::across(c(Date, GFP_Strain, RFP_Strain), as.character)
    )
  
  excl <- excl %>%
    dplyr::mutate(
      dplyr::across(c(Date, GFP_Strain, RFP_Strain), as.character)
    )
  
  drop_days <- excl %>%
    dplyr::filter(is.na(GFP_Strain) & is.na(RFP_Strain)) %>%
    dplyr::distinct(Date)
  
  if (nrow(drop_days)) {
    df <- df %>%
      dplyr::filter(!(Date %in% drop_days$Date))
  }
  
  gfp_only <- excl %>%
    dplyr::filter(!is.na(GFP_Strain) & is.na(RFP_Strain)) %>%
    dplyr::distinct(Date, GFP_Strain)
  
  if (nrow(gfp_only)) {
    df <- df %>%
      dplyr::anti_join(gfp_only, by = c("Date", "GFP_Strain"))
  }
  
  rfp_only <- excl %>%
    dplyr::filter(is.na(GFP_Strain) & !is.na(RFP_Strain)) %>%
    dplyr::distinct(Date, RFP_Strain)
  
  if (nrow(rfp_only)) {
    df <- df %>%
      dplyr::anti_join(rfp_only, by = c("Date", "RFP_Strain"))
  }
  
  triples <- excl %>%
    dplyr::filter(!is.na(GFP_Strain) & !is.na(RFP_Strain)) %>%
    dplyr::distinct(Date, GFP_Strain, RFP_Strain)
  
  if (nrow(triples)) {
    df <- df %>%
      dplyr::anti_join(triples, by = c("Date", "GFP_Strain", "RFP_Strain"))
  }
  
  df
}

fix_condition <- function(rfp, gfp) {
  dplyr::case_when(
    rfp %in% c("None", "YG") & gfp %in% c("None", "YG") ~ "Media Control",
    xor(rfp %in% c("None", "YG"), gfp %in% c("None", "YG")) ~ "Alone",
    rfp %in% c("wt", "wt2") | gfp %in% c("wt", "wt2") ~ "Competing_wt",
    TRUE ~ "Competing"
  )
}

mad_z <- function(x) {
  m <- median(x, na.rm = TRUE)
  md <- stats::mad(x, center = m, constant = 1.4826, na.rm = TRUE)
  (x - m) / pmax(md, 1e-9)
}

is_wt_pair <- function(rfp, gfp) {
  rfp %in% c("wt", "wt2") | gfp %in% c("wt", "wt2")
}

# Corrected fluorescence from the data import

if (!exists("full_data_corrected", inherits = TRUE)) {
  stop("Run the competition data import first (full_data_corrected is missing).")
}

need_cols <- c(
  "Time",
  "Well_ID",
  "Date",
  "Replicate",
  "Condition",
  "RFP_Strain",
  "GFP_Strain",
  "Corrected_RFP",
  "Corrected_GFP"
)

miss <- setdiff(need_cols, names(full_data_corrected))

if (length(miss)) {
  stop(
    "full_data_corrected is missing required columns: ",
    paste(miss, collapse = ", ")
  )
}

# 1) Experimental-date selection
# Use the same accepted strain x competitor x date combinations as the
# AUC-ratio analysis. Trait-specific QC remains independent below.
# Monoculture and WT competition controls are not restricted to these keys.

INPUT_CANDIDATES <- c(
  "analysis_data_day_all_valid.rds",
  file.path("data", "processed", "analysis_data_day_all_valid.rds")
)

input_found <- INPUT_CANDIDATES[file.exists(INPUT_CANDIDATES)]
if (!length(input_found)) {
  stop("Missing analysis_data_day_all_valid.rds. Run the final AUC import first.")
}

accepted_dates <- readRDS(input_found[[1]])
required_keys <- c("Date", "RFP_Strain", "GFP_Strain")

if (!all(required_keys %in% names(accepted_dates))) {
  stop("The final AUC dataset is missing accepted-date identifiers.")
}

accepted_dates <- accepted_dates %>%
  dplyr::mutate(dplyr::across(dplyr::all_of(required_keys), as.character)) %>%
  dplyr::distinct(dplyr::across(dplyr::all_of(required_keys)))

if (nrow(accepted_dates) != 282L ||
    dplyr::n_distinct(accepted_dates$RFP_Strain) != 15L ||
    nrow(dplyr::distinct(accepted_dates, RFP_Strain, GFP_Strain)) != 90L ||
    anyNA(accepted_dates)) {
  stop("Expected 282 accepted dates across 90 strain/competitor pairs and 15 RFP strains.")
}

# These controls are excluded independently of the competition-date selection.
# The 0805 plate is excluded globally, as in the AUC import.

control_exclusions <- tibble::tribble(
  ~Date, ~GFP_Strain, ~RFP_Strain,
  "0805", NA_character_, NA_character_,
  "1305", "Pk", "wt2",
  "1305", "Pkwt1", NA_character_,
  "1305", "Pkwt2", NA_character_
)

ignore_gfp_strains <- c(
  "At3",
  "Pkwt1",
  "Pkwt2",
  "Bh11",
  "Ba31",
  "Ss_wt",
  "Ps_wt"
)

ignore_rfp_strains <- c(
  "wt2",
  "Pf_wt",
  "Pu_wt",
  "Ps_wt"
)

filtered_curves <- full_data_corrected %>%
  dplyr::mutate(
    Date = as.character(Date),
    RFP_Strain = as.character(RFP_Strain),
    GFP_Strain = as.character(GFP_Strain),
    RFP_Strain = dplyr::if_else(RFP_Strain == "YG", "None", RFP_Strain),
    GFP_Strain = dplyr::if_else(GFP_Strain == "YG", "None", GFP_Strain)
  ) %>%
  dplyr::filter(!(GFP_Strain %in% ignore_gfp_strains)) %>%
  dplyr::filter(!(RFP_Strain %in% ignore_rfp_strains)) %>%
  apply_exclusions(control_exclusions) %>%
  dplyr::mutate(
    Condition_fixed = fix_condition(RFP_Strain, GFP_Strain)
  )

# Restrict only competition curves to accepted experimental dates.
# Keep the original monoculture and WT control curves for normalization.

competition_curves <- filtered_curves %>%
  dplyr::filter(Condition_fixed == "Competing") %>%
  dplyr::semi_join(accepted_dates, by = required_keys)

retained_dates <- competition_curves %>%
  dplyr::distinct(dplyr::across(dplyr::all_of(required_keys)))

missing_dates <- dplyr::anti_join(
  accepted_dates, retained_dates, by = required_keys
)

if (nrow(missing_dates) > 0L) {
  cat("\nAccepted dates missing from the growth-curve data:\n")
  print_tbl(missing_dates)
  stop("Cannot calculate K and mu_max for all 282 accepted experimental dates.")
}

filtered_curves <- dplyr::bind_rows(
  filtered_curves %>% dplyr::filter(Condition_fixed != "Competing"),
  competition_curves
)

cat("\n================ COMPETITION DATE CHECK ================\n")
cat("Accepted strain/competitor/date combinations:", nrow(retained_dates), "\n")
cat("RFP strains:", dplyr::n_distinct(retained_dates$RFP_Strain), "\n")
cat("Strain/competitor pairs:",
    nrow(dplyr::distinct(retained_dates, RFP_Strain, GFP_Strain)), "\n")

cat("\n================ CONDITION MISMATCH CHECK ================\n")

filtered_curves %>%
  dplyr::mutate(mislabel = Condition != Condition_fixed) %>%
  dplyr::filter(mislabel) %>%
  dplyr::count(
    Condition,
    Condition_fixed,
    RFP_Strain,
    GFP_Strain,
    Date,
    sort = TRUE
  ) %>%
  print_tbl(n = 50)

# 2) Build long tables

make_curves_long <- function(dat, t_end) {
  dat %>%
    dplyr::filter(Time <= t_end) %>%
    dplyr::filter(Condition_fixed != "Media Control") %>%
    dplyr::select(
      Time,
      Well_ID,
      Date,
      Condition = Condition_fixed,
      Replicate,
      RFP_Strain,
      GFP_Strain,
      RFP = Corrected_RFP,
      GFP = Corrected_GFP
    ) %>%
    tidyr::pivot_longer(
      cols = c(RFP, GFP),
      names_to = "Channel",
      values_to = "Signal"
    ) %>%
    dplyr::arrange(Date, Well_ID, Channel, Time) %>%
    dplyr::group_by(Date, Well_ID, Channel) %>%
    dplyr::mutate(
      Signal = pmax(as.numeric(Signal), THR$fluor_floor)
    ) %>%
    dplyr::ungroup()
}

curves_K <- make_curves_long(
  filtered_curves,
  t_end = T_END_K
)

curves_umax <- make_curves_long(
  filtered_curves,
  t_end = T_END_UMAX
)

curves_K_nonempty <- curves_K %>%
  dplyr::group_by(Date, Well_ID, Channel) %>%
  dplyr::filter(any(is.finite(Time) & is.finite(Signal))) %>%
  dplyr::ungroup()

curves_umax_nonempty <- curves_umax %>%
  dplyr::group_by(Date, Well_ID, Channel) %>%
  dplyr::filter(any(is.finite(Time) & is.finite(Signal))) %>%
  dplyr::ungroup()

cat("\n================ TIME WINDOW CHECK ================\n")

cat("\nK/K* curves max time:\n")
curves_K_nonempty %>%
  dplyr::summarise(max_time = max(Time, na.rm = TRUE)) %>%
  print_tbl()

cat("\numax curves max time:\n")
curves_umax_nonempty %>%
  dplyr::summarise(max_time = max(Time, na.rm = TRUE)) %>%
  print_tbl()

# 3) Fit K with Growthcurver and μmax by selected 0-24 h estimator

dbg <- curves_K_nonempty %>%
  dplyr::group_by(Date, Well_ID, Channel) %>%
  dplyr::summarise(
    res = list(fit_one_gc_dbg(Time, Signal)),
    .groups = "drop"
  )

fits_core <- dbg %>%
  dplyr::mutate(
    gc_ok = purrr::map_lgl(res, ~ isTRUE(.x$ok)),
    K = purrr::map_dbl(
      res,
      ~ if (isTRUE(.x$ok)) get_num1(.x$vals, "k") * .x$y_max else NA_real_
    ),
    sigma = purrr::map_dbl(
      res,
      ~ if (isTRUE(.x$ok)) get_num1(.x$vals, "sigma") else NA_real_
    )
  ) %>%
  dplyr::select(-res)

umax_core <- curves_umax_nonempty %>%
  dplyr::group_by(Date, Well_ID, Channel) %>%
  dplyr::summarise(
    calc = list(calc_umax_rolling_points_0_24(Time, Signal)),
    .groups = "drop"
  ) %>%
  tidyr::unnest_wider(calc)

fits_core <- fits_core %>%
  dplyr::left_join(
    umax_core,
    by = c("Date", "Well_ID", "Channel")
  )

cat("\n================ CURVE FIT SUMMARY ================\n")

fits_core %>%
  dplyr::summarise(
    n = dplyr::n(),
    gc_ok = sum(gc_ok),
    gc_ok_pct = 100 * mean(gc_ok),
    finite_K = sum(is.finite(K)),
    finite_umax = sum(is.finite(umax)),
    finite_umax_pct = 100 * mean(is.finite(umax)),
    median_t_umax = median(t_umax, na.rm = TRUE),
    median_umax_candidates = median(n_candidates, na.rm = TRUE)
  ) %>%
  print_tbl()

# 4) QC metrics for K/K* curves, 0-36 h

qc_tbl <- curves_K_nonempty %>%
  dplyr::group_by(Date, Well_ID, Channel) %>%
  dplyr::arrange(Time, .by_group = TRUE) %>%
  dplyr::summarise(
    dyn_range = diff(range(Signal, na.rm = TRUE)),
    noise = median(abs(diff(Signal)), na.rm = TRUE),
    snr = dyn_range / pmax(noise, 1e-6),
    death_flag = {
      ys <- try(
        stats::predict(stats::loess(Signal ~ Time, span = THR$loess_span)),
        silent = TRUE
      )
      ys <- if (inherits(ys, "try-error") || any(!is.finite(ys))) Signal else ys
      dy <- c(NA, diff(ys)) / c(NA, diff(Time))
      n <- length(dy)
      pk <- which.max(ys)
      post <- if (pk < n) dy[(pk + 1):n] else numeric(0)
      rr <- if (length(post) > 0) {
        rle(post < THR$neg_slope_thr)
      } else {
        list(values = FALSE, lengths = 0)
      }
      any(rr$values & rr$lengths >= THR$neg_run_len)
    },
    .groups = "drop"
  )

meta <- curves_K_nonempty %>%
  dplyr::distinct(
    Date,
    Well_ID,
    Channel,
    Condition,
    Replicate,
    RFP_Strain,
    GFP_Strain
  )

fits <- fits_core %>%
  dplyr::left_join(qc_tbl, by = c("Date", "Well_ID", "Channel")) %>%
  dplyr::left_join(meta, by = c("Date", "Well_ID", "Channel")) %>%
  dplyr::mutate(
    sigma_norm = sigma / pmax(dyn_range, 1e-6),
    logistic_ok_loose = gc_ok & is.finite(K) &
      (sigma_norm <= THR$loose_sig) &
      (snr >= THR$loose_snr),
    logistic_ok_strict = gc_ok & is.finite(K) &
      (sigma_norm <= THR$strict_sig) &
      (snr >= THR$strict_snr) &
      !death_flag
  )

cat("\n================ LOGISTIC QC SUMMARY ================\n")

fits %>%
  dplyr::summarise(
    strict = mean(logistic_ok_strict, na.rm = TRUE),
    loose = mean(logistic_ok_loose, na.rm = TRUE),
    p_death = mean(death_flag, na.rm = TRUE)
  ) %>%
  print_tbl()

fits %>%
  dplyr::summarise(
    dplyr::across(
      c(sigma_norm, snr),
      list(
        med = ~ median(.x, na.rm = TRUE),
        q10 = ~ quantile(.x, 0.10, na.rm = TRUE),
        q90 = ~ quantile(.x, 0.90, na.rm = TRUE)
      )
    )
  ) %>%
  print_tbl()

fits %>%
  dplyr::mutate(
    reason = dplyr::case_when(
      logistic_ok_strict ~ "strict_ok",
      gc_ok & !death_flag & (sigma_norm <= THR$strict_sig) & (snr < THR$strict_snr) ~ "low_snr",
      gc_ok & !death_flag & (sigma_norm > THR$strict_sig) ~ "high_sigma",
      gc_ok & death_flag ~ "death_flag",
      TRUE ~ "other"
    )
  ) %>%
  dplyr::count(reason) %>%
  print_tbl()

# 5) K fallback and final parameters

fallbacks <- curves_K_nonempty %>%
  dplyr::group_by(Date, Well_ID, Channel) %>%
  dplyr::summarise(
    K_fb = K_star(Time, Signal),
    K_vis = K_visual_top(Time, Signal, q = 0.98),
    inc_flag = still_increasing_end(Time, Signal),
    .groups = "drop"
  )

fits_final <- fits %>%
  dplyr::left_join(fallbacks, by = c("Date", "Well_ID", "Channel")) %>%
  dplyr::mutate(
    K_unconstrained_logistic = K,
    K_final = dplyr::case_when(
      logistic_ok_strict & is.finite(K_vis) ~ K_vis,
      is.finite(K_fb)                       ~ K_fb,
      TRUE                                  ~ NA_real_
    ),
    umax_final = umax,
    K_method = dplyr::case_when(
      logistic_ok_strict & is.finite(K_vis) ~ "logistic_supported_Kvis",
      is.finite(K_fb)                       ~ "K_star",
      TRUE                                  ~ "missing_K"
    ),
    umax_method = "rolling_5pt_log_slope_0_24h_baseline_plus_2MAD_3of5"
  )

cat("\n================ FINAL PARAMETER SUMMARY ================\n")

fits_final %>%
  dplyr::summarise(
    n = dplyr::n(),
    finite_K_final = sum(is.finite(K_final)),
    finite_umax_final = sum(is.finite(umax_final)),
    finite_umax_pct = 100 * mean(is.finite(umax_final)),
    logistic_supported_Kvis = sum(K_method == "logistic_supported_Kvis", na.rm = TRUE),
    K_star = sum(K_method == "K_star", na.rm = TRUE),
    missing_K = sum(K_method == "missing_K", na.rm = TRUE),
    unconstrained_logistic_used = sum(K_method == "logistic", na.rm = TRUE),
    logistic_supported_Kvis_pct = 100 * logistic_supported_Kvis / n,
    K_star_pct = 100 * K_star / n
  ) %>%
  print_tbl()

stopifnot(
  sum(fits_final$K_method == "logistic", na.rm = TRUE) == 0
)

cat("\n================ UMAX ROLE CHECK ================\n")

fits_final %>%
  dplyr::mutate(
    failed_umax = !is.finite(umax_final),
    curve_role = dplyr::case_when(
      Condition == "Alone" &
        Channel == "RFP" &
        GFP_Strain == "None" &
        RFP_Strain != "None" ~ "RFP monoculture baseline",
      
      Condition == "Alone" &
        Channel == "GFP" &
        RFP_Strain == "None" &
        GFP_Strain != "None" ~ "GFP monoculture baseline",
      
      Condition == "Alone" ~ "non-focal alone channel",
      Condition == "Competing" ~ "competition curve",
      Condition == "Competing_wt" ~ "WT correction curve",
      TRUE ~ "other"
    )
  ) %>%
  dplyr::count(curve_role, failed_umax) %>%
  dplyr::group_by(curve_role) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  dplyr::ungroup() %>%
  print_tbl(n = Inf)

# 6) Aggregate technical replicates to day level

fits_day <- fits_final %>%
  dplyr::group_by(Date, Channel, Condition, RFP_Strain, GFP_Strain) %>%
  dplyr::summarise(
    K_day = safe_median(K_final),
    umax_day = safe_median(umax_final),
    K_sd = safe_sd(K_final),
    umax_sd = safe_sd(umax_final),
    t_umax_day = safe_median(t_umax),
    n_tech = dplyr::n(),
    n_umax_finite_tech = sum(is.finite(umax_final)),
    .groups = "drop"
  )

fits_day_clean <- fits_day %>%
  dplyr::filter(
    (RFP_Strain != "wt" & GFP_Strain != "wt") |
      (GFP_Strain == "wt" & Channel == "RFP") |
      (RFP_Strain == "wt" & Channel == "GFP")
  )

# 7) Same-day monoculture baselines

mono_day <- fits_day_clean %>%
  dplyr::filter(
    Condition == "Alone",
    (Channel == "RFP" & GFP_Strain == "None") |
      (Channel == "GFP" & RFP_Strain == "None")
  ) %>%
  dplyr::mutate(
    Focal_Strain = ifelse(Channel == "RFP", RFP_Strain, GFP_Strain)
  ) %>%
  dplyr::select(
    Date,
    Channel,
    Focal_Strain,
    K_mono = K_day,
    umax_mono = umax_day
  )

fits_day_tag <- fits_day_clean %>%
  dplyr::mutate(is_wt = is_wt_pair(RFP_Strain, GFP_Strain))

fits_norm <- fits_day_tag %>%
  dplyr::mutate(
    Focal_Strain = ifelse(Channel == "RFP", RFP_Strain, GFP_Strain)
  ) %>%
  dplyr::left_join(
    mono_day,
    by = c("Date", "Channel", "Focal_Strain")
  ) %>%
  dplyr::mutate(
    K_ratio = K_day / pmax(K_mono, THR$K_MIN),
    umax_ratio = pmax(umax_day, THR$umax_rate_min) /
      pmax(umax_mono, THR$umax_rate_min),
    dK_pct = 100 * (K_ratio - 1),
    dumax_pct = 100 * (umax_ratio - 1),
    log2K = log2(pmax(K_ratio, THR$ratio_floor)),
    log2umax = log2(pmax(umax_ratio, THR$ratio_floor))
  )

cat("\n================ BASELINE CHECK ================\n")

fits_norm %>%
  dplyr::mutate(
    K_floor = !is.na(K_mono) & K_mono < THR$K_MIN,
    umax_rate_floor = !is.na(umax_mono) & umax_mono < THR$umax_rate_min
  ) %>%
  dplyr::summarise(
    p_K_floor = mean(K_floor, na.rm = TRUE),
    p_umax_rate_floor = mean(umax_rate_floor, na.rm = TRUE),
    n_missing_baseline = sum(is.na(K_mono) | is.na(umax_mono))
  ) %>%
  print_tbl()

cat("\n================ MONOCULTURE BASELINE COMPLETENESS ================\n")

mono_day %>%
  dplyr::summarise(
    n_mono = dplyr::n(),
    n_K_mono_missing = sum(!is.finite(K_mono)),
    n_umax_mono_missing = sum(!is.finite(umax_mono))
  ) %>%
  print_tbl()

# 8) WT baseline construction on log2 scale

wt_source_debug <- fits_day_tag %>%
  dplyr::filter(Condition == "Competing_wt") %>%
  dplyr::filter(
    (Channel == "RFP" & GFP_Strain == "wt" & RFP_Strain != "wt") |
      (Channel == "GFP" & RFP_Strain == "wt" & GFP_Strain != "wt")
  ) %>%
  dplyr::mutate(
    Focal_Strain = ifelse(Channel == "RFP", RFP_Strain, GFP_Strain)
  ) %>%
  dplyr::left_join(
    mono_day,
    by = c("Date", "Channel", "Focal_Strain")
  ) %>%
  dplyr::mutate(
    K_ratio = K_day / pmax(K_mono, THR$K_MIN),
    umax_ratio = pmax(umax_day, THR$umax_rate_min) /
      pmax(umax_mono, THR$umax_rate_min),
    dK_pct = 100 * (K_ratio - 1),
    dumax_pct = 100 * (umax_ratio - 1),
    log2K = log2(pmax(K_ratio, THR$ratio_floor)),
    log2umax = log2(pmax(umax_ratio, THR$ratio_floor))
  )

wt_tagged <- wt_source_debug %>%
  dplyr::group_by(Channel, Focal_Strain) %>%
  dplyr::mutate(
    n_group = dplyr::n(),
    z_log2K = ifelse(n_group >= 5, mad_z(log2K), NA_real_),
    z_log2umax = ifelse(n_group >= 5, mad_z(log2umax), NA_real_)
  ) %>%
  dplyr::ungroup()

wt_source_filtered <- wt_tagged %>%
  dplyr::filter(
    dplyr::between(K_ratio, 1e-3, 1e3),
    dplyr::between(umax_ratio, 1e-3, 1e3),
    abs(z_log2K) <= 3 | is.na(z_log2K),
    abs(z_log2umax) <= 3 | is.na(z_log2umax)
  )

wt_day <- wt_source_filtered %>%
  dplyr::group_by(Channel, Focal_Strain) %>%
  dplyr::summarise(
    E_log2K = median(log2K, na.rm = TRUE),
    E_log2umax = median(log2umax, na.rm = TRUE),
    E_K_pct = median(dK_pct, na.rm = TRUE),
    E_umax_pct = median(dumax_pct, na.rm = TRUE),
    n_days = dplyr::n(),
    .groups = "drop"
  )

# 9) Final WT-benchmarked log2 indices

indices <- fits_norm %>%
  dplyr::filter(
    Condition == "Competing",
    RFP_Strain != "wt",
    GFP_Strain != "wt"
  ) %>%
  dplyr::left_join(
    wt_day,
    by = c("Channel", "Focal_Strain")
  ) %>%
  dplyr::mutate(
    I_log2K = log2K - E_log2K,
    I_log2umax = log2umax - E_log2umax,
    I_K_pct = dK_pct - E_K_pct,
    I_umax_pct = dumax_pct - E_umax_pct,
    I_K = I_log2K,
    I_umax = I_log2umax
  )

indices_keys <- indices %>%
  dplyr::group_by(Date, RFP_Strain, GFP_Strain) %>%
  dplyr::summarise(
    n_channels = dplyr::n_distinct(Channel),
    n_rows = dplyr::n(),
    .groups = "drop"
  )

if (nrow(indices) != 564L ||
    nrow(indices_keys) != 282L ||
    any(indices_keys$n_channels != 2L) ||
    any(indices_keys$n_rows != 2L) ||
    nrow(dplyr::anti_join(accepted_dates, indices_keys,
                          by = required_keys)) != 0L) {
  stop("Expected 282 accepted experimental dates, each with RFP and GFP indices.")
}

cat("\n================ FINAL LOG2 INDEX CHECK ================\n")

indices %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_Ilog2K_na = sum(is.na(I_log2K)),
    n_Ilog2umax_na = sum(is.na(I_log2umax)),
    n_IK_alias_na = sum(is.na(I_K)),
    n_Iumax_alias_na = sum(is.na(I_umax))
  ) %>%
  print_tbl()

cat("\n================ LOG2 INDEX SUMMARY ================\n")

indices %>%
  dplyr::select(I_log2K, I_log2umax, I_K_pct, I_umax_pct) %>%
  tidyr::pivot_longer(dplyr::everything()) %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(
    n = sum(is.finite(value)),
    median = median(value, na.rm = TRUE),
    q10 = quantile(value, 0.10, na.rm = TRUE),
    q90 = quantile(value, 0.90, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    max = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print_tbl()

cat("\n================ MODEL-RELEVANT UMAX COMPLETENESS ================\n")

indices %>%
  dplyr::mutate(
    Species6 = dplyr::case_when(
      RFP_Strain %in% c("Pf") ~ "Pf",
      RFP_Strain %in% c("Pu") ~ "Pu",
      RFP_Strain %in% c("Ps") ~ "Ps",
      stringr::str_starts(RFP_Strain, "Ba") ~ "Ba",
      stringr::str_starts(RFP_Strain, "Bb") ~ "Bb",
      stringr::str_starts(RFP_Strain, "Bh") ~ "Bh",
      TRUE ~ as.character(RFP_Strain)
    ),
    Competitor = GFP_Strain
  ) %>%
  dplyr::group_by(Channel, Species6, Competitor) %>%
  dplyr::summarise(
    n = dplyr::n(),
    finite_I_log2umax = sum(is.finite(I_log2umax)),
    missing_I_log2umax = sum(!is.finite(I_log2umax)),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Channel, Species6, Competitor) %>%
  print_tbl(n = Inf)

# 10) Sanity checks

stopifnot(
  indices %>%
    dplyr::filter(RFP_Strain == "wt" | GFP_Strain == "wt") %>%
    nrow() == 0
)

cat("\n================ WT LOG2 BASELINES ================\n")

wt_day %>%
  dplyr::arrange(Channel, Focal_Strain) %>%
  print_tbl(n = Inf)

# 11) Optional WT baseline plot data

wt_long <- wt_source_filtered %>%
  dplyr::select(
    Channel,
    Focal_Strain,
    log2K,
    log2umax,
    dK_pct,
    dumax_pct,
    Date
  ) %>%
  tidyr::pivot_longer(
    c(log2K, log2umax, dK_pct, dumax_pct),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = factor(
      metric,
      levels = c("log2K", "log2umax", "dK_pct", "dumax_pct")
    )
  )

# 12) How often is K* used?

cat("\n================ K METHOD USE ================\n")

fits_final %>%
  dplyr::count(K_method) %>%
  dplyr::mutate(
    pct = 100 * n / sum(n)
  ) %>%
  dplyr::arrange(desc(n)) %>%
  print_tbl(n = Inf)

cat("\n================ K METHOD USE BY CHANNEL / CONDITION ================\n")

fits_final %>%
  dplyr::count(Channel, Condition, K_method) %>%
  dplyr::group_by(Channel, Condition) %>%
  dplyr::mutate(
    pct = 100 * n / sum(n)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Channel, Condition, desc(n)) %>%
  print_tbl(n = Inf)

cat("\n================ K METHOD USE FOR BIOLOGICALLY RELEVANT COMPETITION CURVES ================\n")

fits_final %>%
  dplyr::filter(
    Condition == "Competing",
    RFP_Strain != "wt",
    GFP_Strain != "wt",
    GFP_Strain %in% c("Sl", "Ss", "At", "Bp", "Ck", "Pk")
  ) %>%
  dplyr::count(Channel, RFP_Strain, GFP_Strain, K_method) %>%
  dplyr::group_by(Channel, RFP_Strain, GFP_Strain) %>%
  dplyr::mutate(
    pct = 100 * n / sum(n)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(pct), Channel, RFP_Strain, GFP_Strain) %>%
  print_tbl(n = Inf)

cat("\nK / mu_max processing complete. Model input: indices.\n")
