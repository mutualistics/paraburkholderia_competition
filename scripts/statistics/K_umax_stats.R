
# K and mu_max models, inference and diagnostics
# Run after 3_Kr_data_processing in the same R session.
# Model selection, EMMs and Holm tests are calculated here.
# No files are written.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(glmmTMB)
  library(tibble)
  library(emmeans)
  library(DHARMa)
})

if (!exists("indices", inherits = TRUE)) {
  stop("Run 01_K_umax_processing.R first (indices is missing).")
}

need_cols <- c(
  "Date",
  "Channel",
  "RFP_Strain",
  "GFP_Strain",
  "I_log2K",
  "I_log2umax"
)

miss <- setdiff(need_cols, names(indices))

if (length(miss) > 0) {
  stop("indices is missing required columns: ", paste(miss, collapse = ", "))
}

# The input contains one RFP and one GFP response for each of the
# 282 accepted experimental dates. Missing traits remain NA until modeling.

kr_input_dates <- indices %>%
  dplyr::group_by(Date, RFP_Strain, GFP_Strain) %>%
  dplyr::summarise(
    n_channels = dplyr::n_distinct(Channel),
    n_rows = dplyr::n(),
    .groups = "drop"
  )

if (nrow(indices) != 564L ||
    nrow(kr_input_dates) != 282L ||
    any(kr_input_dates$n_channels != 2L) ||
    any(kr_input_dates$n_rows != 2L)) {
  stop("Expected the K/mu_max indices for 282 accepted dates (564 channel records).")
}

old_objs <- c(
  "dat_rfp",
  "dat_gfp",
  "res_rfp_k",
  "res_gfp_k",
  "res_rfp_umax",
  "res_gfp_umax",
  "trait_results",
  "final_species6_rfp_k",
  "final_species6_gfp_k",
  "final_species6_rfp_umax",
  "final_species6_gfp_umax",
  "all_aic",
  "all_complexity",
  "best_summary",
  "winner_table",
  "final_convergence_check"
)

rm(list = intersect(old_objs, ls()))

# Model settings

FINAL_CHOICE <- "best"

species6_levels <- c("Pf", "Pu", "Ps", "Ba", "Bb", "Bh")
plot_group_levels <- c("Other", "Ba", "Bb", "Bh")
comp_levels <- c("Sl", "Ss", "At", "Bp", "Ck", "Pk")

umax_estimator_label <- "rolling_5pt_log_slope_0_24h_baseline_plus_2MAD_3of5"

print_tbl <- function(x, n = Inf) {
  x <- as.data.frame(x)
  
  if (is.finite(n) && nrow(x) > n) {
    x <- x[seq_len(n), , drop = FALSE]
  }
  
  print(x, row.names = FALSE)
  invisible(x)
}

safe_chr <- function(x) {
  if (is.null(x) || length(x) == 0) {
    NA_character_
  } else {
    as.character(x[1])
  }
}

make_species6 <- function(rfp_strain) {
  dplyr::case_when(
    rfp_strain == "Pf"                     ~ "Pf",
    rfp_strain == "Pu"                     ~ "Pu",
    rfp_strain == "Ps"                     ~ "Ps",
    stringr::str_detect(rfp_strain, "^Ba") ~ "Ba",
    stringr::str_detect(rfp_strain, "^Bb") ~ "Bb",
    stringr::str_detect(rfp_strain, "^Bh") ~ "Bh",
    TRUE                                   ~ NA_character_
  )
}

model_stats <- function(model,
                        model_name = NA_character_) {
  
  if (is.null(model)) {
    return(
      tibble(
        Model = model_name,
        npar = NA_integer_,
        AIC = NA_real_,
        BIC = NA_real_,
        logLik = NA_real_,
        convergence = NA_integer_,
        pdHess = FALSE,
        stable = FALSE
      )
    )
  }
  
  conv <- tryCatch(
    as.integer(model$fit$convergence),
    error = function(e) NA_integer_
  )
  
  pd <- tryCatch(
    isTRUE(model$sdr$pdHess),
    error = function(e) FALSE
  )
  
  ll <- tryCatch(
    as.numeric(logLik(model)),
    error = function(e) NA_real_
  )
  
  aic <- tryCatch(
    AIC(model),
    error = function(e) NA_real_
  )
  
  bic <- tryCatch(
    BIC(model),
    error = function(e) NA_real_
  )
  
  npar <- tryCatch(
    attr(logLik(model), "df"),
    error = function(e) NA_integer_
  )
  
  tibble(
    Model = model_name,
    npar = npar,
    AIC = aic,
    BIC = bic,
    logLik = ll,
    convergence = conv,
    pdHess = pd,
    stable = isTRUE(conv == 0) && isTRUE(pd) && is.finite(ll)
  )
}

fit_tmb_safe <- function(cond_formula,
                         disp_formula,
                         dat,
                         reml = FALSE,
                         model_name = "") {
  tryCatch(
    glmmTMB(
      formula = cond_formula,
      dispformula = disp_formula,
      family = gaussian(),
      data = dat,
      REML = reml,
      control = glmmTMBControl(
        optCtrl = list(iter.max = 1e4, eval.max = 1e4)
      )
    ),
    error = function(e) {
      message("Fit failed for ", model_name, ": ", conditionMessage(e))
      NULL
    }
  )
}

score_models <- function(model_list,
                         spec_list) {
  
  keep <- !vapply(model_list, is.null, logical(1))
  model_list <- model_list[keep]
  
  raw_score <- bind_rows(
    lapply(names(model_list), function(nm) {
      model_stats(model_list[[nm]], nm) %>%
        mutate(
          Description = spec_list[[nm]]$desc,
          Complexity_rank = spec_list[[nm]]$complexity_rank,
          .after = Model
        )
    })
  )
  
  score_tab <- raw_score %>%
    filter(
      is.finite(AIC),
      is.finite(BIC),
      is.finite(logLik)
    ) %>%
    mutate(
      delta_AIC = AIC - min(AIC, na.rm = TRUE),
      delta_BIC = BIC - min(BIC, na.rm = TRUE),
      within_2_AIC = delta_AIC <= 2
    ) %>%
    arrange(AIC)
  
  if (nrow(score_tab) == 0) {
    stop("No finite candidate models remained after scoring.")
  }
  
  score_tab
}

print_nested_if_possible <- function(mod_small,
                                     mod_big,
                                     label) {
  
  small_ok <- !is.null(mod_small) && is.finite(suppressWarnings(AIC(mod_small)))
  big_ok <- !is.null(mod_big) && is.finite(suppressWarnings(AIC(mod_big)))
  
  cat("\n---", label, "---\n")
  
  if (small_ok && big_ok) {
    print(anova(mod_small, mod_big))
  } else {
    cat("Skipped: at least one model failed or had nonfinite AIC.\n")
  }
}

make_model_data <- function(dat,
                            channel_keep) {
  dat %>%
    filter(
      Channel == channel_keep,
      GFP_Strain %in% comp_levels,
      RFP_Strain != "wt",
      GFP_Strain != "wt"
    ) %>%
    mutate(
      Species6 = make_species6(RFP_Strain),
      Competitor = GFP_Strain
    ) %>%
    filter(!is.na(Species6)) %>%
    select(
      Date,
      RFP_Strain,
      Species6,
      Competitor,
      I_log2K,
      I_log2umax
    ) %>%
    pivot_longer(
      cols = c(I_log2K, I_log2umax),
      names_to = "trait",
      values_to = "y"
    ) %>%
    mutate(
      trait = recode(
        trait,
        "I_log2K" = "K",
        "I_log2umax" = "umax"
      )
    ) %>%
    filter(is.finite(y)) %>%
    mutate(
      Date = factor(Date),
      RFP_Strain = factor(RFP_Strain),
      Species6 = factor(Species6, levels = species6_levels),
      Competitor = factor(Competitor, levels = comp_levels),
      trait = factor(trait, levels = c("K", "umax"))
    ) %>%
    droplevels()
}

cat("\n================ UPSTREAM METHOD CHECK ================\n")

cat("\nExpected μmax estimator:\n")
cat(umax_estimator_label, "\n")

if (exists("fits_final") && "umax_method" %in% names(fits_final)) {
  cat("\nObserved umax_method values in fits_final:\n")
  
  fits_final %>%
    count(umax_method) %>%
    print_tbl(n = Inf)
}

dat_rfp <- make_model_data(indices, "RFP")
dat_gfp <- make_model_data(indices, "GFP")

cat("\n================ INPUT DATA CHECK ================\n")

cat("\nRows in dat_rfp:\n")
print(nrow(dat_rfp))

cat("\nRows in dat_gfp:\n")
print(nrow(dat_gfp))

cat("\nRFP rows by trait:\n")
print(table(dat_rfp$trait))

cat("\nGFP rows by trait:\n")
print(table(dat_gfp$trait))

cat("\nRFP rows by Species6 and Competitor:\n")
print(table(dat_rfp$Species6, dat_rfp$Competitor))

cat("\nGFP rows by Species6 and Competitor:\n")
print(table(dat_gfp$Species6, dat_gfp$Competitor))

cat("\nResponse scale check:\n")

bind_rows(
  dat_rfp %>% mutate(Channel = "RFP"),
  dat_gfp %>% mutate(Channel = "GFP")
) %>%
  group_by(Channel, trait) %>%
  summarise(
    n = n(),
    median = median(y, na.rm = TRUE),
    q10 = quantile(y, 0.10, na.rm = TRUE),
    q90 = quantile(y, 0.90, na.rm = TRUE),
    min = min(y, na.rm = TRUE),
    max = max(y, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print_tbl()

# Candidate models

trait_model_specs <- list(
  m0_main_date_homo = list(
    cond = y ~ Species6 + Competitor + (1 | Date),
    disp = ~ 1,
    desc = "Main effects + Date; homogeneous variance",
    complexity_rank = 1
  ),
  m1_int_date_homo = list(
    cond = y ~ Species6 * Competitor + (1 | Date),
    disp = ~ 1,
    desc = "Species6 * Competitor + Date; homogeneous variance",
    complexity_rank = 2
  ),
  m2_int_date_strain_homo = list(
    cond = y ~ Species6 * Competitor + (1 | Date) + (1 | RFP_Strain),
    disp = ~ 1,
    desc = "Species6 * Competitor + Date + RFP_Strain; homogeneous variance",
    complexity_rank = 3
  ),
  m3_int_date_varS = list(
    cond = y ~ Species6 * Competitor + (1 | Date),
    disp = ~ Species6,
    desc = "Species6 * Competitor + Date; dispersion ~ Species6",
    complexity_rank = 4
  ),
  m4_int_date_varC = list(
    cond = y ~ Species6 * Competitor + (1 | Date),
    disp = ~ Competitor,
    desc = "Species6 * Competitor + Date; dispersion ~ Competitor",
    complexity_rank = 5
  ),
  m5_int_date_varSC = list(
    cond = y ~ Species6 * Competitor + (1 | Date),
    disp = ~ Species6 + Competitor,
    desc = "Species6 * Competitor + Date; dispersion ~ Species6 + Competitor",
    complexity_rank = 6
  ),
  m6_int_date_strain_varSC = list(
    cond = y ~ Species6 * Competitor + (1 | Date) + (1 | RFP_Strain),
    disp = ~ Species6 + Competitor,
    desc = "Species6 * Competitor + Date + RFP_Strain; dispersion ~ Species6 + Competitor",
    complexity_rank = 7
  )
)

make_selection_order <- function(score_tab,
                                 final_choice = "best") {
  
  stable_tab <- score_tab %>%
    filter(stable)
  
  if (nrow(stable_tab) == 0) {
    warning("No ML candidates had stable Hessian/convergence; using finite AIC candidates as fallback.")
    stable_tab <- score_tab
  }
  
  best_stable <- stable_tab %>%
    arrange(AIC) %>%
    slice(1) %>%
    pull(Model)
  
  simplest_likely <- stable_tab %>%
    filter(within_2_AIC) %>%
    arrange(Complexity_rank, AIC) %>%
    slice(1) %>%
    pull(Model)
  
  if (length(simplest_likely) == 0) {
    simplest_likely <- best_stable
  }
  
  primary <- if (final_choice == "simplest_likely") {
    simplest_likely
  } else {
    best_stable
  }
  
  unique(
    c(
      primary,
      stable_tab %>%
        arrange(AIC) %>%
        pull(Model),
      score_tab %>%
        arrange(AIC) %>%
        pull(Model)
    )
  )
}

fit_final_reml_stable <- function(selection_order,
                                  spec_list,
                                  dat,
                                  channel_name,
                                  trait_name) {
  
  attempts <- tibble()
  attempt_models <- list()
  
  for (nm in selection_order) {
    
    cat("\nTrying final REML candidate:", nm, "\n")
    
    mod <- fit_tmb_safe(
      cond_formula = spec_list[[nm]]$cond,
      disp_formula = spec_list[[nm]]$disp,
      dat = dat,
      reml = TRUE,
      model_name = paste(channel_name, trait_name, nm, "REML", sep = "_")
    )
    
    attempt_models[[nm]] <- mod
    
    chk <- model_stats(mod, nm)
    attempts <- bind_rows(attempts, chk)
    
    print_tbl(chk)
    
    if (isTRUE(chk$stable[1])) {
      return(
        list(
          final_name = nm,
          final_reml = mod,
          final_check = attempts,
          final_stable = TRUE,
          final_was_fallback = nm != selection_order[1]
        )
      )
    }
    
    cat("Skipping ", nm, " as final model because REML fit was unstable.\n", sep = "")
  }
  
  usable <- attempts %>%
    filter(is.finite(AIC), is.finite(logLik)) %>%
    arrange(AIC)
  
  if (nrow(usable) == 0) {
    stop("No usable REML final model could be fit for ", channel_name, " ", trait_name)
  }
  
  fallback_name <- usable$Model[1]
  
  warning(
    "No stable REML model found for ",
    channel_name,
    " ",
    trait_name,
    ". Using lowest-AIC usable REML fit anyway: ",
    fallback_name
  )
  
  list(
    final_name = fallback_name,
    final_reml = attempt_models[[fallback_name]],
    final_check = attempts,
    final_stable = FALSE,
    final_was_fallback = fallback_name != selection_order[1]
  )
}

fit_trait_model_set <- function(dat,
                                channel_name,
                                trait_name) {
  
  dat_sub <- dat %>%
    filter(trait == trait_name) %>%
    droplevels() %>%
    mutate(
      Date = factor(Date),
      RFP_Strain = factor(RFP_Strain),
      Species6 = factor(Species6, levels = species6_levels),
      Competitor = factor(Competitor, levels = comp_levels)
    )
  
  contrasts(dat_sub$Species6) <- contr.sum(nlevels(dat_sub$Species6))
  contrasts(dat_sub$Competitor) <- contr.sum(nlevels(dat_sub$Competitor))
  
  if (nrow(dat_sub) == 0) {
    stop("No rows for channel = ", channel_name, ", trait = ", trait_name)
  }
  
  cat("\n==========================================================\n")
  cat("CHANNEL:", channel_name, "| TRAIT:", trait_name, "\n")
  cat("Scale: WT-benchmarked log2 fold-change\n")
  cat("Rows:", nrow(dat_sub), "\n")
  cat("==========================================================\n")
  
  ml_models <- lapply(
    names(trait_model_specs),
    function(nm) {
      fit_tmb_safe(
        cond_formula = trait_model_specs[[nm]]$cond,
        disp_formula = trait_model_specs[[nm]]$disp,
        dat = dat_sub,
        reml = FALSE,
        model_name = paste(channel_name, trait_name, nm, sep = "_")
      )
    }
  )
  
  names(ml_models) <- names(trait_model_specs)
  
  score_tab_aic <- score_models(ml_models, trait_model_specs)
  
  score_tab_complexity <- score_tab_aic %>%
    arrange(Complexity_rank, AIC)
  
  cat("\n================ MODEL COMPARISON: SORTED BY AIC ================\n")
  print_tbl(score_tab_aic, n = Inf)
  
  cat("\n================ MODEL COMPARISON: SORTED BY COMPLEXITY ================\n")
  print_tbl(score_tab_complexity, n = Inf)
  
  print_nested_if_possible(
    ml_models$m0_main_date_homo,
    ml_models$m1_int_date_homo,
    paste(channel_name, trait_name, "main vs interaction")
  )
  
  print_nested_if_possible(
    ml_models$m1_int_date_homo,
    ml_models$m2_int_date_strain_homo,
    paste(channel_name, trait_name, "add RFP_Strain under homogeneous variance")
  )
  
  print_nested_if_possible(
    ml_models$m1_int_date_homo,
    ml_models$m3_int_date_varS,
    paste(channel_name, trait_name, "add variance by Species6")
  )
  
  print_nested_if_possible(
    ml_models$m1_int_date_homo,
    ml_models$m4_int_date_varC,
    paste(channel_name, trait_name, "add variance by Competitor")
  )
  
  print_nested_if_possible(
    ml_models$m5_int_date_varSC,
    ml_models$m6_int_date_strain_varSC,
    paste(channel_name, trait_name, "add RFP_Strain on top of variance by Species6 + Competitor")
  )
  
  best_name <- score_tab_aic %>%
    arrange(AIC) %>%
    slice(1) %>%
    pull(Model)
  
  best_stable_name <- score_tab_aic %>%
    filter(stable) %>%
    arrange(AIC) %>%
    slice(1) %>%
    pull(Model)
  
  if (length(best_stable_name) == 0) {
    best_stable_name <- NA_character_
  }
  
  simplest_likely_name <- score_tab_aic %>%
    filter(stable, within_2_AIC) %>%
    arrange(Complexity_rank, AIC) %>%
    slice(1) %>%
    pull(Model)
  
  if (length(simplest_likely_name) == 0) {
    simplest_likely_name <- best_stable_name
  }
  
  cat("\nBest finite model by AIC:", safe_chr(best_name), "\n")
  cat("Best stable model by AIC:", safe_chr(best_stable_name), "\n")
  cat("Simplest stable model within delta AIC <= 2:", safe_chr(simplest_likely_name), "\n")
  
  selection_order <- make_selection_order(
    score_tab = score_tab_aic,
    final_choice = FINAL_CHOICE
  )
  
  cat("\nFinal REML candidate order:\n")
  print(selection_order)
  
  final_fit <- fit_final_reml_stable(
    selection_order = selection_order,
    spec_list = trait_model_specs,
    dat = dat_sub,
    channel_name = channel_name,
    trait_name = trait_name
  )
  
  cat("\nSelected final REML model:", final_fit$final_name, "\n")
  cat("Final model stable:", final_fit$final_stable, "\n")
  cat("Final model was fallback:", final_fit$final_was_fallback, "\n")
  
  cat("\n================ FINAL REML MODEL SUMMARY ================\n")
  print(summary(final_fit$final_reml))
  
  invisible(
    list(
      data = dat_sub,
      ml_models = ml_models,
      score_aic = score_tab_aic,
      score_complexity = score_tab_complexity,
      best_name = safe_chr(best_name),
      best_stable_name = safe_chr(best_stable_name),
      simplest_likely_name = safe_chr(simplest_likely_name),
      final_name = final_fit$final_name,
      final_reml = final_fit$final_reml,
      final_check = final_fit$final_check,
      final_stable = final_fit$final_stable,
      final_was_fallback = final_fit$final_was_fallback,
      scale = "WT-benchmarked log2 fold-change"
    )
  )
}

res_rfp_k <- fit_trait_model_set(dat_rfp, "RFP", "K")
res_gfp_k <- fit_trait_model_set(dat_gfp, "GFP", "K")
res_rfp_umax <- fit_trait_model_set(dat_rfp, "RFP", "umax")
res_gfp_umax <- fit_trait_model_set(dat_gfp, "GFP", "umax")

trait_results <- list(
  RFP_K = res_rfp_k,
  GFP_K = res_gfp_k,
  RFP_umax = res_rfp_umax,
  GFP_umax = res_gfp_umax
)

final_species6_rfp_k <- res_rfp_k$final_reml
final_species6_gfp_k <- res_gfp_k$final_reml
final_species6_rfp_umax <- res_rfp_umax$final_reml
final_species6_gfp_umax <- res_gfp_umax$final_reml

all_aic <- bind_rows(
  lapply(names(trait_results), function(nm) {
    trait_results[[nm]]$score_aic %>%
      mutate(Set = nm)
  })
)

all_complexity <- bind_rows(
  lapply(names(trait_results), function(nm) {
    trait_results[[nm]]$score_complexity %>%
      mutate(Set = nm)
  })
)

best_summary <- bind_rows(
  lapply(names(trait_results), function(nm) {
    x <- trait_results[[nm]]
    
    tibble(
      Set = nm,
      Best_by_AIC = safe_chr(x$best_name),
      Best_stable_by_AIC = safe_chr(x$best_stable_name),
      Simplest_likely = safe_chr(x$simplest_likely_name),
      Final_selected = safe_chr(x$final_name),
      Final_stable = isTRUE(x$final_stable),
      Final_was_fallback = isTRUE(x$final_was_fallback)
    )
  })
)

winner_table <- bind_rows(
  all_aic %>%
    group_by(Set) %>%
    slice_min(AIC, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      Set,
      Criterion = "AIC",
      Winner = Model,
      Value = AIC
    ),
  all_aic %>%
    group_by(Set) %>%
    slice_min(BIC, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      Set,
      Criterion = "BIC",
      Winner = Model,
      Value = BIC
    )
) %>%
  arrange(Set, Criterion)

final_convergence_check <- bind_rows(
  lapply(names(trait_results), function(nm) {
    model_stats(
      trait_results[[nm]]$final_reml,
      model_name = nm
    ) %>%
      mutate(
        Final_selected = trait_results[[nm]]$final_name,
        Final_was_fallback = trait_results[[nm]]$final_was_fallback,
        .after = Model
      )
  })
)

cat("\n================ SPECIES6 TRAIT-SPECIFIC BEST MODELS ================\n")
print_tbl(best_summary, n = Inf)

cat("\n================ AIC / BIC WINNERS ================\n")
print_tbl(winner_table, n = Inf)

cat("\n================ FINAL CONVERGENCE CHECK ================\n")
print_tbl(final_convergence_check, n = Inf)


# Estimated marginal means and contrasts
# Each competitor is weighted equally; Non-AA is the equal mean of Pf, Pu and Ps.

kr_groups <- c("Non-AA", "Ba", "Bb", "Bh")
kr_responses <- c("Paraburkholderia response", "Competitor response")
kr_traits <- c("K", "umax")

kr_models <- list(
  RFP_K = res_rfp_k,
  GFP_K = res_gfp_k,
  RFP_umax = res_rfp_umax,
  GFP_umax = res_gfp_umax
)

kr_keys <- tibble::tibble(
  key = names(kr_models),
  Response = rep(kr_responses, times = 2),
  Trait = rep(kr_traits, each = 2)
)

kr_clean_label <- function(x) {
  sub("^\\((.*)\\)$", "\\1", trimws(as.character(x)))
}

kr_ci <- function(x) {
  lo <- intersect(c("lower.CL", "asymp.LCL", "LCL"), names(x))[1]
  hi <- intersect(c("upper.CL", "asymp.UCL", "UCL"), names(x))[1]
  if (is.na(lo) || is.na(hi)) stop("Model EMM confidence limits are missing.")
  list(lo = x[[lo]], hi = x[[hi]])
}

kr_grid <- function(mod, within = FALSE) {
  em6 <- if (within) {
    emmeans(mod, ~ Species6 | Competitor)
  } else {
    emmeans(mod, ~ Species6, weights = "equal")
  }
  
  lev <- unique(as.character(as.data.frame(em6)$Species6))
  
  if (length(lev) != 6L || !setequal(lev, species6_levels)) {
    stop("Six RFP species required in model EMM grid.")
  }
  
  w <- function(members) {
    x <- setNames(rep(0, length(lev)), lev)
    x[members] <- 1 / length(members)
    unname(x)
  }
  
  methods <- list(
    "Non-AA" = w(c("Pf", "Pu", "Ps")),
    "Ba" = w("Ba"),
    "Bb" = w("Bb"),
    "Bh" = w("Bh")
  )
  
  contrast(em6, method = methods, by = if (within) "Competitor" else NULL)
}

kr_inference <- function(result, response, trait, within = FALSE) {
  em4 <- kr_grid(result$final_reml, within)
  e <- as.data.frame(summary(em4, infer = c(TRUE, TRUE), null = 0, adjust = "none"))
  ci <- kr_ci(e)
  
  n_emm <- if (within) 24L else 4L
  
  if (nrow(e) != n_emm || any(!is.finite(e$p.value))) {
    stop("Unexpected EMM count or non-finite p-values for ", response, "/", trait)
  }
  
  emm <- e %>%
    mutate(
      Plot_Group = kr_clean_label(contrast),
      emmean_log2 = estimate,
      lower_log2 = ci$lo,
      upper_log2 = ci$hi,
      p_vs_zero_raw = p.value,
      p_vs_zero_holm = p.adjust(p.value, method = "holm"),
      sig_vs_zero = p_vs_zero_holm < 0.05,
      Response = response,
      Trait = trait
    )
  
  if (!setequal(emm$Plot_Group, kr_groups) ||
      anyDuplicated(emm[c(intersect("Competitor", names(emm)), "Plot_Group")])) {
    stop("Unexpected group labels or duplicates in EMMs.")
  }
  
  pa <- as.data.frame(
    pairs(
      em4,
      by = if (within) "Competitor" else NULL,
      adjust = "none"
    )
  )
  
  n_pairs <- if (within) 36L else 6L
  
  if (nrow(pa) != n_pairs || any(!is.finite(pa$p.value))) {
    stop("Unexpected number of pairwise contrasts for ", response, "/", trait)
  }
  
  if (within) {
    pa <- pa %>%
      group_by(Competitor) %>%
      mutate(p_holm = p.adjust(p.value, method = "holm")) %>%
      ungroup()
    
    if (any(table(pa$Competitor) != 6L)) {
      stop("Expected six pairs per competitor.")
    }
  } else {
    pa$p_holm <- p.adjust(pa$p.value, method = "holm")
  }
  
  parts <- strsplit(as.character(pa$contrast), " - ", fixed = TRUE)
  
  if (any(lengths(parts) != 2L)) {
    stop("Cannot parse pairwise contrast labels.")
  }
  
  pa <- pa %>%
    mutate(
      group1 = kr_clean_label(vapply(parts, `[`, character(1), 1L)),
      group2 = kr_clean_label(vapply(parts, `[`, character(1), 2L)),
      label = case_when(
        p_holm < 0.001 ~ "***",
        p_holm < 0.01 ~ "**",
        p_holm < 0.05 ~ "*",
        TRUE ~ NA_character_
      ),
      Response = response,
      Trait = trait
    )
  
  if (any(!pa$group1 %in% kr_groups) || any(!pa$group2 %in% kr_groups)) {
    stop("Unrecognised group name in pairwise contrast.")
  }
  
  list(emm = emm, pairs = pa)
}

kr_raw_data <- function(dat, response, trait) {
  dat %>%
    mutate(
      Plot_Group = case_when(
        as.character(Species6) %in% c("Pf", "Pu", "Ps") ~ "Non-AA",
        as.character(Species6) %in% c("Ba", "Bb", "Bh") ~ as.character(Species6),
        TRUE ~ NA_character_
      ),
      y_log2 = y,
      Response = response,
      Trait = trait
    ) %>%
    filter(!is.na(Plot_Group), is.finite(y_log2)) %>%
    mutate(
      Plot_Group = factor(Plot_Group, levels = kr_groups),
      Competitor = factor(as.character(Competitor), levels = comp_levels),
      Response = factor(Response, levels = kr_responses),
      Trait = factor(Trait, levels = kr_traits)
    )
}

kr_emm_results <- lapply(
  seq_len(nrow(kr_keys)),
  function(i) {
    key <- kr_keys$key[i]
    
    x <- kr_inference(
      kr_models[[key]],
      kr_keys$Response[i],
      kr_keys$Trait[i],
      FALSE
    )
    
    y <- kr_inference(
      kr_models[[key]],
      kr_keys$Response[i],
      kr_keys$Trait[i],
      TRUE
    )
    
    list(overall = x, within = y)
  }
)

names(kr_emm_results) <- kr_keys$key

kr_overall_emm <- bind_rows(
  lapply(kr_emm_results, function(x) x$overall$emm)
)

kr_overall_pairs <- bind_rows(
  lapply(kr_emm_results, function(x) x$overall$pairs)
)

kr_within_emm <- bind_rows(
  lapply(kr_emm_results, function(x) x$within$emm)
)

kr_within_pairs <- bind_rows(
  lapply(kr_emm_results, function(x) x$within$pairs)
)

kr_raw <- bind_rows(
  lapply(seq_len(nrow(kr_keys)), function(i) {
    kr_raw_data(
      kr_models[[kr_keys$key[i]]]$data,
      kr_keys$Response[i],
      kr_keys$Trait[i]
    )
  })
)

# Descriptive competitor-specific means for Figure 4's underlaid points.

kr_raw_overall <- kr_raw %>%
  group_by(Trait, Response, Plot_Group, Competitor) %>%
  summarise(
    y_log2 = mean(y_log2),
    n_obs = n(),
    .groups = "drop"
  )

kr_raw_within <- kr_raw

cat("\nOVERALL EMMs (Holm across 4 per response/trait)\n")
print(
  as.data.frame(
    kr_overall_emm[, c(
      "Trait", "Response", "Plot_Group",
      "emmean_log2", "lower_log2", "upper_log2", "p_vs_zero_holm"
    )]
  ),
  row.names = FALSE
)

cat("\nOVERALL PAIRWISE (Holm across 6 per response/trait)\n")
print(
  as.data.frame(
    kr_overall_pairs[, c(
      "Trait", "Response", "contrast", "estimate", "SE", "p_holm"
    )]
  ),
  row.names = FALSE
)

cat("\nWITHIN-COMPETITOR TESTS: Holm across 24 vs-zero and 6 pairs/competitor\n")
print(
  kr_within_emm %>%
    count(Trait, Response, sig_vs_zero)
)

print(
  kr_within_pairs %>%
    count(Trait, Response, significant = !is.na(label))
)

# Residual diagnostics for the four selected models

RUN_DIAGNOSTICS <- TRUE
SHOW_DIAGNOSTIC_PLOTS <- FALSE
DIAG_NSIM <- 1000L
PEARSON_CUTOFF <- 3

diag_table <- NULL
extreme_residual_rows <- NULL
diag_objects <- NULL

if (RUN_DIAGNOSTICS) {
  diag_objects <- lapply(
    seq_along(kr_models),
    function(i) {
      r <- kr_models[[i]]
      set.seed(123 + i)
      
      sim <- DHARMa::simulateResiduals(
        r$final_reml,
        n = DIAG_NSIM,
        plot = FALSE
      )
      
      pearson <- as.numeric(
        residuals(r$final_reml, type = "pearson")
      )
      
      fitted <- as.numeric(
        fitted(r$final_reml)
      )
      
      if (length(pearson) != nrow(r$data) ||
          length(fitted) != nrow(r$data)) {
        stop(
          "Model residual lengths differ from fitted data: ",
          names(kr_models)[i]
        )
      }
      
      dat <- r$data %>%
        mutate(
          .pearson = pearson,
          .fitted = fitted,
          .dharma = as.numeric(sim$scaledResiduals)
        )
      
      list(
        sim = sim,
        dat = dat,
        table = tibble(
          Model = names(kr_models)[i],
          Selected = r$final_name,
          n = nrow(dat),
          Uniformity_p = DHARMa::testUniformity(sim)$p.value,
          Dispersion_p = DHARMa::testDispersion(sim)$p.value,
          Outlier_p = DHARMa::testOutliers(sim)$p.value,
          N_abs_Pearson_gt_3 = sum(
            abs(pearson) > PEARSON_CUTOFF,
            na.rm = TRUE
          )
        )
      )
    }
  )
  
  names(diag_objects) <- names(kr_models)
  
  diag_table <- bind_rows(
    lapply(diag_objects, function(x) x$table)
  )
  
  extreme_residual_rows <- bind_rows(
    lapply(names(diag_objects), function(k) {
      diag_objects[[k]]$dat %>%
        filter(abs(.pearson) > PEARSON_CUTOFF) %>%
        mutate(Model = k)
    })
  )
  
  cat("\nMODEL DIAGNOSTICS\n")
  print(as.data.frame(diag_table), row.names = FALSE)
  
  if (nrow(extreme_residual_rows)) {
    cat("\nExtreme Pearson residuals\n")
    
    print(
      as.data.frame(
        extreme_residual_rows[, c(
          "Model", "Date", "RFP_Strain",
          "Species6", "Competitor", "y", ".pearson"
        )]
      ),
      row.names = FALSE
    )
  }
  
  if (SHOW_DIAGNOSTIC_PLOTS) {
    for (x in diag_objects) plot(x$sim)
  }
}

cat("\nK / mu_max model analysis complete.\n")
