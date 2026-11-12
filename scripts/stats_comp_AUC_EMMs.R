# AUC-ratio competition model and estimated marginal means
# Run after 1_FullDataImport_all_valid.R.
# Input: analysis_data_day_all_valid.rds (282 date-level observations).
# Fits the mixed model and calculates contrasts for Figures 2 and 3.
# Figure2.R and Figure3_plot_only_FINAL.R use the objects created here.
# No files are saved unless SAVE_TABLES or SAVE_MODEL_RDS is TRUE.

suppressPackageStartupMessages({
  library(dplyr)
  library(glmmTMB)
  library(emmeans)
  library(DHARMa)
})

# Settings

INPUT_CANDIDATES <- c(
  "analysis_data_day_all_valid.rds",
  file.path("data", "processed", "analysis_data_day_all_valid.rds")
)

CHECK_DISPERSION_MODELS <- TRUE
RUN_DHARMA_PLOTS        <- FALSE
RUN_FIG3_SENSITIVITIES <- TRUE
RUN_NO_BA_REFIT         <- TRUE
SAVE_TABLES             <- FALSE
SAVE_MODEL_RDS          <- FALSE
OUTPUT_DIR              <- file.path("outputs", "competition_model")
DHARMA_SEED             <- 123
DHARMA_SIMULATIONS      <- 1000L

species6 <- c("Pf", "Pu", "Ps", "Ba", "Bb", "Bh")
competitors <- c("Sl", "Ss", "At", "Bp", "Ck", "Pk")
plot_groups <- c("Non-AA", "Ba", "Bb", "Bh")


# 1. Load data

input_found <- INPUT_CANDIDATES[file.exists(INPUT_CANDIDATES)]

if (!length(input_found)) {
  stop(
    "Could not find analysis_data_day_all_valid.rds. ",
    "Run the data import first."
  )
}

INPUT_FILE <- input_found[[1]]
analysis_data_day <- readRDS(INPUT_FILE)

cat("\nInput data:", normalizePath(INPUT_FILE), "\n")

required_cols <- c(
  "Date", "RFP_Strain", "RFP_Species6",
  "GFP_Strain", "log_ratio"
)

missing_cols <- setdiff(required_cols, names(analysis_data_day))

if (length(missing_cols)) {
  stop(
    "Missing required column(s): ",
    paste(missing_cols, collapse = ", ")
  )
}

# Use all six RFP species.
# Non-AA is averaged at the EMM stage.

dat <- analysis_data_day %>%
  mutate(
    Date = factor(Date),
    RFP_Strain = factor(RFP_Strain),
    RFP_Species6 = factor(
      as.character(RFP_Species6),
      levels = species6
    ),
    GFP_Strain = factor(
      as.character(GFP_Strain),
      levels = competitors
    )
  ) %>%
  droplevels()

if (anyNA(dat[, required_cols]) ||
    any(!is.finite(dat$log_ratio))) {
  stop("Missing or non-finite model observations.")
}

if (!setequal(levels(dat$RFP_Species6), species6) ||
    !setequal(levels(dat$GFP_Strain), competitors)) {
  stop("Not all six RFP species and six GFP competitors are represented.")
}

if (any(table(dat$RFP_Species6, dat$GFP_Strain) == 0L)) {
  stop("At least one species/competitor combination is missing.")
}

if (anyDuplicated(dat[c(
  "Date", "RFP_Strain", "GFP_Strain"
)])) {
  stop("Expected one observation per strain/competitor/date.")
}

# One biological observation per date/strain/competitor.
# The dataset contains 282 observations from 15 RFP strains.

if (nrow(dat) != 282L ||
    nlevels(dat$RFP_Strain) != 15L) {
  stop(
    "Expected 282 date-level observations and 15 RFP strains. ",
    "Check the import."
  )
}

replicate_counts <- with(
  dat,
  table(RFP_Strain, GFP_Strain)
)

if (any(replicate_counts < 3L |
        replicate_counts > 5L)) {
  stop("Expected three to five dates per strain/competitor.")
}

# Sum contrasts for Type III tests.

contrasts(dat$RFP_Species6) <-
  contr.sum(nlevels(dat$RFP_Species6))

contrasts(dat$GFP_Strain) <-
  contr.sum(nlevels(dat$GFP_Strain))

cat("\nDATA SUMMARY\n")
cat("Date-level observations:", nrow(dat), "\n")
cat("RFP strains:", nlevels(dat$RFP_Strain), "\n")
cat("Dates:", nlevels(dat$Date), "\n")

print(with(dat, table(RFP_Species6, GFP_Strain)))

cat("\nDates per strain and competitor:\n")
print(replicate_counts)

# Used by Figure2.R.
analysis_data_day <- dat


# 2. Model specification

model_formula <- log_ratio ~
  RFP_Species6 * GFP_Strain +
  (1 | RFP_Strain) +
  (1 | Date)

# Compare four residual-variance structures.

dispersion_formulas <- list(
  homogeneous = ~ 1,
  by_competitor = ~ 0 + GFP_Strain,
  by_rfp_species = ~ 0 + RFP_Species6,
  additive = ~ GFP_Strain + RFP_Species6
)

optim_control <- glmmTMBControl(
  optCtrl = list(
    iter.max = 10000,
    eval.max = 10000
  )
)

fit_model <- function(dispformula, reml, data = dat) {
  glmmTMB(
    formula = model_formula,
    dispformula = dispformula,
    family = gaussian(),
    data = data,
    REML = reml,
    control = optim_control
  )
}

model_is_stable <- function(m) {
  !is.null(m) &&
    isTRUE(m$fit$convergence == 0L) &&
    isTRUE(m$sdr$pdHess) &&
    is.finite(as.numeric(logLik(m)))
}

model_stats <- function(model, name) {
  valid <- model_is_stable(model)
  
  data.frame(
    model = name,
    converged = if (is.null(model)) FALSE else
      isTRUE(model$fit$convergence == 0L),
    pdHess = if (is.null(model)) FALSE else
      isTRUE(model$sdr$pdHess),
    stable = valid,
    n_parameters = if (is.null(model)) NA_integer_ else
      as.integer(attr(logLik(model), "df")),
    AIC = if (is.null(model)) NA_real_ else AIC(model),
    BIC = if (is.null(model)) NA_real_ else BIC(model),
    logLik = if (is.null(model)) NA_real_ else
      as.numeric(logLik(model))
  )
}


# 3. Compare dispersion models (ML)

model_comparison <- NULL
ml_models <- NULL

if (CHECK_DISPERSION_MODELS) {
  
  ml_models <- lapply(
    names(dispersion_formulas),
    function(name) {
      tryCatch(
        fit_model(
          dispersion_formulas[[name]],
          reml = FALSE
        ),
        error = function(e) {
          warning(
            "ML fit failed for ",
            name, ": ",
            conditionMessage(e)
          )
          NULL
        }
      )
    }
  )
  
  names(ml_models) <- names(dispersion_formulas)
  
  model_comparison <- do.call(
    rbind,
    lapply(
      names(ml_models),
      function(name) {
        model_stats(ml_models[[name]], name)
      }
    )
  )
  
  model_comparison$delta_AIC <- NA_real_
  model_comparison$delta_BIC <- NA_real_
  
  eligible <- which(
    model_comparison$stable &
      is.finite(model_comparison$AIC) &
      is.finite(model_comparison$BIC)
  )
  
  if (!length(eligible)) {
    stop("No stable ML model available.")
  }
  
  model_comparison$delta_AIC[eligible] <-
    model_comparison$AIC[eligible] -
    min(model_comparison$AIC[eligible])
  
  model_comparison$delta_BIC[eligible] <-
    model_comparison$BIC[eligible] -
    min(model_comparison$BIC[eligible])
  
  model_comparison <- model_comparison[
    order(model_comparison$AIC),
    ,
    drop = FALSE
  ]
  
  rownames(model_comparison) <- NULL
  
  cat("\nML RESIDUAL-DISPERSION COMPARISON\n")
  
  print(
    model_comparison,
    row.names = FALSE,
    digits = 6
  )
  
  best <- model_comparison$model[
    which(model_comparison$stable)[1]
  ]
  
  if (!identical(best, "additive")) {
    warning(
      "Lowest-AIC dispersion model is '",
      best,
      "'. The fitted model below retains additive dispersion."
    )
  }
}


# 4. Final model (REML)

m_final <- fit_model(
  dispersion_formulas$additive,
  reml = TRUE
)

if (!model_is_stable(m_final)) {
  stop(
    "Final additive-variance REML model did not ",
    "converge with a positive-definite Hessian."
  )
}

if (nobs(m_final) != nrow(dat)) {
  stop("The final model dropped observations.")
}

cat("\nFINAL GAUSSIAN GLMMTMB MODEL (REML)\n")
print(summary(m_final))

cat("\nRANDOM EFFECTS\n")
print(VarCorr(m_final))

# Type III fixed-effect tests.

type3_tests <- NULL

if (requireNamespace("car", quietly = TRUE)) {
  type3_tests <- car::Anova(
    m_final,
    type = 3,
    component = "cond"
  )
  
  cat("\nTYPE III CONDITIONAL FIXED-EFFECT WALD TESTS\n")
  print(type3_tests)
  
} else {
  warning("Install 'car' for Type III fixed-effect tests.")
}


# 5. Model diagnostics

set.seed(DHARMA_SEED)

res_final <- DHARMa::simulateResiduals(
  fittedModel = m_final,
  n = DHARMA_SIMULATIONS,
  plot = FALSE
)

qc_tests <- data.frame(
  test = c(
    "DHARMa uniformity",
    "DHARMa dispersion",
    "DHARMa outliers"
  ),
  p_value = c(
    DHARMa::testUniformity(res_final)$p.value,
    DHARMa::testDispersion(res_final)$p.value,
    DHARMa::testOutliers(res_final)$p.value
  )
)

cat("\nFINAL-MODEL RESIDUAL DIAGNOSTICS\n")
print(qc_tests, row.names = FALSE)

if (RUN_DHARMA_PLOTS) {
  plot(res_final)
  
  DHARMa::plotResiduals(
    res_final,
    form = dat$GFP_Strain
  )
  
  DHARMa::plotResiduals(
    res_final,
    form = dat$RFP_Species6
  )
  
  DHARMa::plotResiduals(
    res_final,
    form = dat$Date
  )
}


# 6. Figure 2: marginal means

em_species <- emmeans::emmeans(
  m_final,
  ~ RFP_Species6 | GFP_Strain
)

species_em_order <- unique(
  as.character(
    as.data.frame(em_species)$RFP_Species6
  )
)

if (length(species_em_order) != 6L ||
    !setequal(species_em_order, species6)) {
  stop("Unexpected EMM species ordering or missing species.")
}

make_weights <- function(species_names) {
  w <- setNames(
    rep(0, length(species_em_order)),
    species_em_order
  )
  
  w[species_names] <- 1 / length(species_names)
  unname(w)
}

collapse_species <- list(
  "Non-AA" = make_weights(c("Pf", "Pu", "Ps")),
  "Ba" = make_weights("Ba"),
  "Bb" = make_weights("Bb"),
  "Bh" = make_weights("Bh")
)

em_group <- contrast(
  em_species,
  method = collapse_species,
  by = "GFP_Strain"
)

clean_group_names <- function(x) {
  z <- trimws(as.character(x))
  gsub("^\\((.*)\\)$", "\\1", z)
}

get_ci <- function(df) {
  lo <- intersect(
    c("lower.CL", "asymp.LCL", "LCL"),
    names(df)
  )[1]
  
  hi <- intersect(
    c("upper.CL", "asymp.UCL", "UCL"),
    names(df)
  )[1]
  
  if (is.na(lo) || is.na(hi)) {
    stop("EMM summary has no identifiable CI columns.")
  }
  
  data.frame(
    lower_ci = df[[lo]],
    upper_ci = df[[hi]]
  )
}

# Four groups x six competitors:
# Holm correction across 24 comparisons against zero.

cell_raw <- as.data.frame(summary(
  em_group,
  infer = c(TRUE, TRUE),
  null = 0,
  adjust = "none"
))

cell_raw$contrast <- clean_group_names(
  cell_raw$contrast
)

if (nrow(cell_raw) != 24L ||
    anyDuplicated(cell_raw[c("GFP_Strain", "contrast")]) ||
    !setequal(unique(cell_raw$contrast), plot_groups) ||
    !setequal(
      unique(as.character(cell_raw$GFP_Strain)),
      competitors
    ) ||
    any(!is.finite(cell_raw$p.value))) {
  stop("Expected 24 unique, finite EMM-vs-zero tests.")
}

emm_cell <- cbind(
  cell_raw,
  get_ci(cell_raw)
)

emm_cell$p_holm_24 <- p.adjust(
  emm_cell$p.value,
  method = "holm"
)

emm_cell$significant_vs_zero <-
  emm_cell$p_holm_24 < 0.05

emm_cell <- emm_cell[
  order(
    match(emm_cell$GFP_Strain, competitors),
    match(emm_cell$contrast, plot_groups)
  ),
  ,
  drop = FALSE
]

rownames(emm_cell) <- NULL

cat("\nFOUR-GROUP EMMs BY COMPETITOR: HOLM ACROSS 24 VS-ZERO TESTS\n")

print(
  emm_cell,
  row.names = FALSE,
  digits = 5
)


# 7. Figure 2: pairwise contrasts

em_group_order <- clean_group_names(
  unique(as.data.frame(em_group)$contrast)
)

if (length(em_group_order) != 4L ||
    !setequal(em_group_order, plot_groups)) {
  stop("Unexpected four-group EMM ordering.")
}

pair_spec <- list(
  "Ba - Bb" = c("Ba", "Bb"),
  "Ba - Bh" = c("Ba", "Bh"),
  "Bb - Bh" = c("Bb", "Bh"),
  "Ba - Non-AA" = c("Ba", "Non-AA"),
  "Bb - Non-AA" = c("Bb", "Non-AA"),
  "Bh - Non-AA" = c("Bh", "Non-AA")
)

pair_weights <- lapply(
  pair_spec,
  function(groups) {
    w <- setNames(rep(0, 4L), em_group_order)
    w[groups[[1]]] <- 1
    w[groups[[2]]] <- -1
    unname(w)
  }
)

pair_obj <- contrast(
  em_group,
  method = pair_weights,
  by = "GFP_Strain"
)

pair_raw <- as.data.frame(summary(
  pair_obj,
  infer = c(TRUE, TRUE),
  adjust = "none"
))

pair_raw$contrast <- gsub(
  "\\(Non-AA\\)",
  "Non-AA",
  pair_raw$contrast
)

if (nrow(pair_raw) != 36L ||
    anyDuplicated(pair_raw[c("GFP_Strain", "contrast")]) ||
    !setequal(unique(pair_raw$contrast), names(pair_spec)) ||
    any(!is.finite(pair_raw$p.value))) {
  stop("Expected 36 unique, finite pairwise tests.")
}

pairwise_by_comp <- cbind(
  pair_raw,
  get_ci(pair_raw)
) %>%
  group_by(GFP_Strain) %>%
  mutate(
    n_tests = n(),
    p_holm_6 = p.adjust(
      p.value,
      method = "holm"
    ),
    significant = p_holm_6 < 0.05,
    fold_difference = 10^estimate
  ) %>%
  ungroup() %>%
  arrange(GFP_Strain, p_holm_6)

if (any(pairwise_by_comp$n_tests != 6L)) {
  stop("At least one competitor has fewer than six tests.")
}

cat("\nPAIRWISE GROUP CONTRASTS: HOLM ACROSS SIX WITHIN EACH COMPETITOR\n")

print(
  as.data.frame(pairwise_by_comp),
  row.names = FALSE,
  digits = 5
)


# 8. Overall group means

# Average equally across six competitors.

em_species_overall <- emmeans(
  m_final,
  ~ RFP_Species6,
  weights = "equal"
)

order_overall <- as.character(
  as.data.frame(em_species_overall)$RFP_Species6
)

if (!setequal(order_overall, species6)) {
  stop("Missing species in overall EMMs.")
}

overall_weights <- lapply(
  list(
    "Non-AA" = c("Pf", "Pu", "Ps"),
    "Ba" = "Ba",
    "Bb" = "Bb",
    "Bh" = "Bh"
  ),
  function(members) {
    w <- setNames(rep(0, 6L), order_overall)
    w[members] <- 1 / length(members)
    unname(w)
  }
)

em_group_overall <- contrast(
  em_species_overall,
  method = overall_weights
)

overall_raw <- as.data.frame(summary(
  em_group_overall,
  infer = c(TRUE, TRUE),
  null = 0,
  adjust = "none"
))

overall_raw$contrast <- clean_group_names(
  overall_raw$contrast
)

if (nrow(overall_raw) != 4L ||
    any(!is.finite(overall_raw$p.value))) {
  stop("Expected four overall group means.")
}

overall_emm <- cbind(
  overall_raw,
  get_ci(overall_raw)
)

overall_emm$p_holm_4 <- p.adjust(
  overall_emm$p.value,
  method = "holm"
)

# Six overall pairwise comparisons.

overall_pairs_obj <- contrast(
  em_group_overall,
  method = pair_weights,
  adjust = "none"
)

overall_order <- clean_group_names(
  as.character(
    as.data.frame(em_group_overall)$contrast
  )
)

if (!identical(overall_order, em_group_order)) {
  stop("Overall and by-GFP EMM group orders differ.")
}

overall_pair_raw <- as.data.frame(summary(
  overall_pairs_obj,
  infer = c(TRUE, TRUE),
  adjust = "none"
))

overall_pair_raw$contrast <- gsub(
  "\\(Non-AA\\)",
  "Non-AA",
  overall_pair_raw$contrast
)

if (nrow(overall_pair_raw) != 6L ||
    any(!is.finite(overall_pair_raw$p.value))) {
  stop("Expected six overall pairwise comparisons.")
}

overall_pairwise <- cbind(
  overall_pair_raw,
  get_ci(overall_pair_raw)
)

overall_pairwise$p_holm_6 <- p.adjust(
  overall_pairwise$p.value,
  method = "holm"
)

cat("\nSECONDARY: OVERALL GROUP EMMs (EQUAL COMPETITOR WEIGHTS)\n")
print(overall_emm, row.names = FALSE, digits = 5)

cat("\nSECONDARY: OVERALL PAIRWISE GROUP CONTRASTS\n")
print(overall_pairwise, row.names = FALSE, digits = 5)


# 9. Figure 3: pooled contrasts

# Each RFP species and competitor receives equal weight.
# Positive estimates favour non-reduced species.

fig3_weights <- list(
  "4 non-reduced-genome species vs Bb/Bh" =
    c(
      Pf = 1/4, Pu = 1/4, Ps = 1/4,
      Ba = 1/4, Bb = -1/2, Bh = -1/2
    ),
  
  "3 non-AA species vs Bb/Bh" =
    c(
      Pf = 1/3, Pu = 1/3, Ps = 1/3,
      Ba = 0, Bb = -1/2, Bh = -1/2
    )
)

if (!all(vapply(
  fig3_weights,
  function(x) abs(sum(x)) < 1e-12,
  logical(1)
))) {
  stop("Figure 3 contrast coefficients do not sum to zero.")
}

fig3_grid_weights <- function(em, named_weights) {
  grid_species <- unique(
    as.character(as.data.frame(em)$RFP_Species6)
  )
  
  if (length(grid_species) != length(named_weights) ||
      !setequal(grid_species, names(named_weights))) {
    stop("Figure 3 weights do not match the EMM grid.")
  }
  
  unname(named_weights[grid_species])
}

fig3_calculate <- function(em, named_weights, label) {
  obj <- emmeans::contrast(
    em,
    method = setNames(
      list(fig3_grid_weights(em, named_weights)),
      label
    )
  )
  
  ans <- as.data.frame(summary(
    obj,
    infer = c(TRUE, TRUE),
    adjust = "none"
  ))
  
  ci <- get_ci(ans)
  
  ans$CI_lower <- ci$lower_ci
  ans$CI_upper <- ci$upper_ci
  ans$fold_difference <- 10^ans$estimate
  ans$fold_CI_lower <- 10^ans$CI_lower
  ans$fold_CI_upper <- 10^ans$CI_upper
  
  ans
}

pooled_em <- em_species_overall

pooled <- do.call(
  rbind,
  lapply(
    names(fig3_weights),
    function(label) {
      fig3_calculate(
        pooled_em,
        fig3_weights[[label]],
        label
      )
    }
  )
)

rownames(pooled) <- NULL

if (nrow(pooled) != 2L ||
    !identical(
      as.character(pooled$contrast),
      names(fig3_weights)
    ) ||
    any(!is.finite(pooled$estimate)) ||
    any(!is.finite(pooled$p.value))) {
  stop("Expected two primary Figure 3 contrasts.")
}

# Holm correction across two planned contrasts.

pooled$p_holm_2 <- stats::p.adjust(
  pooled$p.value,
  method = "holm"
)

cat("\nFIGURE 3: TWO PRIMARY POOLED CONTRASTS (HOLM ACROSS TWO)\n")

print(
  pooled[, c(
    "contrast",
    "estimate",
    "SE",
    "CI_lower",
    "CI_upper",
    "fold_difference",
    "p.value",
    "p_holm_2"
  )],
  row.names = FALSE,
  digits = 7
)


# 10. Figure 3: competitor-specific contrasts

by_em <- em_species

by_gfp <- do.call(
  rbind,
  lapply(
    names(fig3_weights),
    function(label) {
      fig3_calculate(
        by_em,
        fig3_weights[[label]],
        label
      )
    }
  )
)

rownames(by_gfp) <- NULL
by_gfp$GFP_Strain <- as.character(by_gfp$GFP_Strain)

if (nrow(by_gfp) != 12L ||
    anyDuplicated(by_gfp[c("contrast", "GFP_Strain")]) ||
    !setequal(by_gfp$GFP_Strain, competitors) ||
    any(!is.finite(by_gfp$p.value)) ||
    any(table(by_gfp$contrast) != 6L)) {
  stop("Expected two contrasts across six competitors.")
}

by_gfp$p_holm_12 <- stats::p.adjust(
  by_gfp$p.value,
  method = "holm"
)

fig3_joint_p <- stats::p.adjust(
  c(pooled$p.value, by_gfp$p.value),
  method = "holm"
)

pooled$p_holm_14 <- fig3_joint_p[
  seq_len(nrow(pooled))
]

by_gfp$p_holm_14 <- fig3_joint_p[
  nrow(pooled) + seq_len(nrow(by_gfp))
]

cat("\nFIGURE 3: COMPETITOR-SPECIFIC FOLLOW-UPS (HOLM ACROSS 12)\n")

print(
  by_gfp[, c(
    "contrast",
    "GFP_Strain",
    "estimate",
    "SE",
    "CI_lower",
    "CI_upper",
    "p.value",
    "p_holm_12"
  )],
  row.names = FALSE,
  digits = 5
)

# Wald tests for heterogeneity across competitors.

heterogeneity <- do.call(
  rbind,
  lapply(
    names(fig3_weights),
    function(label) {
      
      obj <- emmeans::contrast(
        by_em,
        method = setNames(
          list(
            fig3_grid_weights(
              by_em,
              fig3_weights[[label]]
            )
          ),
          label
        )
      )
      
      e <- as.data.frame(
        summary(obj, adjust = "none")
      )$estimate
      
      V <- as.matrix(vcov(obj))
      
      if (length(e) != 6L ||
          !identical(dim(V), c(6L, 6L)) ||
          any(!is.finite(e)) ||
          any(!is.finite(V))) {
        stop("Invalid Figure 3 heterogeneity covariance.")
      }
      
      H <- cbind(-1, diag(5))
      d <- drop(H %*% e)
      
      stat <- as.numeric(
        crossprod(
          d,
          solve(H %*% V %*% t(H), d)
        )
      )
      
      data.frame(
        contrast = label,
        chi_square = stat,
        df = 5L,
        p.value = stats::pchisq(
          stat,
          df = 5L,
          lower.tail = FALSE
        )
      )
    }
  )
)

rownames(heterogeneity) <- NULL

cat("\nFIGURE 3: CONTRAST HETEROGENEITY ACROSS SIX GFP COMPETITORS\n")

print(
  heterogeneity,
  row.names = FALSE,
  digits = 5
)


# 11. Figure 3: sensitivity analyses

leave_gfp <- NULL
fig3_no_ba_refit <- NULL
m_no_ba <- NULL

if (RUN_FIG3_SENSITIVITIES) {
  
  # Leave-one-competitor-out marginal averages.
  
  leave_gfp <- do.call(
    rbind,
    lapply(
      c("At", "Bp", "Ck", "Pk", "Sl", "Ss"),
      function(omit) {
        
        sub_em <- suppressMessages(
          emmeans::emmeans(
            m_final,
            ~ RFP_Species6,
            at = list(
              GFP_Strain = setdiff(competitors, omit)
            ),
            weights = "equal"
          )
        )
        
        ans <- do.call(
          rbind,
          lapply(
            names(fig3_weights),
            function(label) {
              fig3_calculate(
                sub_em,
                fig3_weights[[label]],
                label
              )
            }
          )
        )
        
        ans$omitted_GFP <- omit
        ans
      }
    )
  )
  
  rownames(leave_gfp) <- NULL
  
  if (nrow(leave_gfp) != 12L ||
      any(!is.finite(leave_gfp$estimate))) {
    stop("Figure 3 leave-one-competitor-out audit failed.")
  }
  
  cat("\nFIGURE 3: LEAVE-ONE-COMPETITOR-OUT AVERAGES\n")
  
  print(
    leave_gfp[, c(
      "omitted_GFP",
      "contrast",
      "estimate",
      "SE",
      "p.value"
    )],
    row.names = FALSE,
    digits = 5
  )
  
  if (RUN_NO_BA_REFIT) {
    
    # Fit the model without P. agricolaris.
    
    no_ba <- droplevels(
      dat[
        as.character(dat$RFP_Species6) != "Ba",
        ,
        drop = FALSE
      ]
    )
    
    contrasts(no_ba$RFP_Species6) <-
      contr.sum(nlevels(no_ba$RFP_Species6))
    
    contrasts(no_ba$GFP_Strain) <-
      contr.sum(nlevels(no_ba$GFP_Strain))
    
    m_no_ba <- tryCatch(
      fit_model(
        dispersion_formulas$additive,
        reml = TRUE,
        data = no_ba
      ),
      error = function(e) {
        warning(
          "Figure 3 no-Ba refit failed: ",
          conditionMessage(e)
        )
        NULL
      }
    )
    
    if (model_is_stable(m_no_ba)) {
      
      em_no_ba <- suppressMessages(
        emmeans::emmeans(
          m_no_ba,
          ~ RFP_Species6,
          weights = "equal"
        )
      )
      
      weights_without_ba <-
        fig3_weights[["3 non-AA species vs Bb/Bh"]]
      
      weights_without_ba <- weights_without_ba[
        names(weights_without_ba) != "Ba"
      ]
      
      fig3_no_ba_refit <- fig3_calculate(
        em_no_ba,
        weights_without_ba,
        "3 non-AA vs Bb/Bh: model refitted without Ba"
      )
      
      cat("\nFIGURE 3: NO-BA MODEL REFIT\n")
      
      print(
        fig3_no_ba_refit[, c(
          "contrast",
          "estimate",
          "SE",
          "CI_lower",
          "CI_upper",
          "p.value"
        )],
        row.names = FALSE,
        digits = 5
      )
      
    } else {
      warning("No stable Figure 3 model refit without Ba.")
    }
  }
}


# 12. Optional output files

if (SAVE_TABLES || SAVE_MODEL_RDS) {
  dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

if (SAVE_TABLES) {
  
  tables <- list(
    model_comparison = model_comparison,
    residual_diagnostics = qc_tests,
    EMMs_by_competitor = emm_cell,
    pairwise_within_competitor = as.data.frame(pairwise_by_comp),
    overall_EMMs = overall_emm,
    overall_pairwise = overall_pairwise,
    Fig3_primary_pooled = pooled,
    Fig3_by_competitor = by_gfp,
    Fig3_heterogeneity = heterogeneity,
    Fig3_leave_one_competitor_out = leave_gfp,
    Fig3_no_Ba_refit = fig3_no_ba_refit
  )
  
  for (name in names(tables)) {
    if (!is.null(tables[[name]])) {
      write.csv(
        as.data.frame(tables[[name]]),
        file.path(OUTPUT_DIR, paste0(name, ".csv")),
        row.names = FALSE
      )
    }
  }
  
  if (!is.null(type3_tests)) {
    write.csv(
      as.data.frame(type3_tests),
      file.path(OUTPUT_DIR, "type3_fixed_effects.csv")
    )
  }
  
  writeLines(
    capture.output(sessionInfo()),
    file.path(OUTPUT_DIR, "R_session_info.txt")
  )
}

if (SAVE_MODEL_RDS) {
  saveRDS(
    m_final,
    file.path(OUTPUT_DIR, "competition_model_final.rds")
  )
}

cat("\nAUC-ratio analysis complete.\n")
