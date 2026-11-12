
# Figure S2b: average change in AUC relative to monoculture
# Uses the final replicate-level AUC dataset from the competition import.
# RFP and GFP responses are fitted separately; they are not derived from
# the log-ratio competition model.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
  library(emmeans)
})

INPUT_CANDIDATES <- c(
  "analysis_data_rep_all_valid.rds",
  file.path("data", "processed", "analysis_data_rep_all_valid.rds")
)
SAVE_FIGS2B <- FALSE
OUTPUT_DIR <- file.path("outputs", "figures")

species_levels <- c("Pf", "Pu", "Ps", "Ba", "Bb", "Bh")
competitor_levels <- c("Sl", "Ss", "At", "Bp", "Ck", "Pk")

input_found <- INPUT_CANDIDATES[file.exists(INPUT_CANDIDATES)]
if (!length(input_found)) {
  stop("Missing analysis_data_rep_all_valid.rds. Run the final import first.")
}

auc_rep <- readRDS(input_found[[1]])
needed <- c(
  "Date", "RFP_Strain", "RFP_Species6", "GFP_Strain",
  "Replicate", "Normalized_AUC_RFP", "Normalized_AUC_GFP"
)
if (!all(needed %in% names(auc_rep))) {
  stop("The imported AUC data are missing required columns: ",
       paste(setdiff(needed, names(auc_rep)), collapse = ", "))
}

# Average log10(AUC) across technical wells within each experimental date.
auc_day <- auc_rep %>%
  ungroup() %>%
  filter(
    as.character(RFP_Species6) %in% species_levels,
    as.character(GFP_Strain) %in% competitor_levels,
    is.finite(Normalized_AUC_RFP),
    is.finite(Normalized_AUC_GFP),
    Normalized_AUC_RFP > 0, Normalized_AUC_GFP > 0
  ) %>%
  mutate(
    log_RFP = log10(Normalized_AUC_RFP),
    log_GFP = log10(Normalized_AUC_GFP)
  ) %>%
  group_by(Date, RFP_Strain, RFP_Species6, GFP_Strain) %>%
  summarise(
    log_RFP = mean(log_RFP),
    log_GFP = mean(log_GFP),
    .groups = "drop"
  ) %>%
  mutate(
    Date = factor(Date),
    RFP_Strain = factor(RFP_Strain),
    RFP_Species6 = factor(as.character(RFP_Species6), levels = species_levels),
    GFP_Strain = factor(as.character(GFP_Strain), levels = competitor_levels)
  )

if (nrow(auc_day) != 282L ||
    n_distinct(auc_day$RFP_Strain) != 15L ||
    nrow(distinct(auc_day, RFP_Strain, GFP_Strain)) != 90L ||
    any(!is.finite(auc_day$log_RFP)) ||
    any(!is.finite(auc_day$log_GFP))) {
  stop("Expected 282 complete date-level observations from 90 strain/competitor pairs.")
}

cat("\nFigure S2b: ", nrow(auc_day), " experimental dates\n", sep = "")

# Separate Gaussian mixed models for the two fluorescence channels.
fit_auc_model <- function(response) {
  model <- glmmTMB(
    formula = as.formula(paste0(
      response, " ~ RFP_Species6 * GFP_Strain + ",
      "(1 | RFP_Strain) + (1 | Date)"
    )),
    dispformula = ~ RFP_Species6 + GFP_Strain,
    family = gaussian(),
    data = auc_day,
    REML = TRUE
  )
  if (!isTRUE(model$fit$convergence == 0L) ||
      !isTRUE(model$sdr$pdHess) ||
      nobs(model) != nrow(auc_day)) {
    stop("Unstable fit or dropped observations for ", response)
  }
  model
}

model_RFP <- fit_auc_model("log_RFP")
model_GFP <- fit_auc_model("log_GFP")

# Average equally across six competitor species, then six RFP species.
extract_estimate <- function(model, response) {
  em_species <- emmeans(
    model, ~ RFP_Species6, weights = "equal"
  )
  actual_order <- as.character(as.data.frame(em_species)$RFP_Species6)
  if (!identical(actual_order, species_levels)) {
    stop("Unexpected species order in marginal means.")
  }
  
  em_all <- contrast(
    em_species,
    method = list(Grand_mean = rep(1 / length(species_levels),
                                   length(species_levels)))
  )
  result <- as.data.frame(summary(em_all, infer = FALSE, adjust = "none"))
  estimate <- result$estimate[[1]]
  se <- result$SE[[1]]
  lower <- estimate - qnorm(0.975) * se
  upper <- estimate + qnorm(0.975) * se
  
  data.frame(
    Response = response,
    Estimate_log10 = estimate,
    SE_log10 = se,
    Percent = 100 * (10^estimate - 1),
    Lower = 100 * (10^lower - 1),
    Upper = 100 * (10^upper - 1),
    p_raw = 2 * pnorm(-abs(estimate / se))
  )
}

plot_df <- bind_rows(
  extract_estimate(model_RFP, "Paraburkholderia"),
  extract_estimate(model_GFP, "Competitor")
) %>%
  mutate(
    p_holm = p.adjust(p_raw, method = "holm"),
    Response = factor(Response, levels = c("Paraburkholderia", "Competitor")),
    Effect_label = sprintf("%.1f%%", Percent),
    Label_y = Upper + 5
  )

cat("\nFIGURE S2b: MODEL-ESTIMATED AUC CHANGE\n")
print(as.data.frame(plot_df), row.names = FALSE, digits = 7)

plot_colours <- c(
  "Paraburkholderia" = "#AD5876",
  "Competitor" = "#408E83"
)

# Expand the display only if the updated intervals or labels need room.
y_lower <- min(-68, floor((min(c(plot_df$Lower, 0)) - 5) / 10) * 10)
y_upper <- max(8, ceiling((max(c(plot_df$Label_y, 0)) + 2) / 5) * 5)
y_breaks <- seq(ceiling(y_lower / 10) * 10,
                floor(y_upper / 10) * 10, by = 10)

p_s2b_pub <- ggplot(plot_df, aes(Response, Percent, colour = Response)) +
  geom_hline(
    yintercept = 0, linetype = "dashed", linewidth = 0.5,
    colour = "grey50"
  ) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.14, linewidth = 0.85) +
  geom_point(shape = 21, size = 4, stroke = 1, fill = "white") +
  geom_text(
    aes(y = Label_y, label = Effect_label),
    colour = "black", fontface = "bold", size = 4.1,
    vjust = 0, show.legend = FALSE
  ) +
  scale_x_discrete(
    labels = expression(
      atop(italic(Paraburkholderia), "(RFP)"),
      atop("Competitor", "(GFP)")
    ),
    expand = expansion(add = 0.55)
  ) +
  scale_colour_manual(values = plot_colours, guide = "none") +
  scale_y_continuous(
    name = "Change in AUC relative to monoculture (%)",
    breaks = y_breaks,
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(ylim = c(y_lower, y_upper), clip = "off") +
  labs(x = NULL) +
  theme_classic(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 11, colour = "black", margin = margin(t = 9)),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.y = element_line(colour = "black", linewidth = 0.45),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.4),
    legend.position = "none",
    plot.margin = margin(t = 10, r = 8, b = 5, l = 5)
  )

print(p_s2b_pub)

if (SAVE_FIGS2B) {
  if (!requireNamespace("svglite", quietly = TRUE)) {
    stop("Install 'svglite' to export the SVG figure.")
  }
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = file.path(OUTPUT_DIR, "Figure_S2b.svg"),
    plot = p_s2b_pub, device = svglite::svglite,
    width = 4, height = 4, units = "in", bg = "white"
  )
}
