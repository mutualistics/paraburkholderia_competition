
# Figure S1: normalized fluorescence trajectories
# Run after 1_FullDataImport_GitHub_no_auc_warnings.R in the same R session.
# Uses the final accepted AUC wells; no additional exclusions are applied.
# Curves are normalized to the same-day monoculture maximum, not AUC.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

SAVE_FIGS1 <- FALSE
OUTPUT_DIR <- file.path("outputs", "figures")
FIG_WIDTH <- 10
FIG_HEIGHT <- 8

Y_LOWER <- 0
Y_UPPER <- 1.25
PANEL_GAP_X <- 0.6
PANEL_GAP_Y <- 0.6

colours <- c(
  "Competitor (GFP)" = "#318978",
  "Paraburkholderia (RFP)" = "#B65E85"
)

competitors <- c("Pk", "Ck", "Bp", "Sl", "At", "Ss")
competitor_names <- c(
  Pk = "Pseudomonas", Ck = "Comamonas", Bp = "Brucella",
  Sl = "Serratia", At = "Agrobacterium", Ss = "Staphylococcus"
)
species_order <- c("Ba", "Bb", "Bh", "Pf", "Pu", "Ps")

strain_map <- tibble::tibble(
  RFP_Strain = c(
    "Ba80", "Ba161", "Ba31", "Ba159",
    "Bb862", "Bb822", "Bb859", "Bb433",
    "Bh11", "Bh171", "Bh22", "Bh21",
    "Pf", "Pu", "Ps"
  ),
  Species_plot = c(rep(species_order[1:3], each = 4), species_order[4:6])
)

if (!exists("full_data_corrected", inherits = TRUE) ||
    !exists("analysis_data_rep", inherits = TRUE)) {
  stop("Run the final competition import in this R session first.")
}

raw_columns <- c(
  "Date", "RFP_Strain", "GFP_Strain", "Replicate", "Well_ID",
  "Time", "Corrected_RFP", "Corrected_GFP"
)
rep_columns <- c(
  "Date", "RFP_Strain", "GFP_Strain", "Replicate", "RFP_Species6"
)

if (!all(raw_columns %in% names(full_data_corrected)) ||
    !all(rep_columns %in% names(analysis_data_rep))) {
  stop("Required corrected-fluorescence or accepted-replicate columns are missing.")
}

key_cols <- c("Date", "RFP_Strain", "GFP_Strain", "Replicate")

accepted <- analysis_data_rep %>%
  mutate(across(all_of(key_cols), as.character)) %>%
  filter(
    RFP_Strain %in% strain_map$RFP_Strain,
    GFP_Strain %in% competitors
  ) %>%
  distinct(across(all_of(key_cols)))

if (nrow(accepted) == 0L ||
    nrow(accepted %>% distinct(Date, RFP_Strain, GFP_Strain)) != 282L ||
    nrow(accepted %>% distinct(RFP_Strain, GFP_Strain)) != 90L) {
  stop("Expected 282 accepted experimental dates across 90 strain/competitor pairs.")
}

raw <- full_data_corrected %>%
  mutate(
    across(all_of(key_cols), as.character),
    Time = as.numeric(as.character(Time))
  )

# The original growth-curve denominator is the mean of the maximum
# corrected fluorescence for each same-day monoculture well.
max_controls <- raw %>%
  filter(GFP_Strain == "None" | RFP_Strain == "None") %>%
  filter(!(GFP_Strain == "Ss" & Time > 36 & Time <= 48)) %>%
  group_by(Date, RFP_Strain, GFP_Strain, Replicate) %>%
  summarise(
    max_RFP = if (any(is.finite(Corrected_RFP)))
      max(Corrected_RFP[is.finite(Corrected_RFP)]) else NA_real_,
    max_GFP = if (any(is.finite(Corrected_GFP)))
      max(Corrected_GFP[is.finite(Corrected_GFP)]) else NA_real_,
    .groups = "drop"
  ) %>%
  group_by(Date, RFP_Strain, GFP_Strain) %>%
  summarise(
    reference_RFP = mean(max_RFP, na.rm = TRUE),
    reference_GFP = mean(max_GFP, na.rm = TRUE),
    .groups = "drop"
  )

rfp_controls <- max_controls %>%
  filter(GFP_Strain == "None") %>%
  select(Date, RFP_Strain, reference_RFP)

gfp_controls <- max_controls %>%
  filter(RFP_Strain == "None") %>%
  select(Date, GFP_Strain, reference_GFP)

# Only wells retained by the final AUC import can enter Figure S1.
# Bh21 x Ck on 1307 is retained; failed runs are excluded upstream.
plot_data_wide <- raw %>%
  semi_join(accepted, by = key_cols) %>%
  inner_join(strain_map, by = "RFP_Strain") %>%
  filter(
    GFP_Strain %in% competitors,
    is.finite(Time), Time >= 0,
    if_else(GFP_Strain == "Ss", Time <= 36, Time <= 48)
  ) %>%
  left_join(rfp_controls, by = c("Date", "RFP_Strain")) %>%
  left_join(gfp_controls, by = c("Date", "GFP_Strain")) %>%
  mutate(
    Normalized_RFP = if_else(
      is.finite(reference_RFP) & reference_RFP > 0,
      Corrected_RFP / reference_RFP, NA_real_
    ),
    Normalized_GFP = if_else(
      is.finite(reference_GFP) & reference_GFP > 0,
      Corrected_GFP / reference_GFP, NA_real_
    ),
    Species_plot = factor(Species_plot, levels = species_order),
    Competitor_plot = factor(
      unname(competitor_names[GFP_Strain]),
      levels = unname(competitor_names[competitors])
    ),
    Trace_ID = interaction(
      Date, RFP_Strain, GFP_Strain, Replicate, Well_ID, drop = TRUE
    )
  ) %>%
  arrange(Competitor_plot, Species_plot, RFP_Strain, Date, Replicate, Time)

matched <- plot_data_wide %>% distinct(across(all_of(key_cols)))
missing_wells <- anti_join(accepted, matched, by = key_cols)
if (nrow(missing_wells) > 0L) {
  print(missing_wells, n = Inf)
  stop("Some accepted competition wells have no corrected growth-curve data.")
}

if (any(!is.finite(plot_data_wide$Normalized_RFP)) ||
    any(!is.finite(plot_data_wide$Normalized_GFP))) {
  stop("Missing or non-finite fluorescence normalization in accepted wells.")
}

# Date-level means are averaged within strain, never across strains.
daily_means <- plot_data_wide %>%
  group_by(Date, Time, GFP_Strain, RFP_Strain, Species_plot, Competitor_plot) %>%
  summarise(
    Mean_GFP = mean(Normalized_GFP),
    Mean_RFP = mean(Normalized_RFP),
    .groups = "drop"
  )

summary_curves <- daily_means %>%
  group_by(Time, GFP_Strain, RFP_Strain, Species_plot, Competitor_plot) %>%
  summarise(
    Final_Average_GFP = mean(Mean_GFP),
    Final_SE_GFP = sd(Mean_GFP) / sqrt(n()),
    Final_Average_RFP = mean(Mean_RFP),
    Final_SE_RFP = sd(Mean_RFP) / sqrt(n()),
    n_dates = n(),
    .groups = "drop"
  )

rep_counts <- plot_data_wide %>%
  group_by(GFP_Strain, RFP_Strain) %>%
  summarise(
    Num_Dates = n_distinct(Date),
    Num_Wells = n_distinct(Trace_ID),
    .groups = "drop"
  )

cat("\nFIGURE S1: ACCEPTED GROWTH CURVES\n")
cat("Experimental dates:", nrow(accepted %>% distinct(Date, RFP_Strain, GFP_Strain)), "\n")
cat("Competition wells:", n_distinct(plot_data_wide$Trace_ID), "\n")
cat("Strain/competitor pairs:", nrow(rep_counts), "\n")
cat("Bh21 x Ck, 1307 retained:", any(
  plot_data_wide$RFP_Strain == "Bh21" &
    plot_data_wide$GFP_Strain == "Ck" & plot_data_wide$Date == "1307"
), "\n")

trace_long <- plot_data_wide %>%
  pivot_longer(
    cols = c(Normalized_GFP, Normalized_RFP),
    names_to = "Channel", values_to = "Fluorescence"
  ) %>%
  mutate(Channel = recode(
    Channel,
    Normalized_GFP = "Competitor (GFP)",
    Normalized_RFP = "Paraburkholderia (RFP)"
  )) %>%
  arrange(factor(Channel, levels = names(colours)))

mean_long <- summary_curves %>%
  pivot_longer(
    cols = c(Final_Average_GFP, Final_Average_RFP),
    names_to = "Channel", values_to = "Fluorescence"
  ) %>%
  mutate(Channel = recode(
    Channel,
    Final_Average_GFP = "Competitor (GFP)",
    Final_Average_RFP = "Paraburkholderia (RFP)"
  )) %>%
  arrange(factor(Channel, levels = names(colours)))

se_long <- bind_rows(
  summary_curves %>% transmute(
    Time, Competitor_plot, Species_plot, RFP_Strain,
    Channel = "Competitor (GFP)",
    Lower = Final_Average_GFP - Final_SE_GFP,
    Upper = Final_Average_GFP + Final_SE_GFP
  ),
  summary_curves %>% transmute(
    Time, Competitor_plot, Species_plot, RFP_Strain,
    Channel = "Paraburkholderia (RFP)",
    Lower = Final_Average_RFP - Final_SE_RFP,
    Upper = Final_Average_RFP + Final_SE_RFP
  )
) %>%
  filter(Time %in% c(12, 24, 36, 48), is.finite(Lower), is.finite(Upper))

plot_theme <- theme_classic(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = "#D5D5D5", linewidth = 0.4),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "#555555", linewidth = 0.4),
    axis.ticks.length = grid::unit(2, "pt"),
    axis.text.x = element_text(size = 9, colour = "#222222", margin = margin(t = 6)),
    axis.text.y = element_text(size = 10, colour = "#222222", margin = margin(r = 6)),
    axis.title.x = element_text(size = 14, margin = margin(t = 15)),
    axis.title.y = element_text(size = 14, margin = margin(r = 13)),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.x = element_text(size = 14, face = "bold", margin = margin(b = 10)),
    strip.text.y.left = element_text(
      size = 11.5, face = "bold", angle = 0, hjust = 1, margin = margin(r = 10)
    ),
    panel.spacing.x = grid::unit(PANEL_GAP_X, "lines"),
    panel.spacing.y = grid::unit(PANEL_GAP_Y, "lines"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.text = element_text(size = 12),
    legend.key.width = grid::unit(1.2, "cm"),
    plot.margin = margin(t = 10, r = 16, b = 12, l = 12)
  )

make_base <- function() {
  ggplot() +
    geom_hline(yintercept = 1, colour = "#A6A6A6", linetype = "dashed", linewidth = 0.43) +
    facet_grid(
      rows = vars(Competitor_plot), cols = vars(Species_plot),
      switch = "y", drop = FALSE
    ) +
    scale_x_continuous(
      limits = c(0, 48), breaks = c(0, 24, 48), minor_breaks = NULL,
      expand = expansion(add = c(1.5, 1.5))
    ) +
    scale_y_continuous(
      breaks = c(0, 0.5, 1), labels = c("0.0", "0.5", "1.0"),
      minor_breaks = NULL, expand = expansion(mult = c(0, 0))
    ) +
    scale_colour_manual(
      values = colours, breaks = names(colours),
      labels = expression(
        "Competitor (GFP)", italic(Paraburkholderia) ~ "(RFP)"
      ),
      drop = FALSE
    ) +
    coord_cartesian(ylim = c(Y_LOWER, Y_UPPER), clip = "on") +
    labs(
      x = "Time (hours)",
      y = "Normalized fluorescence\nrelative to monoculture maximum",
      colour = NULL
    ) +
    plot_theme
}

all_panels_traces <- make_base() +
  geom_line(
    data = trace_long,
    aes(Time, Fluorescence, group = interaction(Trace_ID, Channel), colour = Channel),
    linewidth = 0.42, alpha = 0.30, na.rm = TRUE
  ) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, linewidth = 1.6)))

all_panels_summary <- make_base() +
  geom_line(
    data = mean_long,
    aes(Time, Fluorescence, group = interaction(RFP_Strain, Channel), colour = Channel),
    linewidth = 0.78, alpha = 0.88, na.rm = TRUE
  ) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, linewidth = 1.6)))

all_panels_summary_se <- all_panels_summary +
  geom_errorbar(
    data = se_long,
    aes(Time, ymin = Lower, ymax = Upper, colour = Channel, group = RFP_Strain),
    inherit.aes = FALSE, width = 0, linewidth = 0.45, alpha = 0.4
  )

clipping_audit <- plot_data_wide %>%
  summarise(
    GFP_above = sum(Normalized_GFP > Y_UPPER),
    RFP_above = sum(Normalized_RFP > Y_UPPER),
    max_GFP = max(Normalized_GFP),
    max_RFP = max(Normalized_RFP)
  )
cat("\nValues above displayed y-axis:\n")
print(clipping_audit)

print(all_panels_traces)
# print(all_panels_summary)
# print(all_panels_summary_se)

if (SAVE_FIGS1) {
  if (!requireNamespace("svglite", quietly = TRUE)) {
    stop("Install 'svglite' to export the SVG figure.")
  }
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  ggsave(
    filename = file.path(OUTPUT_DIR, "Fig_S1.svg"),
    plot = all_panels_traces, device = svglite::svglite,
    fix_text_size = FALSE,
    width = FIG_WIDTH, height = FIG_HEIGHT, units = "in", bg = "white"
  )
}
