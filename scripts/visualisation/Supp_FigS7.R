
# Supplementary figure: both species across all 15 starting ratios
# Plotting only. Run the corrected 5a_data_summary.R first.
# Technical replicates and cross-date uncertainty are calculated upstream.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

SAVE_FIGURE <- FALSE
OUTPUT_DIR <- file.path("outputs", "figures")
SHOW_DATE_CURVES <- TRUE
FIGURE_WIDTH <- 11.5
FIGURE_HEIGHT <- 12.5

required <- c(
  "fig5a_all_day", "fig5a_all_mean", "fig5a_ratio_key",
  "TIME_MAX", "TIME_SHIFT_HOURS", "UNCERTAINTY"
)

missing <- required[
  !vapply(required, exists, logical(1), inherits = TRUE)
]

if (length(missing)) {
  stop(
    "Run the corrected 5a_data_summary.R first. Missing: ",
    paste(missing, collapse = ", ")
  )
}

need_day <- c(
  "Ratio_label", "Species", "Date", "Time_plot", "day_mean"
)

need_mean <- c(
  "Ratio_label", "Species", "Time_plot",
  "mean_fluorescence", "ymin", "ymax"
)

if (!all(need_day %in% names(fig5a_all_day)) ||
    !all(need_mean %in% names(fig5a_all_mean))) {
  stop("Expected the updated two-channel starting-ratio summary tables.")
}

if (nrow(fig5a_ratio_key) != 15L ||
    n_distinct(fig5a_all_day$Ratio_label) != 15L ||
    n_distinct(fig5a_all_mean$Ratio_label) != 15L) {
  stop("Expected exactly 15 starting ratios.")
}

# Match the fluorescence-channel palette of Figure S1.

channel_levels <- c("Pseudomonas", "P. agricolaris")

channel_colours <- c(
  "Pseudomonas" = "#318978",
  "P. agricolaris" = "#B65E85"
)

channel_labels <- expression(
  italic("Pseudomonas") ~ "(GFP)",
  italic("P. agricolaris") ~ "(RFP)"
)

ratio_order <- as.character(fig5a_ratio_key$Ratio_label)

plot_day <- fig5a_all_day %>%
  mutate(
    Ratio_label = factor(
      as.character(Ratio_label),
      levels = ratio_order
    ),
    Species = factor(
      as.character(Species),
      levels = channel_levels
    )
  ) %>%
  arrange(Ratio_label, Species, Date, Time_plot)

plot_mean <- fig5a_all_mean %>%
  mutate(
    Ratio_label = factor(
      as.character(Ratio_label),
      levels = ratio_order
    ),
    Species = factor(
      as.character(Species),
      levels = channel_levels
    )
  ) %>%
  arrange(Ratio_label, Species, Time_plot)

if (anyNA(plot_day$Species) ||
    anyNA(plot_mean$Species) ||
    anyNA(plot_day$Ratio_label) ||
    anyNA(plot_mean$Ratio_label)) {
  stop("Unexpected channel or starting-ratio labels.")
}

# Sparse uncertainty bars at four time positions.
# SE/CI values are supplied by the processing script.

x_end <- TIME_MAX - TIME_SHIFT_HOURS

if (!is.finite(x_end) || x_end <= 0) {
  stop("Invalid time window.")
}

bar_times <- c(6, 12, 18, x_end)

bar_times <- unique(
  bar_times[bar_times >= 0 & bar_times <= x_end]
)

error_bars <- plot_mean %>%
  filter(
    is.finite(ymin),
    is.finite(ymax)
  ) %>%
  tidyr::crossing(target_time = bar_times) %>%
  group_by(Ratio_label, Species, target_time) %>%
  slice_min(
    abs(Time_plot - target_time),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  distinct(
    Ratio_label,
    Species,
    Time_plot,
    .keep_all = TRUE
  )

# Identical axes across all panels, accommodating every observation.

y_values <- c(
  plot_day$day_mean,
  plot_mean$mean_fluorescence,
  plot_mean$ymin,
  plot_mean$ymax,
  1
)

y_values <- y_values[is.finite(y_values)]

if (!length(y_values)) {
  stop("No finite fluorescence measurements to plot.")
}

y_upper <- max(
  1.25,
  ceiling((max(y_values) + 0.03) * 4) / 4
)

y_breaks <- seq(0, y_upper, by = 0.5)

# Reference: same-day monoculture maximum = 1.

p_starting_ratios_all <- ggplot() +
  
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    colour = "#A6A6A6",
    linewidth = 0.43
  )

# Individual biological-date means.

if (SHOW_DATE_CURVES) {
  
  p_starting_ratios_all <- p_starting_ratios_all +
    
    geom_line(
      data = plot_day,
      aes(
        x = Time_plot,
        y = day_mean,
        colour = Species,
        group = interaction(Date, Species)
      ),
      linewidth = 0.42,
      alpha = 0.27,
      lineend = "round",
      na.rm = TRUE
    )
}

# Across-date means and uncertainty.

p_starting_ratios_all <- p_starting_ratios_all +
  
  geom_line(
    data = plot_mean,
    aes(
      x = Time_plot,
      y = mean_fluorescence,
      colour = Species,
      group = Species
    ),
    linewidth = 0.92,
    alpha = 0.95,
    lineend = "round",
    na.rm = TRUE
  ) +
  
  geom_errorbar(
    data = error_bars,
    aes(
      x = Time_plot,
      ymin = ymin,
      ymax = ymax,
      colour = Species,
      group = Species
    ),
    width = 0,
    linewidth = 0.48,
    alpha = 0.65,
    na.rm = TRUE
  ) +
  
  facet_wrap(
    ~ Ratio_label,
    ncol = 3,
    scales = "fixed",
    drop = FALSE
  ) +
  
  scale_colour_manual(
    name = NULL,
    values = channel_colours,
    breaks = channel_levels,
    labels = channel_labels,
    drop = FALSE
  ) +
  
  scale_x_continuous(
    breaks = sort(
      unique(c(seq(0, x_end, by = 6), x_end))
    ),
    minor_breaks = NULL,
    expand = expansion(add = c(0.4, 0.4))
  ) +
  
  scale_y_continuous(
    breaks = y_breaks,
    minor_breaks = NULL,
    expand = expansion(mult = c(0, 0))
  ) +
  
  coord_cartesian(
    xlim = c(0, x_end),
    ylim = c(0, y_upper),
    clip = "on"
  ) +
  
  labs(
    title = expression(
      "Initial " * italic("P. agricolaris") * " (%)"
    ),
    x = "Time (hours)",
    y = "Normalized fluorescence\nrelative to monoculture maximum"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    panel.grid = element_blank(),
    
    panel.border = element_rect(
      fill = NA,
      colour = "#D5D5D5",
      linewidth = 0.4
    ),
    
    axis.line = element_blank(),
    
    axis.ticks = element_line(
      colour = "#555555",
      linewidth = 0.4
    ),
    
    axis.ticks.length = grid::unit(2, "pt"),
    
    axis.text.x = element_text(
      size = 10,
      colour = "#222222",
      margin = margin(t = 6)
    ),
    
    axis.text.y = element_text(
      size = 10,
      colour = "#222222",
      margin = margin(r = 6)
    ),
    
    axis.title.x = element_text(
      size = 14,
      margin = margin(t = 12)
    ),
    
    axis.title.y = element_text(
      size = 14,
      margin = margin(r = 12)
    ),
    
    strip.background = element_blank(),
    strip.placement = "outside",
    
    strip.text.x = element_text(
      face = "bold",
      size = 12,
      colour = "#222222",
      margin = margin(b = 7)
    ),
    
    panel.spacing.x = grid::unit(0.65, "lines"),
    panel.spacing.y = grid::unit(0.75, "lines"),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    
    legend.text = element_text(size = 12),
    legend.key.width = grid::unit(1.2, "cm"),
    
    plot.title = element_text(
      size = 14,
      margin = margin(b = 12)
    ),
    
    plot.margin = margin(
      t = 10,
      r = 14,
      b = 12,
      l = 12
    )
  ) +
  
  guides(
    colour = guide_legend(
      nrow = 1,
      override.aes = list(
        alpha = 1,
        linewidth = 1.6
      )
    )
  )

print(p_starting_ratios_all)

cat("\nSupplementary starting-ratio figure:\n")

cat(
  "Starting ratios:",
  n_distinct(plot_mean$Ratio_label),
  "\n"
)

cat(
  "Species:",
  paste(channel_levels, collapse = ", "),
  "\n"
)

cat(
  "Biological dates:",
  n_distinct(plot_day$Date),
  "\n"
)

cat(
  "Uncertainty:",
  UNCERTAINTY,
  "(shown at four time positions)\n"
)

# Optional publication-quality exports.

if (SAVE_FIGURE) {
  
  if (!requireNamespace("svglite", quietly = TRUE)) {
    stop("Install 'svglite' to export the SVG figure.")
  }
  
  dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  ggsave(
    file.path(OUTPUT_DIR, "FigureS_starting_ratios.svg"),
    p_starting_ratios_all,
    device = svglite::svglite,
    width = FIGURE_WIDTH,
    height = FIGURE_HEIGHT,
    units = "in",
    bg = "white"
  )
  
  ggsave(
    file.path(OUTPUT_DIR, "FigureS_starting_ratios.png"),
    p_starting_ratios_all,
    width = FIGURE_WIDTH,
    height = FIGURE_HEIGHT,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}
