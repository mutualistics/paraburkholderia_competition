
# Figure 5a: starting-ratio experiment — visualization only
# Run after 5a_data_summary.R in the same R session.
# No data import, normalization or statistical estimation occurs here.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

SAVE_FIGURE <- FALSE
OUTPUT_DIR <- file.path("outputs", "figures")

required_objects <- c("fig5a_mean", "error_data", "BA_LABELS",
                      "TIME_MAX", "TIME_SHIFT_HOURS")
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]
if (length(missing_objects)) {
  stop("Run 5a_data_summary.R first. Missing: ",
       paste(missing_objects, collapse = ", "))
}

BA_COLOURS <- c(
  "0%" = "#111111", "0.01%" = "#4B358F", "0.1%" = "#A42A63",
  "1%" = "#C7361B", "10%" = "#D98200", "25%" = "#8C8C8C",
  "40%" = "#B2B2B2", "50%" = "#C9C9C9"
)

all_y <- c(fig5a_mean$mean_Pk, error_data$ymin, error_data$ymax)
all_y <- all_y[is.finite(all_y)]
if (!length(all_y)) stop("No finite plotted values.")
y_limits <- c(min(0, floor((min(all_y) - 0.015) / 0.05) * 0.05),
              max(1.05, ceiling((max(all_y) + 0.015) / 0.05) * 0.05))
y_breaks <- seq(0, floor(y_limits[2] * 4) / 4, by = 0.25)
x_max <- TIME_MAX - TIME_SHIFT_HOURS

p_pk_curves_points <- ggplot() +
  geom_errorbar(data = error_data,
                aes(Time_plot, ymin = ymin, ymax = ymax, colour = Curve_group),
                width = 0.15, linewidth = 0.30, alpha = 0.50, na.rm = TRUE) +
  geom_line(data = filter(fig5a_mean, Curve_role == "high_agricolaris"),
            aes(Time_plot, mean_Pk, colour = Curve_group, group = Curve_group),
            linewidth = 1.20, alpha = 0.78, lineend = "round", na.rm = TRUE) +
  geom_line(data = filter(fig5a_mean, Curve_role != "high_agricolaris"),
            aes(Time_plot, mean_Pk, colour = Curve_group, group = Curve_group),
            linewidth = 1.20, alpha = 1, lineend = "round", na.rm = TRUE) +
  geom_point(data = fig5a_mean,
             aes(Time_plot, mean_Pk, colour = Curve_group),
             shape = 16, size = 2.40, na.rm = TRUE) +
  scale_colour_manual(values = BA_COLOURS, limits = BA_LABELS,
                      breaks = BA_LABELS, drop = FALSE,
                      name = "Starting P. agricolaris (%)") +
  guides(colour = guide_legend(ncol = 1, byrow = TRUE,
                               override.aes = list(linewidth = 1, shape = 16, size = 2.6, alpha = 1))) +
  scale_x_continuous(breaks = sort(unique(c(seq(0, x_max, by = 6), x_max))),
                     expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = y_breaks, expand = expansion(mult = 0)) +
  coord_cartesian(xlim = c(-0.6, x_max + 0.6), ylim = y_limits, expand = FALSE) +
  labs(x = "Time (hours)", y = "Normalized Pseudomonas fluorescence") +
  theme_classic(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.70),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", linewidth = 0.55),
    axis.text = element_text(colour = "grey20", size = 13.5),
    axis.title = element_text(colour = "black", size = 16),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    legend.position = "inside", legend.position.inside = c(0.040, 0.955),
    legend.justification = c(0, 1), legend.direction = "vertical",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11.5),
    legend.key.width = grid::unit(1.05, "lines"),
    legend.key.height = grid::unit(0.64, "lines"),
    legend.spacing.x = grid::unit(0.25, "lines"),
    legend.spacing.y = grid::unit(0.08, "lines"),
    legend.background = element_rect(fill = scales::alpha("white", 0.88), colour = NA),
    legend.margin = margin(4, 5, 4, 5),
    plot.margin = margin(12, 12, 12, 12)
  )

print(p_pk_curves_points)

if (SAVE_FIGURE) {
  if (!requireNamespace("svglite", quietly = TRUE)) stop("Install svglite for SVG output.")
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(OUTPUT_DIR, "Fig_5a.svg"), p_pk_curves_points,
         device = svglite::svglite, fix_text_size = FALSE,
         width = 8.5, height = 6.5, units = "in", bg = "white")
}
