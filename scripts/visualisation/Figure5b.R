
# Figure 5b: supernatant experiment — visualization only
# Run after 5b_data_and_stats.R in the same R session.
# No raw-data import, model fitting, EMMs or p-value adjustments here.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

SAVE_FIGURE <- FALSE
OUTPUT_DIR <- file.path("outputs", "figures")

required_objects <- c("treatment_dat", "treatment_levels", "fig5b_plot_emm",
                      "fig5b_bracket_labels")
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]
if (length(missing_objects)) {
  stop("Run 5b_data_and_stats.R first. Missing: ",
       paste(missing_objects, collapse = ", "))
}

# Plot layout and aesthetics retained from the existing Figure 5b.
source_cols <- c("Pseudomonas control" = "#C7C7C7",
                 "P. agricolaris" = "#F0A202", "P. hayleyella" = "#5CC8A1")
context_shapes <- c("Monoculture" = 23, "Coculture" = 21)
treatment_layout <- tibble::tribble(
  ~Treatment,              ~Treatment_x, ~Axis_label,   ~Source,               ~Context,
  "Pseudomonas_mono",      1.00,          "Monoculture", "Pseudomonas control", "Monoculture",
  "Agricolaris_mono",      2.40,          "Monoculture", "P. agricolaris",     "Monoculture",
  "Agricolaris_coculture", 3.40,          "Coculture",   "P. agricolaris",     "Coculture",
  "Hayleyella_mono",       5.00,          "Monoculture", "P. hayleyella",      "Monoculture",
  "Hayleyella_coculture",  6.00,          "Coculture",   "P. hayleyella",      "Coculture"
) %>%
  mutate(Treatment = factor(Treatment, levels = treatment_levels),
         Source = factor(Source, levels = names(source_cols)),
         Context = factor(Context, levels = names(context_shapes)))

source_labels <- tibble::tribble(
  ~x, ~label,
  1.00, "Pseudomonas", 2.90, "P. agricolaris", 5.50, "P. hayleyella"
)
separators <- tibble::tibble(xintercept = c(1.7, 4.2))

plot_date <- treatment_dat %>%
  mutate(Treatment = factor(Treatment, levels = treatment_levels)) %>%
  left_join(treatment_layout, by = "Treatment")
plot_emm <- fig5b_plot_emm %>%
  mutate(Treatment = factor(Treatment, levels = treatment_levels)) %>%
  left_join(treatment_layout, by = "Treatment")
if (nrow(plot_emm) != 5L || anyNA(plot_date$Treatment_x) ||
    anyNA(plot_emm$Treatment_x)) {
  stop("Treatment plotting positions are incomplete.")
}

sig_brackets <- tibble::tibble(
  x_start = c(1.00, 2.40, 5.00),
  x_end = c(3.40, 3.40, 6.00),
  y = c(1.00, 0.82, 0.68),
  label = unname(fig5b_bracket_labels[
    c("baseline_Ba_coculture", "within_Ba", "within_Bh")
  ])
)
if (anyNA(sig_brackets$label)) stop("Incomplete significance labels.")

cat("\nFigure 5b brackets\n")
print(sig_brackets)

p_treatment_final <- ggplot(plot_date, aes(x = Treatment_x, y = log10_Fitness)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.72) +
  geom_vline(data = separators, aes(xintercept = xintercept),
             linetype = "dotted", colour = "grey55", linewidth = 0.95,
             inherit.aes = FALSE) +
  geom_point(position = position_jitter(width = 0.08, height = 0, seed = 1),
             shape = 16, size = 3, alpha = 0.44, colour = "grey55") +
  geom_errorbar(data = plot_emm,
                aes(x = Treatment_x, ymin = lower_ci, ymax = upper_ci),
                width = 0.13, linewidth = 1.05, colour = "black", inherit.aes = FALSE) +
  geom_point(data = plot_emm,
             aes(x = Treatment_x, y = emmean, fill = Source, shape = Context),
             size = 8.5, stroke = 1.35, colour = "black", inherit.aes = FALSE) +
  geom_segment(data = sig_brackets,
               aes(x = x_start, xend = x_end, y = y, yend = y),
               linewidth = 0.58, inherit.aes = FALSE) +
  geom_segment(data = sig_brackets,
               aes(x = x_start, xend = x_start, y = y, yend = y - 0.045),
               linewidth = 0.58, inherit.aes = FALSE) +
  geom_segment(data = sig_brackets,
               aes(x = x_end, xend = x_end, y = y, yend = y - 0.045),
               linewidth = 0.58, inherit.aes = FALSE) +
  geom_text(data = sig_brackets,
            aes(x = (x_start + x_end) / 2, y = y + 0.055, label = label),
            size = 5.4, fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = source_labels, aes(x = x, y = -1.92, label = label),
            size = 4.75, colour = "grey20", vjust = 0.5, inherit.aes = FALSE) +
  scale_x_continuous(breaks = treatment_layout$Treatment_x,
                     labels = treatment_layout$Axis_label,
                     limits = c(0.40, 6.55), expand = expansion(mult = 0)) +
  scale_y_continuous(breaks = seq(-1.5, 1.0, by = 0.5),
                     expand = expansion(mult = 0)) +
  scale_fill_manual(values = source_cols, guide = "none") +
  scale_shape_manual(values = context_shapes, guide = "none") +
  coord_cartesian(ylim = c(-1.60, 1.18), clip = "off") +
  labs(x = NULL, y = expression(log[10] * "(normalized Pseudomonas AUC)")) +
  theme_classic(base_size = 16) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.35),
    panel.grid = element_blank(), axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", linewidth = 0.55),
    axis.text.x = element_text(size = 13.5, colour = "grey20", margin = margin(t = 10)),
    axis.text.y = element_text(size = 13.5, colour = "grey20"),
    axis.title.y = element_text(size = 16, colour = "black", margin = margin(r = 8)),
    axis.title.x = element_blank(), legend.position = "none",
    plot.margin = margin(t = 12, r = 12, b = 58, l = 12)
  )

print(p_treatment_final)

if (SAVE_FIGURE) {
  if (!requireNamespace("svglite", quietly = TRUE)) stop("Install svglite for SVG output.")
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(OUTPUT_DIR, "Fig_5b.svg"), p_treatment_final,
         device = svglite::svglite, fix_text_size = FALSE,
         width = 8.5, height = 6.5, units = "in", bg = "white")
}
