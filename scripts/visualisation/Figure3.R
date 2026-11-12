# Figure 3: average competitive performance by genome group
# Run after 02_Competition_model_and_EMMs_final.R in the same R session.
# Uses the two pooled model contrasts and their Holm-adjusted p-values.

suppressPackageStartupMessages(library(ggplot2))

SAVE_FIG3 <- FALSE
OUTPUT_DIR <- file.path("outputs", "figures")

if (!exists("pooled", inherits = TRUE)) {
  stop("Run 02_Competition_model_and_EMMs_final.R first: 'pooled' is missing.")
}

needed <- c("contrast", "estimate", "CI_lower", "CI_upper", "p_holm_2")
if (!all(needed %in% names(pooled))) {
  stop("Missing pooled contrast columns: ",
       paste(setdiff(needed, names(pooled)), collapse = ", "))
}

expected <- c(
  "4 non-reduced-genome species vs Bb/Bh",
  "3 non-AA species vs Bb/Bh"
)
if (nrow(pooled) != 2L ||
    !identical(as.character(pooled$contrast), expected) ||
    any(!is.finite(as.matrix(pooled[c("estimate", "CI_lower", "CI_upper", "p_holm_2")]))) ||
    any(pooled$CI_lower > pooled$estimate) ||
    any(pooled$CI_upper < pooled$estimate)) {
  stop("Figure 3 requires the two complete pooled contrasts in model order.")
}

fig3_data <- pooled
fig3_data$comparison_label <- c(
  "All non-reduced-genome species\nvs reduced-genome symbionts",
  "Non-AA species\nvs reduced-genome symbionts"
)
fig3_data$comparison_key <- c("Non-reduced", "Non-AA")
fig3_data$y <- c(2, 1)
fig3_data$stars <- ifelse(
  fig3_data$p_holm_2 < 0.001, "***",
  ifelse(fig3_data$p_holm_2 < 0.01, "**",
         ifelse(fig3_data$p_holm_2 < 0.05, "*", ""))
)
fig3_data$cap_low <- fig3_data$y - 0.07
fig3_data$cap_high <- fig3_data$y + 0.07
fig3_data$star_x <- fig3_data$CI_upper + c(0.014, 0.012)

fill_values <- c("Non-reduced" = "#E4B41B", "Non-AA" = "#BEBEBE")

p_fig3 <- ggplot(fig3_data, aes(y = y)) +
  geom_vline(
    xintercept = 0, colour = "black", linetype = "dashed", linewidth = 0.8
  ) +
  geom_segment(
    aes(x = CI_lower, xend = CI_upper, yend = y),
    linewidth = 1.55, colour = "black"
  ) +
  geom_segment(
    aes(x = CI_lower, xend = CI_lower, y = cap_low, yend = cap_high),
    linewidth = 1.55, colour = "black"
  ) +
  geom_segment(
    aes(x = CI_upper, xend = CI_upper, y = cap_low, yend = cap_high),
    linewidth = 1.55, colour = "black"
  ) +
  geom_point(
    aes(x = estimate, fill = comparison_key),
    shape = 23, size = 8.2, stroke = 1.5, colour = "black"
  ) +
  geom_text(
    aes(x = star_x, label = stars),
    hjust = 0, vjust = 0.5, fontface = "bold", size = 8.3, colour = "black"
  ) +
  scale_fill_manual(values = fill_values, guide = "none") +
  scale_y_continuous(
    breaks = fig3_data$y, labels = fig3_data$comparison_label,
    limits = c(0.45, 2.55), expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(-0.06, 0.255), breaks = seq(-0.05, 0.25, by = 0.05),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = expression(
      Delta ~ log[10] * "(relative " * italic(Paraburkholderia) *
        " competitive performance)"
    ),
    y = NULL
  ) +
  theme_classic(base_size = 20) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.9),
    axis.line = element_blank(),
    axis.ticks.x = element_line(colour = "black", linewidth = 0.85),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(
      size = 18, colour = "black", lineheight = 1.03,
      margin = margin(r = 14)
    ),
    axis.title.x = element_text(size = 22, colour = "black", margin = margin(t = 12)),
    plot.margin = margin(t = 18, r = 44, b = 14, l = 18)
  )

print(p_fig3)
cat("\nFigure 3: two pooled contrasts, Holm-adjusted across two tests.\n")

if (SAVE_FIG3) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(OUTPUT_DIR, "Fig_3.png"), p_fig3,
         width = 10.2, height = 6.0, dpi = 600, bg = "white")
  ggsave(file.path(OUTPUT_DIR, "Fig_3.pdf"), p_fig3,
         width = 10.2, height = 6.0, bg = "white")
  ggsave(file.path(OUTPUT_DIR, "Fig_3.svg"), p_fig3,
         device = svglite::svglite, fix_text_size = FALSE,
         width = 10.2, height = 6.0, bg = "white")
}
