
# Figure 2: relative Paraburkholderia competitive performance
# Run after 02_Competition_model_and_EMMs_final.R in the same R session.
# The model script supplies the date-level data, EMMs and Holm-adjusted tests.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Settings
SAVE_FIG2 <- FALSE
OUTPUT_DIR <- file.path("outputs", "figures")

plot_groups <- c("Non-AA", "Ba", "Bb", "Bh")
gfp_order <- c("At", "Bp", "Ck", "Pk", "Sl", "Ss")

colours <- c(
  "Non-AA" = "#C7C7C7",
  "Ba" = "#F0A202",
  "Bb" = "#6EC6F2",
  "Bh" = "#5CC8A1"
)

competitor_names <- c(
  "At" = "Agrobacterium", "Bp" = "Brucella",
  "Ck" = "Comamonas", "Pk" = "Pseudomonas",
  "Sl" = "Serratia", "Ss" = "Staphylococcus"
)

# Required model outputs
required_objects <- c("analysis_data_day", "m_final", "emm_cell", "pairwise_by_comp")
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]
if (length(missing_objects)) {
  stop("Run 02_Competition_model_and_EMMs_final.R first. Missing: ",
       paste(missing_objects, collapse = ", "))
}

if (!isTRUE(m_final$fit$convergence == 0L) ||
    !isTRUE(m_final$sdr$pdHess) ||
    nrow(analysis_data_day) != nobs(m_final)) {
  stop("The fitted model and date-level data must match.")
}

# Model means and the 24 tests against zero (one Holm family)
needed_emm <- c("GFP_Strain", "contrast", "estimate", "lower_ci",
                "upper_ci", "p_holm_24", "significant_vs_zero")
if (!all(needed_emm %in% names(emm_cell))) {
  stop("Missing EMM columns: ", paste(setdiff(needed_emm, names(emm_cell)), collapse = ", "))
}

fig2_emm <- emm_cell %>%
  mutate(
    GFP_Strain = factor(as.character(GFP_Strain), levels = gfp_order),
    Group = factor(as.character(contrast), levels = plot_groups),
    x = match(as.character(Group), plot_groups)
  )

if (nrow(fig2_emm) != 24L ||
    anyDuplicated(fig2_emm[c("GFP_Strain", "Group")]) ||
    any(table(fig2_emm$GFP_Strain) != 4L) ||
    anyNA(fig2_emm[c("GFP_Strain", "Group", "x")]) ||
    any(!is.finite(as.matrix(fig2_emm[c("estimate", "lower_ci", "upper_ci", "p_holm_24")]))) ||
    any(is.na(fig2_emm$significant_vs_zero))) {
  stop("Figure 2 requires 24 complete group-by-competitor EMMs.")
}

# Underlying date-level observations; Pf, Pu and Ps remain separate points
fig2_raw <- analysis_data_day %>%
  mutate(
    GFP_Strain = factor(as.character(GFP_Strain), levels = gfp_order),
    Group = case_when(
      as.character(RFP_Species6) %in% c("Pf", "Pu", "Ps") ~ "Non-AA",
      as.character(RFP_Species6) %in% c("Ba", "Bb", "Bh") ~
        as.character(RFP_Species6),
      TRUE ~ NA_character_
    ),
    x = match(Group, plot_groups)
  ) %>%
  filter(!is.na(GFP_Strain), !is.na(x), is.finite(log_ratio))

if (!nrow(fig2_raw)) stop("No date-level data for Figure 2.")

# Six pairwise contrasts per competitor, Holm-adjusted within competitor
needed_pairs <- c("GFP_Strain", "contrast", "estimate", "p_holm_6")
if (!all(needed_pairs %in% names(pairwise_by_comp))) {
  stop("Missing pairwise columns: ",
       paste(setdiff(needed_pairs, names(pairwise_by_comp)), collapse = ", "))
}

fig2_pairs <- pairwise_by_comp %>%
  mutate(
    GFP_Strain = factor(as.character(GFP_Strain), levels = gfp_order),
    contrast = gsub("\\(Non-AA\\)", "Non-AA", as.character(contrast))
  )

expected_pairs <- c(
  "Ba - Bb", "Ba - Bh", "Bb - Bh",
  "Ba - Non-AA", "Bb - Non-AA", "Bh - Non-AA"
)
if (nrow(fig2_pairs) != 36L ||
    anyDuplicated(fig2_pairs[c("GFP_Strain", "contrast")]) ||
    !setequal(unique(fig2_pairs$contrast), expected_pairs) ||
    any(table(fig2_pairs$GFP_Strain) != 6L) ||
    any(!is.finite(fig2_pairs$p_holm_6))) {
  stop("Figure 2 requires six complete Holm-adjusted pairwise tests per competitor.")
}

parts <- strsplit(fig2_pairs$contrast, " - ", fixed = TRUE)
if (any(lengths(parts) != 2L)) stop("Unrecognised pairwise contrast labels.")
fig2_pairs$group1 <- vapply(parts, `[`, character(1), 1L)
fig2_pairs$group2 <- vapply(parts, `[`, character(1), 2L)
fig2_pairs$x1 <- match(fig2_pairs$group1, plot_groups)
fig2_pairs$x2 <- match(fig2_pairs$group2, plot_groups)
if (anyNA(fig2_pairs[c("x1", "x2")])) stop("Unknown group in pairwise comparisons.")

fig2_pairs <- fig2_pairs %>%
  mutate(
    x_left = pmin(x1, x2),
    x_right = pmax(x1, x2),
    x_middle = (x_left + x_right) / 2,
    span = x_right - x_left,
    stars = case_when(
      p_holm_6 < 0.001 ~ "***",
      p_holm_6 < 0.01 ~ "**",
      p_holm_6 < 0.05 ~ "*",
      TRUE ~ NA_character_
    )
  )

# Place shorter brackets on lower tiers. Non-overlapping brackets can share a tier.
assign_tiers <- function(d) {
  d <- d[order(d$span, d$p_holm_6), , drop = FALSE]
  occupied <- list()
  d$tier <- integer(nrow(d))
  
  for (i in seq_len(nrow(d))) {
    tier <- 1L
    repeat {
      if (tier > length(occupied)) {
        occupied[[tier]] <- matrix(numeric(), ncol = 2L)
      }
      spans <- occupied[[tier]]
      conflict <- nrow(spans) > 0L &&
        any(d$x_left[i] <= spans[, 2] & d$x_right[i] >= spans[, 1])
      if (!conflict) break
      tier <- tier + 1L
    }
    occupied[[tier]] <- rbind(occupied[[tier]], c(d$x_left[i], d$x_right[i]))
    d$tier[i] <- tier
  }
  d
}

fig2_tops <- bind_rows(
  fig2_raw %>% transmute(GFP_Strain, top = log_ratio),
  fig2_emm %>% transmute(GFP_Strain, top = upper_ci)
) %>%
  group_by(GFP_Strain) %>%
  summarise(top = max(top), .groups = "drop")

fig2_sig <- filter(fig2_pairs, !is.na(stars))
if (nrow(fig2_sig)) {
  fig2_brackets <- fig2_sig %>%
    split(.$GFP_Strain, drop = TRUE) %>%
    lapply(assign_tiers) %>%
    bind_rows() %>%
    left_join(fig2_tops, by = "GFP_Strain") %>%
    mutate(
      y = top + 0.14 + (tier - 1L) * 0.14,
      star_y = y + 0.035
    )
} else {
  fig2_brackets <- fig2_sig %>%
    mutate(tier = integer(), top = numeric(),
           y = numeric(), star_y = numeric())
}

if (nrow(fig2_brackets) != sum(fig2_pairs$p_holm_6 < 0.05) ||
    any(!is.finite(fig2_brackets$y))) {
  stop("Significance bracket positions are incomplete.")
}

# Shared limits across panels. Keep the previous 1.30 upper limit where possible;
# expand only if observations, intervals or significance brackets need more room.
y_values <- c(
  fig2_raw$log_ratio, fig2_emm$lower_ci, fig2_emm$upper_ci,
  fig2_brackets$star_y, 0
)
y_values <- y_values[is.finite(y_values)]
y_min <- floor((min(y_values) - 0.13) / 0.25) * 0.25
y_max <- max(1.30, ceiling((max(y_values) + 0.09) * 10) / 10)
y_ticks <- seq(ceiling(y_min * 2) / 2, floor(y_max * 2) / 2, by = 0.5)

# Plot
p_fig2 <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.65) +
  geom_point(
    data = fig2_raw,
    aes(x = x, y = log_ratio),
    position = position_jitter(width = 0.11, height = 0, seed = 2026),
    size = 2.0, alpha = 0.42, colour = "grey55"
  ) +
  geom_errorbar(
    data = fig2_emm,
    aes(x = x, ymin = lower_ci, ymax = upper_ci),
    width = 0.16, linewidth = 1.0, colour = "black"
  ) +
  geom_point(
    data = fig2_emm,
    aes(x = x, y = estimate, fill = Group),
    shape = 23, size = 5.7, stroke = 1.0, colour = "black"
  ) +
  geom_point(
    data = filter(fig2_emm, significant_vs_zero),
    aes(x = x, y = estimate),
    shape = 23, size = 2.0, stroke = 0.35,
    fill = "black", colour = "black"
  ) +
  geom_segment(
    data = fig2_brackets,
    aes(x = x_left, xend = x_right, y = y, yend = y),
    linewidth = 0.48, inherit.aes = FALSE
  ) +
  geom_segment(
    data = fig2_brackets,
    aes(x = x_left, xend = x_left, y = y - 0.035, yend = y),
    linewidth = 0.48, inherit.aes = FALSE
  ) +
  geom_segment(
    data = fig2_brackets,
    aes(x = x_right, xend = x_right, y = y - 0.035, yend = y),
    linewidth = 0.48, inherit.aes = FALSE
  ) +
  geom_text(
    data = fig2_brackets,
    aes(x = x_middle, y = star_y, label = stars),
    size = 4.2, vjust = 0, inherit.aes = FALSE
  ) +
  facet_wrap(
    ~ GFP_Strain, ncol = 3, drop = FALSE,
    labeller = labeller(GFP_Strain = competitor_names)
  ) +
  scale_x_continuous(
    limits = c(0.48, 4.52), breaks = 1:4,
    labels = expression(
      "Non-AA", italic("P. agricolaris"),
      italic("P. bonniea"), italic("P. hayleyella")
    ),
    minor_breaks = NULL, expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    breaks = y_ticks, minor_breaks = NULL,
    expand = expansion(mult = 0)
  ) +
  scale_fill_manual(values = colours, guide = "none", drop = FALSE) +
  coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
  labs(
    x = NULL,
    y = expression(
      log[10] * "(relative " * italic(Paraburkholderia) *
        " competitive performance)"
    )
  ) +
  theme_classic(base_size = 15) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text.x = element_text(
      angle = 60, hjust = 1, vjust = 1,
      size = 12, colour = "grey20", margin = margin(t = 7)
    ),
    axis.text.y = element_text(size = 12, colour = "grey20"),
    axis.title.y = element_text(
      size = 15, colour = "black", margin = margin(r = 9)
    ),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 16, face = "bold", colour = "black",
      margin = margin(b = 7)
    ),
    panel.spacing.x = grid::unit(0.9, "lines"),
    panel.spacing.y = grid::unit(1.2, "lines"),
    legend.position = "none",
    plot.margin = margin(t = 12, r = 14, b = 18, l = 12)
  )

print(p_fig2)

cat("\nFigure 2: ", nrow(fig2_emm), " EMMs; ",
    nrow(fig2_pairs), " pairwise tests; ",
    nrow(fig2_brackets), " significant brackets.\n", sep = "")

if (SAVE_FIG2) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(OUTPUT_DIR, "Fig_2.png"), p_fig2,
         width = 14, height = 10, units = "in", dpi = 300, bg = "white")
  ggsave(file.path(OUTPUT_DIR, "Fig_2.svg"), p_fig2,
         device = svglite::svglite, fix_text_size = FALSE,
         width = 14, height = 10, units = "in", bg = "white")
}
