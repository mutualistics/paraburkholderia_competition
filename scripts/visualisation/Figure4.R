
# Figure 4: overall K and maximum early growth rate
# Run after 02_K_umax_models_and_EMMs_282.R in the same R session.
# EMMs, tests and observed means are produced by the statistics script.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
})

SAVE_FIGURE4 <- FALSE
OUTPUT_DIR <- file.path("outputs", "figures")

required <- c("kr_overall_emm", "kr_overall_pairs", "kr_raw_overall")
missing <- required[!vapply(required, exists, logical(1), inherits = TRUE)]

if (length(missing)) {
  stop(
    "Run 02_K_umax_models_and_EMMs_282.R first: ",
    paste(missing, collapse = ", ")
  )
}

plot_group_levels <- c("Non-AA", "Ba", "Bb", "Bh")
response_levels <- c("Paraburkholderia response", "Competitor response")
trait_levels <- c("K", "umax")
comp_levels <- c("Sl", "Ss", "At", "Bp", "Ck", "Pk")

BRACKET_BASE_PAD <- 0.35
BRACKET_STEP <- 0.32

group_colors <- c(
  "Non-AA" = "#C7C7C7",
  "Ba" = "#F0A202",
  "Bb" = "#6EC6F2",
  "Bh" = "#5CC8A1"
)

group_axis_labels <- expression(
  "Non-AA",
  italic("P. agricolaris"),
  italic("P. bonniea"),
  italic("P. hayleyella")
)

competitor_shapes <- c(
  Sl = 21, Ss = 22, At = 8,
  Bp = 24, Ck = 25, Pk = 23
)

# Plot data

raw_plot <- kr_raw_overall %>%
  mutate(
    Trait = factor(as.character(Trait), levels = trait_levels),
    Response = factor(as.character(Response), levels = response_levels),
    Plot_Group = factor(as.character(Plot_Group), levels = plot_group_levels),
    Competitor = factor(as.character(Competitor), levels = comp_levels)
  )

emm_plot <- kr_overall_emm %>%
  mutate(
    Trait = factor(Trait, levels = trait_levels),
    Response = factor(Response, levels = response_levels),
    Plot_Group = factor(Plot_Group, levels = plot_group_levels)
  )

pairwise_all <- kr_overall_pairs %>%
  mutate(
    Trait = factor(Trait, levels = trait_levels),
    Response = factor(Response, levels = response_levels),
    group1 = factor(group1, levels = plot_group_levels),
    group2 = factor(group2, levels = plot_group_levels)
  )

sig_tab <- filter(pairwise_all, !is.na(label))

# Four EMMs and six pairwise tests per response and trait.

if (nrow(emm_plot) != 16L ||
    nrow(pairwise_all) != 24L ||
    anyNA(emm_plot[c("Trait", "Response", "Plot_Group")]) ||
    anyNA(pairwise_all[c("Trait", "Response", "group1", "group2")]) ||
    any(!is.finite(emm_plot$emmean_log2)) ||
    any(!is.finite(pairwise_all$p_holm)) ||
    anyDuplicated(emm_plot[c("Trait", "Response", "Plot_Group")])) {
  stop("Figure 4 requires four EMMs and six contrasts per response and trait.")
}

# Plot one trait, with two response panels.

make_trait_plot <- function(trait_keep, y_lab, show_x_axis = TRUE) {
  
  raw_sub <- raw_plot %>%
    filter(Trait == trait_keep)
  
  emm_sub <- emm_plot %>%
    filter(Trait == trait_keep)
  
  sig_sub <- sig_tab %>%
    filter(Trait == trait_keep)
  
  # Place brackets independently in each response panel.
  
  facet_tops <- bind_rows(
    raw_sub %>% transmute(Response, y_top = y_log2),
    emm_sub %>% transmute(Response, y_top = upper_log2)
  ) %>%
    group_by(Response) %>%
    summarise(
      y_top = max(y_top, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (nrow(sig_sub) > 0) {
    
    if (anyNA(sig_sub$group1) || anyNA(sig_sub$group2)) {
      stop("Significant Figure 4 contrast has missing group coordinates.")
    }
    
    sig_sub <- sig_sub %>%
      left_join(facet_tops, by = "Response") %>%
      group_by(Response) %>%
      arrange(p_holm, .by_group = TRUE) %>%
      mutate(
        y.position = y_top + BRACKET_BASE_PAD +
          BRACKET_STEP * (row_number() - 1L)
      ) %>%
      ungroup()
    
  } else {
    sig_sub$y.position <- numeric(0)
  }
  
  if (nrow(sig_sub) > 0 &&
      any(!is.finite(sig_sub$y.position))) {
    stop("A Figure 4 significance bracket has a non-finite height.")
  }
  
  # Shared symmetric y-axis for both responses.
  
  all_values <- c(
    raw_sub$y_log2,
    emm_sub$lower_log2,
    emm_sub$upper_log2,
    sig_sub$y.position + 0.25,
    0
  )
  
  all_values <- all_values[is.finite(all_values)]
  
  if (!length(all_values)) {
    stop("No finite plotting values for ", trait_keep)
  }
  
  ymax_abs <- max(
    1,
    ceiling(max(abs(all_values)) * 2) / 2
  )
  
  y_limits_use <- c(-ymax_abs, ymax_abs)
  
  y_breaks <- seq(
    ceiling(y_limits_use[1]),
    floor(y_limits_use[2]),
    by = 1
  )
  
  p_base <- ggplot() +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.8,
      colour = "black"
    ) +
    
    # Observed competitor-specific means.
    
    geom_point(
      data = raw_sub,
      aes(
        x = Plot_Group,
        y = y_log2,
        shape = Competitor
      ),
      position = position_jitter(
        width = 0.10,
        height = 0,
        seed = 1
      ),
      size = 4.0,
      alpha = 0.62,
      stroke = 0.55,
      colour = "grey40",
      fill = "grey75"
    ) +
    
    # Pointwise 95% model confidence intervals.
    
    geom_errorbar(
      data = emm_sub,
      aes(
        x = Plot_Group,
        ymin = lower_log2,
        ymax = upper_log2
      ),
      width = 0.16,
      linewidth = 1.35,
      colour = "black"
    ) +
    
    # Model-adjusted EMM diamonds.
    
    geom_point(
      data = emm_sub,
      aes(
        x = Plot_Group,
        y = emmean_log2,
        fill = Plot_Group
      ),
      shape = 23,
      size = 8,
      stroke = 1.2,
      colour = "black"
    ) +
    
    # Black centres indicate Holm-significant differences from WT.
    
    geom_point(
      data = emm_sub %>% filter(sig_vs_zero),
      aes(
        x = Plot_Group,
        y = emmean_log2
      ),
      shape = 23,
      size = 3,
      stroke = 0.7,
      fill = "black",
      colour = "black"
    ) +
    
    facet_grid(
      . ~ Response,
      drop = FALSE
    ) +
    
    scale_fill_manual(
      values = group_colors,
      guide = "none",
      drop = FALSE
    ) +
    
    scale_shape_manual(
      values = competitor_shapes,
      limits = comp_levels,
      breaks = comp_levels,
      name = "Competitor",
      drop = FALSE,
      na.translate = FALSE
    ) +
    
    scale_x_discrete(
      limits = plot_group_levels,
      labels = group_axis_labels,
      drop = FALSE,
      na.translate = FALSE
    ) +
    
    scale_y_continuous(
      breaks = y_breaks,
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    
    coord_cartesian(
      ylim = y_limits_use,
      clip = "off"
    ) +
    
    labs(
      x = NULL,
      y = y_lab
    ) +
    
    theme_classic(base_size = 20) +
    
    theme(
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.7
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.7
      ),
      axis.text.x = if (show_x_axis) {
        element_text(
          angle = 60,
          hjust = 1,
          vjust = 1,
          size = 16,
          colour = "grey20",
          margin = margin(t = 8)
        )
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x_axis) {
        element_line(
          colour = "black",
          linewidth = 0.7
        )
      } else {
        element_blank()
      },
      axis.title.x = element_blank(),
      axis.text.y = element_text(
        size = 16,
        colour = "grey20"
      ),
      axis.title.y = element_text(
        size = 18,
        colour = "black"
      ),
      strip.background = element_blank(),
      strip.text.x = element_text(
        size = 18,
        face = "bold",
        colour = "black",
        margin = margin(b = 6)
      ),
      legend.position = "right",
      legend.title = element_text(
        size = 18,
        face = "bold"
      ),
      legend.text = element_text(size = 16),
      legend.key = element_blank(),
      plot.margin = margin(
        t = 12, r = 12, b = 18, l = 12
      )
    )
  
  if (nrow(sig_sub) > 0) {
    
    p_base <- p_base +
      ggpubr::stat_pvalue_manual(
        data = sig_sub,
        label = "label",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0.01,
        bracket.size = 0.55,
        size = 5.5,
        inherit.aes = FALSE
      )
  }
  
  p_base
}

# Main figure

p_k <- make_trait_plot(
  trait_keep = "K",
  y_lab = expression(log[2] ~ "fold-change in K vs WT"),
  show_x_axis = FALSE
)

p_umax <- make_trait_plot(
  trait_keep = "umax",
  y_lab = expression(
    log[2] ~ "fold-change in" ~ mu[max] ~ "vs WT"
  ),
  show_x_axis = TRUE
)

p_overall <- p_k / p_umax +
  plot_layout(
    heights = c(1, 1),
    guides = "collect"
  ) &
  theme(legend.position = "right")

print(p_overall)

cat(
  "\nFigure 4: ",
  nrow(emm_plot),
  " overall EMMs; ",
  nrow(sig_tab),
  " significant pairwise brackets.\n",
  sep = ""
)

# Optional export

if (SAVE_FIGURE4) {
  
  if (!requireNamespace("svglite", quietly = TRUE)) {
    stop("Install 'svglite' to export the SVG figure.")
  }
  
  dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  ggsave(
    file.path(OUTPUT_DIR, "Figure4_K_umax_overall.png"),
    p_overall,
    width = 12,
    height = 11,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    file.path(OUTPUT_DIR, "Figure4_K_umax_overall.svg"),
    p_overall,
    device = svglite::svglite,
    fix_text_size = FALSE,
    width = 12,
    height = 11,
    units = "in",
    bg = "white"
  )
}
