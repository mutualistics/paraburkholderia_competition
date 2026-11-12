
# Supplementary Figures S5 and S6: K and maximum early growth rate by competitor
# Run after 02_K_umax_models_and_EMMs_282.R in the same R session.
# Within-competitor inference comes from the same models as Figure 4.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
})

SAVE_SUPPLEMENT <- FALSE
OUTPUT_DIR <- file.path("outputs", "figures")

required <- c("kr_within_emm", "kr_within_pairs", "kr_raw_within")
missing <- required[!vapply(required, exists, logical(1), inherits = TRUE)]

if (length(missing)) {
  stop(
    "Run 02_K_umax_models_and_EMMs_282.R first: ",
    paste(missing, collapse = ", ")
  )
}

plot_group_levels <- c("Non-AA", "Ba", "Bb", "Bh")
comp_levels <- c("At", "Bp", "Ck", "Pk", "Sl", "Ss")
response_levels <- c("Paraburkholderia response", "Competitor response")
trait_levels <- c("K", "umax")

BRACKET_BASE_PAD <- 0.26
BRACKET_STEP <- 0.26

group_colors <- c(
  "Non-AA" = "#C7C7C7",
  "Ba" = "#F0A202",
  "Bb" = "#6EC6F2",
  "Bh" = "#5CC8A1"
)

competitor_names <- c(
  At = "Agrobacterium",
  Bp = "Brucella",
  Ck = "Comamonas",
  Pk = "Pseudomonas",
  Sl = "Serratia",
  Ss = "Staphylococcus"
)

group_axis_labels <- expression(
  "Non-AA",
  italic("P. agricolaris"),
  italic("P. bonniea"),
  italic("P. hayleyella")
)

# Plot factors

add_plot_factors <- function(x) {
  x %>%
    mutate(
      Competitor = factor(
        as.character(Competitor),
        levels = comp_levels
      ),
      Competitor_lab = factor(
        competitor_names[as.character(Competitor)],
        levels = unname(competitor_names[comp_levels])
      ),
      Response = factor(
        as.character(Response),
        levels = response_levels
      ),
      Trait = factor(
        as.character(Trait),
        levels = trait_levels
      )
    )
}

# Appearance

supp_theme <- function() {
  theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(
        colour = "black", fill = NA, linewidth = 0.65
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_line(
        colour = "black", linewidth = 0.55
      ),
      axis.text.x = element_text(
        angle = 60,
        hjust = 1,
        vjust = 1,
        size = 10,
        colour = "grey20",
        margin = margin(t = 7)
      ),
      axis.text.y = element_text(
        size = 11, colour = "grey20"
      ),
      axis.title.y = element_text(
        size = 15,
        colour = "black",
        margin = margin(r = 9)
      ),
      strip.background = element_blank(),
      strip.text.x = element_text(
        size = 14,
        face = "bold",
        colour = "black",
        margin = margin(b = 7)
      ),
      strip.text.y.right = element_text(
        size = 12,
        face = "bold",
        colour = "black",
        angle = 270,
        margin = margin(l = 5)
      ),
      panel.spacing.x = grid::unit(0.65, "lines"),
      panel.spacing.y = grid::unit(0.8, "lines"),
      legend.position = "none",
      plot.margin = margin(
        t = 14, r = 18, b = 13, l = 12
      )
    )
}

# Two response rows by six competitor columns.

make_response_plot <- function(trait, y_label) {
  
  raw <- kr_raw_within %>%
    filter(Trait == trait) %>%
    add_plot_factors() %>%
    mutate(
      Plot_Group = factor(
        as.character(Plot_Group),
        levels = plot_group_levels
      )
    )
  
  emm <- kr_within_emm %>%
    filter(Trait == trait) %>%
    add_plot_factors() %>%
    mutate(
      Plot_Group = factor(
        as.character(Plot_Group),
        levels = plot_group_levels
      )
    )
  
  pairs_all <- kr_within_pairs %>%
    filter(Trait == trait) %>%
    add_plot_factors() %>%
    mutate(
      group1 = factor(group1, levels = plot_group_levels),
      group2 = factor(group2, levels = plot_group_levels)
    )
  
  # Keep all tests; brackets show only Holm p < 0.05.
  
  pairs_sig <- pairs_all %>%
    filter(!is.na(label))
  
  # Four EMMs and six pairwise tests per competitor and response.
  
  if (nrow(emm) != 48L ||
      nrow(pairs_all) != 72L ||
      anyNA(emm[c(
        "Response", "Competitor_lab", "Plot_Group"
      )]) ||
      anyNA(pairs_all[c(
        "Response", "Competitor_lab", "group1", "group2"
      )]) ||
      any(!is.finite(emm$emmean_log2)) ||
      any(!is.finite(pairs_all$p_holm)) ||
      anyDuplicated(emm[c(
        "Response", "Competitor_lab", "Plot_Group"
      )])) {
    stop(
      "Expected four EMMs and six contrasts per competitor ",
      "and response for ", trait
    )
  }
  
  if (nrow(pairs_sig) &&
      (anyNA(pairs_sig$group1) ||
       anyNA(pairs_sig$group2) ||
       anyNA(pairs_sig$Competitor_lab) ||
       anyNA(pairs_sig$Response))) {
    stop(
      "Significant contrasts have missing bracket coordinates for ",
      trait
    )
  }
  
  # Brackets positioned independently within each panel.
  
  panel_tops <- bind_rows(
    raw %>%
      transmute(
        Response,
        Competitor_lab,
        top = y_log2
      ),
    emm %>%
      transmute(
        Response,
        Competitor_lab,
        top = upper_log2
      )
  ) %>%
    group_by(Response, Competitor_lab) %>%
    summarise(
      top = max(top, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (nrow(pairs_sig)) {
    
    pairs_sig <- pairs_sig %>%
      left_join(
        panel_tops,
        by = c("Response", "Competitor_lab")
      ) %>%
      group_by(Response, Competitor_lab) %>%
      arrange(p_holm, .by_group = TRUE) %>%
      mutate(
        y.position = top + BRACKET_BASE_PAD +
          BRACKET_STEP * (row_number() - 1L)
      ) %>%
      ungroup()
    
  } else {
    pairs_sig$y.position <- numeric(0)
  }
  
  if (nrow(pairs_sig) &&
      any(!is.finite(pairs_sig$y.position))) {
    stop(
      "Non-finite significance-bracket heights for ",
      trait
    )
  }
  
  cat(
    "\n", trait, ": ",
    nrow(pairs_sig),
    " Holm-significant pairwise brackets.\n",
    sep = ""
  )
  
  # Shared limits across both responses and all competitors.
  
  all_y <- c(
    raw$y_log2,
    emm$lower_log2,
    emm$upper_log2,
    pairs_sig$y.position,
    0
  )
  
  all_y <- all_y[is.finite(all_y)]
  
  if (!length(all_y)) {
    stop("No finite plotting values for ", trait)
  }
  
  y_limits <- c(
    floor((min(all_y) - 0.25) * 2) / 2,
    ceiling((max(all_y) + 0.55) * 2) / 2
  )
  
  y_breaks <- pretty(y_limits, n = 6)
  
  y_breaks <- y_breaks[
    y_breaks >= y_limits[1] &
      y_breaks <= y_limits[2]
  ]
  
  p <- ggplot() +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.75,
      colour = "black"
    ) +
    
    # Date-level observations.
    
    geom_point(
      data = raw,
      aes(
        x = Plot_Group,
        y = y_log2
      ),
      position = position_jitter(
        width = 0.11,
        height = 0,
        seed = 17
      ),
      shape = 21,
      size = 2.1,
      stroke = 0.4,
      alpha = 0.42,
      colour = "grey45",
      fill = "grey75"
    ) +
    
    # Model-adjusted EMMs and 95% confidence intervals.
    
    geom_errorbar(
      data = emm,
      aes(
        x = Plot_Group,
        ymin = lower_log2,
        ymax = upper_log2
      ),
      width = 0.17,
      linewidth = 1.05,
      colour = "black"
    ) +
    
    geom_point(
      data = emm,
      aes(
        x = Plot_Group,
        y = emmean_log2,
        fill = Plot_Group
      ),
      shape = 23,
      size = 5.3,
      stroke = 1.05,
      colour = "black"
    ) +
    
    # Black centres indicate Holm-significant differences from WT.
    
    geom_point(
      data = emm %>% filter(sig_vs_zero),
      aes(
        x = Plot_Group,
        y = emmean_log2
      ),
      shape = 23,
      size = 2.05,
      stroke = 0.45,
      fill = "black",
      colour = "black"
    ) +
    
    facet_grid(
      Response ~ Competitor_lab,
      drop = FALSE
    ) +
    
    scale_fill_manual(
      values = group_colors,
      guide = "none",
      drop = FALSE
    ) +
    
    scale_x_discrete(
      limits = plot_group_levels,
      labels = group_axis_labels,
      drop = FALSE
    ) +
    
    scale_y_continuous(
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0))
    ) +
    
    coord_cartesian(
      ylim = y_limits,
      clip = "off"
    ) +
    
    labs(
      x = NULL,
      y = y_label
    ) +
    
    supp_theme()
  
  # Pairwise significance brackets.
  
  if (nrow(pairs_sig)) {
    
    p <- p +
      ggpubr::stat_pvalue_manual(
        data = pairs_sig,
        label = "label",
        xmin = "group1",
        xmax = "group2",
        y.position = "y.position",
        tip.length = 0.008,
        bracket.size = 0.45,
        size = 4,
        inherit.aes = FALSE
      )
  }
  
  list(
    plot = p,
    raw = raw,
    emm = emm,
    contrast_all = pairs_all,
    sig = pairs_sig,
    ylim = y_limits
  )
}

# Generate Figures S5 and S6.

K_plot_out <- make_response_plot(
  trait = "K",
  y_label = expression(
    log[2] ~ "fold-change in K vs WT"
  )
)

umax_plot_out <- make_response_plot(
  trait = "umax",
  y_label = expression(
    log[2] ~ "fold-change in" ~ mu[max] ~ "vs WT"
  )
)

p_s5 <- K_plot_out$plot
p_s6 <- umax_plot_out$plot

p_s5_s6 <- p_s5 / p_s6 +
  plot_layout(heights = c(1, 1))

print(p_s5)
print(p_s6)

cat(
  "\nSupplement: ",
  nrow(kr_within_emm),
  " EMMs; ",
  sum(!is.na(kr_within_pairs$label)),
  " significant pairwise brackets.\n",
  sep = ""
)

# Optional exports

if (SAVE_SUPPLEMENT) {
  
  if (!requireNamespace("svglite", quietly = TRUE)) {
    stop("Install 'svglite' to export SVG figures.")
  }
  
  dir.create(
    OUTPUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  ggsave(
    file.path(OUTPUT_DIR, "Fig_S5.svg"),
    p_s5,
    device = svglite::svglite,
    fix_text_size = FALSE,
    width = 15,
    height = 7.5,
    units = "in",
    bg = "white"
  )
  
  ggsave(
    file.path(OUTPUT_DIR, "Fig_S6.svg"),
    p_s6,
    device = svglite::svglite,
    fix_text_size = FALSE,
    width = 15,
    height = 7.5,
    units = "in",
    bg = "white"
  )
}
