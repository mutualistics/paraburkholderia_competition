
# Supplementary Figure S4: fluorescence response in WT competition controls
# Run after the final plate-reader import and K/mu_max processing scripts.
# Points are experimental-date observations, open circles are strain medians,
# and error bars show interquartile ranges. No inferential model is fitted.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

SAVE_FIGURE <- FALSE
OUTPUT_DIR <- file.path("outputs", "figures")
EXCLUDE_CK_WT_DATE <- "2710"

for (object in c("filtered_data", "wt_long")) {
  if (!exists(object, inherits = TRUE)) stop("Missing ", object, ". Run the import and K/mu_max processing first.")
}

required_auc <- c("Date", "RFP_Strain", "GFP_Strain",
                  "Normalized_AUC_RFP", "Normalized_AUC_GFP")
required_wt <- c("Date", "Channel", "Focal_Strain", "metric", "value")

if (!all(required_auc %in% names(filtered_data))) {
  stop("filtered_data is missing: ", paste(setdiff(required_auc, names(filtered_data)), collapse = ", "))
}
if (!all(required_wt %in% names(wt_long))) {
  stop("wt_long is missing: ", paste(setdiff(required_wt, names(wt_long)), collapse = ", "))
}
if (!all(c("dK_pct", "dumax_pct") %in% as.character(wt_long$metric))) {
  stop("WT growth-trait metrics dK_pct and dumax_pct are required.")
}

gfp_order <- c("At", "Bp", "Ck", "Pk", "Sl", "Ss")
rfp_order <- c(
  "Ba159", "Ba161", "Ba31", "Ba80",
  "Bb433", "Bb822", "Bb859", "Bb862",
  "Bh11", "Bh171", "Bh21", "Bh22",
  "Pf", "Ps", "Pu"
)
competitor_names <- c(
  At = "Agrobacterium", Bp = "Brucella", Ck = "Comamonas",
  Pk = "Pseudomonas", Sl = "Serratia", Ss = "Staphylococcus"
)
group_colors <- c(
  GFP = "#696969", "Non-AA" = "#C7C7C7", Ba = "#F0A202",
  Bb = "#6EC6F2", Bh = "#5CC8A1"
)

strain_groups <- tibble(Focal_Strain = rfp_order) %>%
  mutate(Group = case_when(
    grepl("^Ba", Focal_Strain) ~ "Ba",
    grepl("^Bb", Focal_Strain) ~ "Bb",
    grepl("^Bh", Focal_Strain) ~ "Bh",
    TRUE ~ "Non-AA"
  ))

# Keep the established exclusion for the Ck + unlabeled WT control on 2710.
# It applies to this control only, not to any accepted competition date.
auc_controls <- filtered_data %>%
  filter(!(
    as.character(Date) == EXCLUDE_CK_WT_DATE &
      as.character(GFP_Strain) == "Ck" &
      as.character(RFP_Strain) == "wt"
  ))

trait_controls <- wt_long %>%
  filter(!(
    as.character(Date) == EXCLUDE_CK_WT_DATE &
      as.character(Channel) == "GFP" &
      as.character(Focal_Strain) == "Ck"
  ))

# Mean technical replicates within each experimental date for normalized AUC.
gfp_auc <- auc_controls %>%
  filter(as.character(RFP_Strain) == "wt",
         as.character(GFP_Strain) %in% gfp_order,
         is.finite(Normalized_AUC_GFP)) %>%
  group_by(GFP_Strain, Date) %>%
  summarise(value = mean(Normalized_AUC_GFP), .groups = "drop") %>%
  transmute(Channel = "GFP", Focal_Strain = as.character(GFP_Strain),
            metric = "AUC_wt", value)

rfp_auc <- auc_controls %>%
  filter(as.character(GFP_Strain) == "wt",
         as.character(RFP_Strain) %in% rfp_order,
         is.finite(Normalized_AUC_RFP)) %>%
  group_by(RFP_Strain, Date) %>%
  summarise(value = mean(Normalized_AUC_RFP), .groups = "drop") %>%
  transmute(Channel = "RFP", Focal_Strain = as.character(RFP_Strain),
            metric = "AUC_wt", value)

trait_values <- trait_controls %>%
  filter(as.character(metric) %in% c("dK_pct", "dumax_pct"), is.finite(value)) %>%
  transmute(
    Channel = as.character(Channel),
    Focal_Strain = as.character(Focal_Strain),
    metric = recode(as.character(metric),
                    dK_pct = "K_wt_pct", dumax_pct = "umax_wt_pct"),
    value
  )

plot_dat <- bind_rows(gfp_auc, rfp_auc, trait_values) %>%
  filter(Channel %in% c("GFP", "RFP"), is.finite(value)) %>%
  left_join(strain_groups, by = "Focal_Strain") %>%
  mutate(
    Group = if_else(Channel == "GFP", "GFP", coalesce(Group, "Non-AA")),
    Group = factor(Group, levels = names(group_colors))
  )

summary_dat <- plot_dat %>%
  group_by(Channel, metric, Focal_Strain, Group) %>%
  summarise(
    median = median(value), lower = quantile(value, 0.25),
    upper = quantile(value, 0.75), n = n(), .groups = "drop"
  )

cat("\nFigure S4 WT-control observation counts:\n")
print(as.data.frame(count(plot_dat, Channel, metric)), row.names = FALSE)

get_values <- function(metric) {
  x <- plot_dat$value[plot_dat$metric == metric]
  x[is.finite(x)]
}

auc_values <- get_values("AUC_wt")
k_values <- get_values("K_wt_pct")
umax_values <- get_values("umax_wt_pct")
if (!length(auc_values) || !length(k_values) || !length(umax_values)) {
  stop("At least one WT-control metric is empty.")
}

auc_limits <- c(
  min(0, min(auc_values) - 0.03),
  max(1.05, ceiling(max(auc_values) * 20) / 20 + 0.05)
)

percentage_limits <- function(x, reference, bottom, top) {
  c(
    min(bottom, floor(min(x) / 25) * 25 - 5, reference - 5),
    max(top, ceiling(max(x) / 25) * 25 + 7, reference + 5)
  )
}
k_limits <- percentage_limits(k_values, -50, -105, 30)
umax_limits <- percentage_limits(umax_values, 0, -105, 40)
percentage_breaks <- function(limits) seq(
  ceiling(limits[1] / 25) * 25, floor(limits[2] / 25) * 25, by = 25
)

plot_theme <- theme_classic(base_size = 13, base_family = "sans") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "#252525", linewidth = 0.6),
    axis.ticks = element_line(colour = "#252525", linewidth = 0.5),
    axis.ticks.length = grid::unit(2.5, "pt"),
    axis.text.y = element_text(colour = "#252525", size = 11),
    axis.title.y = element_text(colour = "#222222", size = 13.5, margin = margin(r = 14)),
    plot.title = element_text(size = 15, colour = "#111111", margin = margin(b = 12)),
    plot.margin = margin(10, 12, 12, 10),
    legend.position = "none"
  )

make_panel <- function(metric, channel, strain_order, reference, ylim, breaks,
                       ylabel = NULL, title = NULL, show_x = FALSE) {
  d <- plot_dat %>%
    filter(metric == .env$metric, Channel == channel, Focal_Strain %in% strain_order) %>%
    mutate(Focal_Strain = factor(Focal_Strain, levels = strain_order))
  s <- summary_dat %>%
    filter(metric == .env$metric, Channel == channel, Focal_Strain %in% strain_order) %>%
    mutate(Focal_Strain = factor(Focal_Strain, levels = strain_order))
  if (!nrow(d)) stop("No observations for ", metric, " / ", channel)
  
  p <- ggplot(d, aes(Focal_Strain, value))
  if (channel == "RFP") {
    p <- p + geom_vline(xintercept = c(4.5, 8.5, 12.5),
                        colour = "#E0E0E0", linewidth = 0.6)
  }
  
  p <- p +
    geom_hline(yintercept = reference, colour = "#777777",
               linewidth = 0.65, linetype = "dashed") +
    geom_point(aes(colour = Group),
               position = position_jitter(width = 0.10, height = 0, seed = 2026),
               size = 2, shape = 16, alpha = 0.50) +
    geom_errorbar(
      data = s, aes(x = Focal_Strain, ymin = lower, ymax = upper, colour = Group),
      inherit.aes = FALSE, width = 0.19, linewidth = 1.0
    ) +
    geom_point(
      data = s, aes(x = Focal_Strain, y = median, colour = Group),
      inherit.aes = FALSE, shape = 21, fill = "white", stroke = 1.2, size = 3.8
    ) +
    scale_x_discrete(
      limits = strain_order,
      labels = if (channel == "GFP") unname(competitor_names[strain_order]) else strain_order,
      drop = FALSE, expand = expansion(add = c(0.65, 0.65))
    ) +
    scale_y_continuous(
      breaks = breaks,
      labels = if (metric == "AUC_wt") {
        function(x) sprintf("%.2f", x)
      } else {
        function(x) paste0(x, "%")
      }
    ) +
    scale_colour_manual(values = group_colors, guide = "none", drop = FALSE) +
    coord_cartesian(ylim = ylim) +
    labs(x = NULL, y = if (channel == "GFP") ylabel else NULL, title = title) +
    plot_theme +
    theme(
      axis.text.x = if (show_x) element_text(
        angle = if (channel == "GFP") 43 else 50,
        hjust = 1, vjust = 1,
        size = if (channel == "GFP") 11 else 10.5,
        colour = "#222222", margin = margin(t = 7)
      ) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      plot.title = if (is.null(title)) element_blank() else element_text()
    )
  p
}

p_auc_gfp <- make_panel(
  "AUC_wt", "GFP", gfp_order, 0.5, auc_limits,
  seq(0, floor(auc_limits[2] / 0.25) * 0.25, by = 0.25),
  ylabel = "Normalized AUC",
  title = expression(bold("(A)") ~ "Competitor strains (GFP)")
)
p_auc_rfp <- make_panel(
  "AUC_wt", "RFP", rfp_order, 0.5, auc_limits,
  seq(0, floor(auc_limits[2] / 0.25) * 0.25, by = 0.25),
  title = expression(bold("(B)") ~ italic(Paraburkholderia) ~ "strains (RFP)")
)
p_k_gfp <- make_panel(
  "K_wt_pct", "GFP", gfp_order, -50, k_limits,
  percentage_breaks(k_limits), ylabel = expression(Delta * K ~ "(%)")
)
p_k_rfp <- make_panel(
  "K_wt_pct", "RFP", rfp_order, -50, k_limits,
  percentage_breaks(k_limits)
)
p_umax_gfp <- make_panel(
  "umax_wt_pct", "GFP", gfp_order, 0, umax_limits,
  percentage_breaks(umax_limits), ylabel = expression(Delta * mu[max] ~ "(%)"),
  show_x = TRUE
)
p_umax_rfp <- make_panel(
  "umax_wt_pct", "RFP", rfp_order, 0, umax_limits,
  percentage_breaks(umax_limits), show_x = TRUE
)

row_auc <- (p_auc_gfp | p_auc_rfp) + plot_layout(widths = c(1, 1.65))
row_k <- (p_k_gfp | p_k_rfp) + plot_layout(widths = c(1, 1.65))
row_umax <- (p_umax_gfp | p_umax_rfp) + plot_layout(widths = c(1, 1.65))
final_combined_plot <- (row_auc / plot_spacer() / row_k / plot_spacer() / row_umax) +
  plot_layout(heights = c(1.05, 0.09, 1.00, 0.09, 1.15))

print(final_combined_plot)

if (SAVE_FIGURE) {
  if (!requireNamespace("svglite", quietly = TRUE)) stop("Install svglite for SVG export.")
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(OUTPUT_DIR, "Fig_S4.svg"), final_combined_plot,
         device = svglite::svglite, width = 12, height = 8, units = "in", bg = "white")
}
