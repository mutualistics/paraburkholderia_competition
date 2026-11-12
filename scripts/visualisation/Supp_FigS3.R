
# Supplementary Figure S3: monoculture maximum OD and growth rate
# Run after the final plate-reader import (full_data_corrected).
# Technical wells are averaged within experimental dates by default.
# Maximum OD is observed; growth rate r is fitted with growthcurver.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(growthcurver)
  library(patchwork)
})

# Settings
SAVE_FIGURE <- FALSE
OUTPUT_DIR <- file.path("outputs", "figures")
BOX_UNIT <- "date"             # "date" or "well"
MAX_OD_SOURCE <- "observed"    # "observed" or "fitted_K"
MIN_BOX_N <- 3
OD_FLOOR <- 1e-4
RATE_YMAX <- 0.70
FONT_FAMILY <- "Arial"

if (!BOX_UNIT %in% c("date", "well") ||
    !MAX_OD_SOURCE %in% c("observed", "fitted_K")) {
  stop("Check BOX_UNIT and MAX_OD_SOURCE.")
}

if (!exists("full_data_corrected", inherits = TRUE)) {
  stop("Run the plate-reader import first (full_data_corrected missing).")
}

required <- c("Date", "Well_ID", "RFP_Strain", "GFP_Strain", "Time", "OD")
missing <- setdiff(required, names(full_data_corrected))
if (length(missing)) stop("Missing corrected-data columns: ", paste(missing, collapse = ", "))

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
  GFP = "#696969", "Non-AA" = "#7D7D7D", Ba = "#F0A202",
  Bb = "#6EC6F2", Bh = "#5CC8A1"
)
group_fills <- c(
  GFP = "#BEBEBE", "Non-AA" = "#C7C7C7", Ba = "#F0A202",
  Bb = "#6EC6F2", Bh = "#5CC8A1"
)

# The figure contains monocultures only. Competition-run exclusion lists do
# not apply to these curves; the globally excluded 0805 plate is omitted.
od_mono <- full_data_corrected %>%
  mutate(
    Date = trimws(as.character(Date)),
    RFP_Strain = as.character(RFP_Strain),
    GFP_Strain = as.character(GFP_Strain),
    RFP_Strain = if_else(RFP_Strain == "YG", "None", RFP_Strain),
    GFP_Strain = if_else(GFP_Strain == "YG", "None", GFP_Strain),
    Time = as.numeric(as.character(Time)),
    OD_use = pmax(as.numeric(OD), OD_FLOOR),
    Channel = case_when(
      RFP_Strain == "None" & GFP_Strain %in% gfp_order ~ "GFP",
      GFP_Strain == "None" & RFP_Strain %in% rfp_order ~ "RFP",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(Date != "0805", !is.na(Channel), is.finite(Time), is.finite(OD_use)) %>%
  mutate(
    Focal_Strain = if_else(Channel == "GFP", GFP_Strain, RFP_Strain),
    Group = case_when(
      Channel == "GFP" ~ "GFP",
      grepl("^Ba", Focal_Strain) ~ "Ba",
      grepl("^Bb", Focal_Strain) ~ "Bb",
      grepl("^Bh", Focal_Strain) ~ "Bh",
      TRUE ~ "Non-AA"
    )
  )

if (!nrow(od_mono)) stop("No monoculture OD curves found.")

# Fit the logistic growth rate for the monocultures plotted in this figure.
fit_od <- function(t, y) {
  t <- as.numeric(t)
  y <- pmax(as.numeric(y), OD_FLOOR)
  ok <- is.finite(t) & is.finite(y)
  t <- t[ok]
  y <- y[ok]
  ord <- order(t)
  t <- t[ord]
  y <- y[ord]
  
  tryCatch({
    v <- growthcurver::SummarizeGrowth(data_t = t, data_n = y)$vals
    list(ok = TRUE, r = as.numeric(v["r"]),
         K = as.numeric(v["k"]), ObservedMax = max(y))
  }, error = function(e) {
    list(ok = FALSE, r = NA_real_, K = NA_real_, ObservedMax = NA_real_)
  })
}

od_fits <- od_mono %>%
  group_by(Date, Well_ID, Channel, Focal_Strain, Group) %>%
  summarise(fit = list(fit_od(Time, OD_use)), .groups = "drop") %>%
  mutate(
    fit_ok = map_lgl(fit, "ok"),
    r = map_dbl(fit, "r"),
    K = map_dbl(fit, "K"),
    ObservedMax = map_dbl(fit, "ObservedMax")
  ) %>%
  select(-fit) %>%
  filter(fit_ok) %>%
  mutate(MaxOD = if (MAX_OD_SOURCE == "observed") ObservedMax else K)

well_data <- od_fits %>%
  select(Date, Well_ID, Channel, Focal_Strain, Group, MaxOD, r) %>%
  pivot_longer(c(MaxOD, r), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric, r = "GrowthRate")) %>%
  filter(is.finite(Value))

date_data <- well_data %>%
  group_by(Date, Channel, Focal_Strain, Group, Metric) %>%
  summarise(Value = mean(Value), .groups = "drop")

plot_data <- if (BOX_UNIT == "date") date_data else well_data
plot_data <- plot_data %>%
  add_count(Channel, Focal_Strain, Metric, name = "n_units") %>%
  mutate(
    x = if_else(
      Channel == "GFP",
      match(Focal_Strain, gfp_order),
      match(Focal_Strain, rfp_order)
    ),
    x_raw = x + if_else(Group == "Non-AA", 0.18, 0)
  )

if (!nrow(plot_data) || anyNA(plot_data$x)) stop("Incomplete OD plotting data.")

cat("\nFigure S3: ", n_distinct(od_mono$Date), " experimental dates; ",
    n_distinct(od_fits$Well_ID, od_fits$Date), " fitted monoculture wells.\n", sep = "")
print(as.data.frame(count(plot_data, Channel, Metric, Focal_Strain)), row.names = FALSE)

above_axis <- plot_data %>% filter(Metric == "GrowthRate", Value > RATE_YMAX)
cat("Growth-rate observations above displayed axis:", nrow(above_axis), "\n")
if (nrow(above_axis)) {
  print(as.data.frame(select(above_axis, Channel, Focal_Strain, Date, Value)),
        row.names = FALSE)
}

max_od <- plot_data$Value[plot_data$Metric == "MaxOD"]
if (!length(max_od)) stop("No maximum OD measurements available.")
max_top <- max(1.8, ceiling((max(max_od) + 0.04) / 0.2) * 0.2)

plot_theme <- theme_classic(base_size = 12, base_family = FONT_FAMILY) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(colour = "#252525", linewidth = 0.60),
    axis.ticks = element_line(colour = "#252525", linewidth = 0.48),
    axis.ticks.length = grid::unit(2.4, "pt"),
    axis.text.y = element_text(size = 10.5, colour = "#252525", margin = margin(r = 6)),
    axis.title.y = element_text(size = 13.5, colour = "#222222", margin = margin(r = 13)),
    plot.title = element_text(size = 14.5, colour = "#111111", margin = margin(b = 10)),
    plot.margin = margin(9, 12, 9, 9),
    legend.position = "none"
  )

make_panel <- function(metric, channel, strain_order, ylim, breaks,
                       ylabel = NULL, title = NULL, show_x = FALSE) {
  d <- filter(plot_data, Metric == metric, Channel == channel,
              Focal_Strain %in% strain_order)
  if (!nrow(d)) stop("No data for ", metric, " / ", channel)
  
  boxes <- filter(d, n_units >= MIN_BOX_N)
  visible <- filter(d, metric != "GrowthRate" | Value <= RATE_YMAX)
  high <- filter(d, metric == "GrowthRate", Value > RATE_YMAX)
  labels <- if (channel == "GFP") {
    parse(text = paste0("italic('", unname(competitor_names[strain_order]), "')"))
  } else strain_order
  
  p <- ggplot()
  if (channel == "RFP") {
    p <- p + geom_vline(xintercept = c(4.5, 8.5, 12.5),
                        colour = "#DEDEDE", linewidth = 0.55)
  }
  if (nrow(boxes)) {
    # All values, including off-axis growth rates, enter the boxplot statistics.
    p <- p + geom_boxplot(
      data = boxes,
      aes(x = x, y = Value, group = Focal_Strain, colour = Group, fill = Group),
      width = 0.48, coef = 1.5, outlier.shape = NA,
      alpha = 0.16, linewidth = 0.72, na.rm = TRUE
    )
  }
  p <- p + geom_point(
    data = visible, aes(x = x_raw, y = Value, colour = Group),
    position = position_jitter(width = 0.12, height = 0, seed = 2026),
    shape = 16, size = 1.75, alpha = 0.52
  )
  if (nrow(high)) {
    # Triangles denote values beyond the plotted growth-rate range.
    high <- mutate(high, display_y = RATE_YMAX - 0.025)
    p <- p + geom_point(
      data = high, aes(x = x_raw, y = display_y, colour = Group),
      position = position_jitter(width = 0.12, height = 0, seed = 2026),
      shape = 2, size = 2.6, stroke = 0.8
    )
  }
  
  p +
    scale_x_continuous(
      limits = c(0.05, length(strain_order) + 0.95),
      breaks = seq_along(strain_order), labels = labels,
      minor_breaks = NULL, expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      breaks = breaks, labels = function(x) sprintf("%.1f", x),
      minor_breaks = NULL, expand = expansion(mult = 0)
    ) +
    scale_colour_manual(values = group_colors, guide = "none") +
    scale_fill_manual(values = group_fills, guide = "none") +
    coord_cartesian(ylim = ylim, clip = "on") +
    labs(x = NULL, y = if (channel == "GFP") ylabel else NULL, title = title) +
    plot_theme +
    theme(
      axis.text.x = if (show_x) element_text(
        angle = if (channel == "GFP") 43 else 50,
        hjust = 1, vjust = 1,
        size = if (channel == "GFP") 10.5 else 10,
        colour = "#252525", margin = margin(t = 7)
      ) else element_blank(),
      axis.ticks.x = if (show_x) element_line(colour = "#252525") else element_blank(),
      plot.title = if (is.null(title)) element_blank() else element_text()
    )
}

max_label <- if (MAX_OD_SOURCE == "observed") "Maximum OD" else
  expression("Fitted carrying capacity " * italic(K))

p_max_gfp <- make_panel(
  "MaxOD", "GFP", gfp_order, c(0, max_top), seq(0, max_top, by = 0.4),
  ylabel = max_label,
  title = expression(bold("(A)") ~ "Competitor strains (GFP)")
)
p_max_rfp <- make_panel(
  "MaxOD", "RFP", rfp_order, c(0, max_top), seq(0, max_top, by = 0.4),
  title = expression(bold("(B)") ~ italic(Paraburkholderia) ~ "strains (RFP)")
)
p_rate_gfp <- make_panel(
  "GrowthRate", "GFP", gfp_order, c(0, RATE_YMAX),
  seq(0, RATE_YMAX + 1e-8, by = 0.2),
  ylabel = expression("Growth rate " * italic(r) ~ "(h"^{-1} * ")"),
  show_x = TRUE
)
p_rate_rfp <- make_panel(
  "GrowthRate", "RFP", rfp_order, c(0, RATE_YMAX),
  seq(0, RATE_YMAX + 1e-8, by = 0.2), show_x = TRUE
)

row_max <- (p_max_gfp | p_max_rfp) + plot_layout(widths = c(1, 1.65))
row_rate <- (p_rate_gfp | p_rate_rfp) + plot_layout(widths = c(1, 1.65))
final_box_plot <- (row_max / plot_spacer() / row_rate) +
  plot_layout(heights = c(1, 0.055, 1.15))

print(final_box_plot)

if (SAVE_FIGURE) {
  if (!requireNamespace("svglite", quietly = TRUE)) stop("Install svglite for SVG export.")
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(OUTPUT_DIR, "Fig_S3.svg"), final_box_plot,
         device = svglite::svglite, width = 12, height = 8, units = "in", bg = "white")
}
