# JAMA-style figure helpers for the PCP combination meta-analysis.
# Usage:  source("analysis/jama_style.R"); ... + theme_jama()
#         scale_color_manual(values = jama_blue_grey)
library(ggplot2)

# Blue-grey palette (light -> dark). Use jama_blue_grey for discrete groups,
# or index specific shades (e.g. jama_blue_grey[["dark"]]) for single-colour plots.
jama_blue_grey <- c(
  light  = "#A9BCCF",
  medium = "#6E8CA8",
  steel  = "#3C5F80",
  dark   = "#243B53",
  accent = "#16324F"
)

# Tier mapping used in the forest plots (primary / sensitivity / pooled).
jama_tier <- c(
  primary     = "#243B53",  # dark blue-grey
  sensitivity = "#6E8CA8",  # medium blue-grey
  pooled      = "#16324F"   # near-navy, for the summary estimate
)

theme_jama <- function(base_size = 12) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line      = element_line(colour = "black", linewidth = 0.4),
      axis.ticks     = element_line(colour = "black", linewidth = 0.4),
      axis.text      = element_text(colour = "black"),
      axis.title     = element_text(colour = "black"),
      plot.title     = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle  = element_text(colour = "grey30"),
      panel.grid     = element_blank(),
      legend.key     = element_blank(),
      legend.title   = element_blank(),
      legend.position = "right"
    )
}
