# PRISMA 2020 flow diagram (databases & registers arm) for the PCP combination review.
# Reads box counts from prisma_pcp.csv and draws a JAMA blue-grey flow diagram.
# Usage:  source("analysis/prisma_pcp.R")  ->  writes figures/prisma_pcp.png and returns the ggplot.

library(tidyverse)
source("analysis/jama_style.R")

d <- read_csv("analysis/prisma_pcp.csv", show_col_types = FALSE)
getn <- function(key) d$n[d$id == key]

# ---- box display text -------------------------------------------------------
disp <- function(row) {
  lab <- str_wrap(row$label, width = 28)
  cnt <- if (!is.na(row$extra) && nzchar(row$extra))
    paste0("(n = ", row$n, "; ", row$extra, ")") else paste0("(n = ", row$n, ")")
  star <- if (isTRUE(row$provisional)) "*" else ""
  paste0(lab, "\n", cnt, star)
}

# Full-text exclusion box: header + reasons
reasons <- d |> filter(role == "reason")
ft_txt <- paste0(
  "Reports excluded (n = ", getn("excluded_ft"), "):\n",
  paste0("\u2013 ", str_wrap(reasons$label, 30) |>
           str_replace_all("\n", "\n   "),
         " (n = ", reasons$n, ")", collapse = "\n")
)

# ---- geometry ---------------------------------------------------------------
mx <- 4; mhw <- 1.95   # main column x-centre / half-width
sx <- 9; shw <- 2.05   # side column x-centre / half-width

box <- tribble(
  ~id,               ~x, ~y,   ~hh,  ~text,
  "identified",      mx, 10.0, 0.60, disp(d[d$id=="identified",]),
  "screened",        mx,  8.2, 0.55, disp(d[d$id=="screened",]),
  "sought",          mx,  6.4, 0.55, disp(d[d$id=="sought",]),
  "assessed",        mx,  4.6, 0.55, disp(d[d$id=="assessed",]),
  "included",        mx,  2.4, 0.62, disp(d[d$id=="included",]),
  "nma",             mx,  0.7, 0.60, disp(d[d$id=="nma",]),
  "duplicates",      sx,  9.1, 0.55, disp(d[d$id=="duplicates",]),
  "excluded_screen", sx,  7.3, 0.55, disp(d[d$id=="excluded_screen",]),
  "not_retrieved",   sx,  5.5, 0.50, disp(d[d$id=="not_retrieved",]),
  "excluded_ft",     sx,  3.45,1.05, ft_txt
) |>
  mutate(
    xmin = ifelse(x == mx, mx - mhw, sx - shw),
    xmax = ifelse(x == mx, mx + mhw, sx + shw),
    ymin = y - hh, ymax = y + hh,
    side = x == sx
  )

# vertical arrows down the main column
vseg <- tibble(
  x = mx, xend = mx,
  y    = c(9.40, 7.65, 5.85, 3.98, 1.78),
  yend = c(8.75, 6.95, 5.15, 3.02, 1.30)
)
# horizontal exclusion arrows branching to the side column
hseg <- tibble(
  x = mx, xend = sx - shw,
  y    = c(9.1, 7.3, 5.5, 3.5),
  yend = c(9.1, 7.3, 5.5, 3.5)
)

# phase bands (left rotated labels)
band <- tribble(
  ~label,           ~ymin, ~ymax,
  "Identification",  8.45, 10.65,
  "Screening",       3.90,  8.80,
  "Included",        0.05,  3.05
) |> mutate(xmin = 0.15, xmax = 0.75, y = (ymin + ymax) / 2)

arr <- arrow(length = unit(0.16, "cm"), type = "closed")

p <- ggplot() +
  # phase bands
  geom_rect(data = band, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = jama_blue_grey[["steel"]]) +
  geom_text(data = band, aes(x = (xmin + xmax) / 2, y = y, label = label),
            angle = 90, colour = "white", fontface = "bold", size = 3.6) +
  # arrows
  geom_segment(data = vseg, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arr, linewidth = 0.4, colour = jama_blue_grey[["dark"]]) +
  geom_segment(data = hseg, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arr, linewidth = 0.4, colour = jama_blue_grey[["dark"]]) +
  # boxes
  geom_rect(data = box,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = side),
            colour = jama_blue_grey[["dark"]], linewidth = 0.4) +
  geom_text(data = box |> filter(id != "excluded_ft"),
            aes(x = x, y = y, label = text),
            size = 3.0, colour = "black", lineheight = 0.95) +
  geom_text(data = box |> filter(id == "excluded_ft"),
            aes(x = xmin + 0.15, y = y, label = text),
            size = 2.7, colour = "black", lineheight = 0.95, hjust = 0) +
  scale_fill_manual(values = c(`FALSE` = "#DCE4EC", `TRUE` = "#EDF1F5"), guide = "none") +
  coord_cartesian(xlim = c(0, 11.2), ylim = c(-0.1, 10.9), expand = FALSE) +
  labs(caption = "*Provisional \u2014 records identified before duplicate removal and the 52\u201348 (reports not retrieved) step pending confirmation from Covidence.") +
  theme_void(base_family = "sans") +
  theme(plot.caption = element_text(hjust = 0, size = 8, colour = "grey30"),
        plot.margin = margin(6, 6, 6, 6))

dir.create("figures", showWarnings = FALSE)
ggsave("figures/prisma_pcp.png", p, width = 8.5, height = 10, dpi = 200, bg = "white")
p
