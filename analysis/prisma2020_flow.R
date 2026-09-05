# PRISMA 2020 flow diagram via the PRISMA2020 package (Haddaway et al.,
# github.com/prisma-flowdiagram/PRISMA2020). Reads analysis/prisma2020_data.csv
# (the package's CSV template filled with our reconciled Covidence counts) and
# renders the databases-&-registers flow using the package's default (canonical) palette.
# Usage:  source("analysis/prisma2020_flow.R")  ->  writes figures/prisma2020_flowdiagram.png
library(PRISMA2020)

d <- read.csv("analysis/prisma2020_data.csv", stringsAsFactors = FALSE)
d <- PRISMA_data(d)

p <- PRISMA_flowdiagram(
  d,
  previous = FALSE, other = FALSE, meta_analysis = TRUE, fontsize = 10
)

dir.create("figures", showWarnings = FALSE)
PRISMA_save(p, filename = "figures/prisma2020_flowdiagram.png",
            filetype = "PNG", overwrite = TRUE, width = 1600)
p
