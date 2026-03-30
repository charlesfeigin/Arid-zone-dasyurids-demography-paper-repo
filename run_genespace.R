#!/usr/bin/env Rscript
library(GENESPACE)

genespace_project <- init_genespace(
  wd = getwd(),
  path2mcscanx = "/data/tigrr/home/cfeigin/src/MCScanX-1.0.0/")

processed_synteny <- run_genespace(genespace_project)

svg("riparian/riparian_plot.svg", width = 10, height = 7)
plot_riparian(
  gsParam = processed_synteny,
  refGenome = "Antechinus",
  backgroundColor = "white",
  chrFill = "grey",
  chrBorderCol = "black",
  chrBorderLwd = 0.5
)
dev.off()
