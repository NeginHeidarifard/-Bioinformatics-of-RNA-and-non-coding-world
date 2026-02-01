#!/usr/bin/env Rscript

# KDE analysis for C3'-C3' distances (Task 3)
# Author: Negin Heidarifard

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript kde_plot.R input.csv output.png")
}

input_csv <- args[1]
output_png <- args[2]

data <- read.csv(input_csv)

png(output_png, width = 800, height = 600)

hist(
  data$distance,
  breaks = seq(0, 20, by = 1),
  freq = FALSE,
  col = "lightgray",
  border = "white",
  xlab = "C3'-C3' distance (Å)",
  main = "Histogram and KDE of C3'-C3' distances (AU)"
)

dens <- density(data$distance, bw = 1.0)
lines(dens, col = "blue", lwd = 2)

dev.off()
