#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2023, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(ggpubr)) # 0.4.0

# ----------------------------------------------------------------------------

# Figure 5
source("scripts/palm_plot_lines_pval.R")
source("scripts/clues_plot_snps.R")

ggarrange(plt_lines, plt_snps, ncol=1, heights = c(3, 1))

# save the plots
ggsave(filename = "figure_hires/figure_5_18x17.pdf", width = 18, height = 17, units = "cm", dpi = 300, device=cairo_pdf, scale=1.2)
ggsave(filename = "figure_hires/figure_5_18x16.pdf", width = 18, height = 16, units = "cm", dpi = 300, device=cairo_pdf, scale=1.2)
ggsave(filename = "figure_hires/figure_5_18x15.pdf", width = 18, height = 15, units = "cm", dpi = 300, device=cairo_pdf, scale=1.2)
ggsave(filename = "figure_hires/figure_5_18x14.pdf", width = 18, height = 14, units = "cm", dpi = 300, device=cairo_pdf, scale=1.2)
ggsave(filename = "figure_hires/figure_5_18x13.pdf", width = 18, height = 13, units = "cm", dpi = 300, device=cairo_pdf, scale=1.2)
ggsave(filename = "figure_hires/figure_5_18x12.pdf", width = 18, height = 12, units = "cm", dpi = 300, device=cairo_pdf, scale=1.2)
ggsave(filename = "figure_hires/figure_5_18x11.pdf", width = 18, height = 11, units = "cm", dpi = 300, device=cairo_pdf, scale=1.2)
ggsave(filename = "figure_hires/figure_5_18x10.pdf", width = 18, height = 10, units = "cm", dpi = 300, device=cairo_pdf, scale=1.2)

# ----------------------------------------------------------------------------

# Supp_Figure_6.1
source("scripts/palm_plot_lines_pval.R")

plt_lines

