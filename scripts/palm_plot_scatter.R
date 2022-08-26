#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(argparser)) # v0.6
quiet(library(tidyverse)) # v1.3.1

# load the helper functions
source("scripts/clues_utils.R")

# get the command line arguments
p <- arg_parser("Plot a scatter plot of -log10(pval) vs delra_prs")
p <- add_argument(p, "--palm", help = "PALM report for all ancestries", default = "results/palm/ancestral_paths_new-ms-all-palm_report_prs.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "results/palm/ancestral_paths_new-ms-all-palm_report_prs-scatter.png")

argv <- parse_args(p)

palm <- read_tsv(argv$palm, col_types = cols()) %>%
    mutate(logp = ifelse(delta_prs > 0, -log10(p.value), log10(p.value)))

ggplot(palm, aes(x=logp, y=delta_prs, colour=ancestry)) +
    geom_point(shape=1) +
    geom_smooth(aes(fill=ancestry), method='lm', formula="y ~ x", fullrange=TRUE) +

    # display as a single row
    facet_wrap(~ancestry, nrow = 1) +

    # set the colour scales
    scale_fill_manual(values = ancestry_colors) +
    scale_color_manual(values = ancestry_colors) +
        
    ylab("Δ PRS per SNP") +
    xlab("Directional -log10(pval)") +
    
    # basic styling
    theme_minimal() +
    theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white", color = "white"),
        panel.spacing = unit(1, "lines")
    )
