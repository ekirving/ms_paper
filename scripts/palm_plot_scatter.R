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
quiet(library(ggrepel)) # v0.9.1

# load the helper functions
source("scripts/clues_utils.R")

# get the command line arguments
p <- arg_parser("Plot a scatter plot of -log10(pval) vs delra_prs")
p <- add_argument(p, "--palm", help = "PALM report for all ancestries", default = "results/palm/ancestral_paths_new-ms-r0.05-kb250-palm_report_prs.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "results/palm/ancestral_paths_new-ms-r0.05-kb250-scatter.png")

argv <- parse_args(p)

# define a delta PRS threshold for labeling SNPs
delta_prs_threshold <- 0.005

palm <- read_tsv(argv$palm, col_types = cols()) %>%
    # calculate the directional p-value
    mutate(logp = ifelse(delta_prs > 0, -log10(p.value), log10(p.value))) %>%
    # label SNPs above the threshold that are not significant
    mutate(label = ifelse(abs(delta_prs) >= delta_prs_threshold & !significant, rsid, ""))

ymin <- floor(min(palm$delta_prs) / delta_prs_threshold) * delta_prs_threshold
ymax <- ceiling(max(palm$delta_prs) / delta_prs_threshold) * delta_prs_threshold

plt <- ggplot(palm, aes(x = logp, y = delta_prs, colour = ancestry)) +
    geom_point(aes(shape = significant)) +
    geom_smooth(aes(fill = ancestry), method = "lm", formula = "y ~ x", fullrange = TRUE) +
    # Add auto-positioned text
    geom_text_repel(
        aes(label = label),
        size = 9 / .pt, # font size 9 pt
        point.padding = 0.3,
        box.padding = 0.7,
        min.segment.length = 0,
        max.overlaps = 1000,
        bg.color = "white"
    ) +

    # display as a single row
    facet_wrap(~ancestry, nrow = 1) +

    # set the colour scales
    scale_fill_manual(values = ancestry_colors) +
    scale_color_manual(values = ancestry_colors) +
    scale_shape_manual(values = c(1, 20)) +

    # set the axis breaks
    scale_y_continuous(limits = c(ymin, ymax), breaks = seq(ymin, ymax, delta_prs_threshold)) +
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

# save the plot
ggsave(argv$output, plot = plt, width = 16, height = 4) %>% suppressWarnings()
