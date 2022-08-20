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
quiet(library(jsonlite)) # v1.8.0
quiet(library(directlabels)) # v2021.1.13
quiet(library(zoo)) # v1.8_10

# load the helper functions
source("scripts/clues_utils.R")

# DRB1*15:01 tag SNPs
# report_tsv <- "results/palm/ancestral_paths_new-ra-all-palm_report_prs.tsv"
# output_png <- "DRB1-15-01.png"
# rsid_list <- c("rs3129934", "rs3135391", "rs3135388", "rs3129889", "rs9271366")

# TB SNPs
report_tsv <- "results/palm/ancestral_paths_new-ms-all-palm_report_prs.tsv"
output_png <- "TB-SNPs.png"
rsid_list <- c("rs422951", "rs3135363", "rs3806156", "rs2395162", "rs9268530", "rs2596472", "rs2844494")

# load the report, with all the results
snps <- read_tsv(report_tsv, col_types = cols()) %>%
    # compose the model prefixes from the PALM metadata
    mutate(prefix = paste0("results/clues/", rsid, "/ancestral_paths_new-", chrom, ":", pos, ":", ancestral_allele, ":", derived_allele, "-", ancestry))

# get the models to plot
snps <- snps %>%
    filter(rsid %in% rsid_list) %>%
    select(rsid, ancestry, prefix)

models <- list()
for (i in 1:nrow(snps)) {
    model <- clues_trajectory(snps[i, ]$rsid, snps[i, ]$ancestry, snps[i, ]$prefix)
    models <- append(models, list(model))
}

# load all the trajectories
traj <- bind_rows(models)

gen_time <- 28
max_age <- 13665
smooth <- 10

# constrain the extent of the plotting
xmin <- round(-max_age / gen_time)
xmax <- max(traj$epoch)
xbreaks <- -seq(-xmax, -xmin, round(2000 / gen_time))
xlabels <- round(xbreaks * gen_time / 1000)

traj %>%
    # plot the heatmap
    ggplot(aes(x = epoch, y = freq)) +

    # plot the maximum posterior trajectory
    geom_line(aes(color = rsid), cex = 1, na.rm = TRUE) +
    facet_wrap(~ancestry, nrow = 1) +

    # set the SNP colors
    scale_color_manual(values = snp_colors) +

    # print the labels
    geom_dl(aes(label = rsid, color = rsid), method = list(dl.trans(x = x + 0.1), "last.qp", cex = 0.8), na.rm = TRUE) +

    # set the axis breaks
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), position = "left", expand = expansion(add = c(0.03, 0))) +
    scale_x_continuous(limits = c(xmin, xmax), breaks = xbreaks, labels = xlabels, expand = expansion(add = c(0, 230))) +

    # labs(title = title) +
    ylab("DAF") +
    xlab("kyr BP") +

    # basic styling
    theme_minimal() +
    theme(
        legend.position = "none",
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white", color = "white")
    )

# save the plot
ggsave(output_png, width = 16, height = 4)
