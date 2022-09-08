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

# # DRB1*15:01 tag SNPs
# description <- "DRB1*15:01"
# report_tsv <- "results/palm/ancestral_paths_new-ra-all-palm_report_prs.tsv"
# output_png <- "DRB1-15-01.png"
# rsid_list <- c("rs3129934", "rs3135391", "rs3135388", "rs3129889", "rs9271366")

# # TB SNPs
# description <- "Tuberculosis (TB)"
# report_tsv <- "results/palm/ancestral_paths_new-ms-all-palm_report_prs.tsv"
# output_png <- "TB-SNPs.png"
# rsid_list <- c("rs422951", "rs3135363", "rs3806156", "rs2395162", "rs9268530", "rs2596472", "rs2844494")

# TB Meta-GWAS SNPs
description <- "Tuberculosis (TB) - MetaGWAS SNPs"
report_tsv <- "results/palm/ancestral_paths_new-ra-all-palm_report_prs.tsv"
output_png <- "TB-metaGWAS-SNPs.png"
rsid_list <- c("rs2858331", "rs2857106", "rs2856993")

# load the palm report, with all the results
palm <- read_tsv(report_tsv, col_types = cols())

# get the models to plot
snps <- palm %>%
    filter(rsid %in% rsid_list) %>%
    # compose the model prefixes from the PALM metadata
    mutate(prefix = paste0("results/clues/", rsid, "/ancestral_paths_new-", chrom, ":", pos, ":", ancestral_allele, ":", derived_allele, "-", ancestry))

models <- list()
for (i in 1:nrow(snps)) {
    # load the CLUES trajectory
    model <- clues_trajectory(snps[i, ]$rsid, snps[i, ]$ancestry, snps[i, ]$prefix)
    models <- append(models, list(model))
}

# load all the trajectories
traj <- bind_rows(models) %>%
    # and join all the metadata
    inner_join(snps, by = c("rsid", "ancestry"), suffix = c("", ".focal"))

# flip `rs3135391`
traj$flip <- traj$rsid %in% c("rs3135391")

traj <- traj %>%
    mutate(
        # use the `flip` flag to polarize the trajectories
        freq = ifelse(flip, 1.0 - freq, freq),
        # add the focal allele to the SNP label
        snp_label = paste0(rsid, ":", ifelse(flip, ancestral_allele, derived_allele)),
        # add the description
        description = description
    )

# display the ancestries in custom sorted order
traj$ancestry <- factor(traj$ancestry, levels = c("ALL", "WHG", "EHG", "CHG", "ANA"))

gen_time <- 28
max_age <- 13665

# constrain the extent of the plotting
xmin <- round(-max_age / gen_time)
xmax <- max(traj$epoch)
xbreaks <- -seq(-xmax, -xmin, round(2000 / gen_time))
xlabels <- round(xbreaks * gen_time / 1000)

# the approximate time (in generations) during which these ancestries existed as discrete populations
ancestry_epochs <- tibble(
    ancestry = c("ALL", "WHG", "EHG", "CHG", "ANA"),
    start = c(-150, -500, -500, -800, -800),
    finish = c(0, -250, -200, -200, -250)
) %>%
    # don't exceed the range of the plot
    rowwise() %>%
    mutate(start = max(start, xmin), finish = min(finish, xmax))

ggplot(traj) +

    # shade the ancestry epoch
    geom_rect(data = ancestry_epochs, aes(xmin = start, xmax = finish, ymin = 0, ymax = 1), alpha = 0.5, fill = "#F4F4F4") +

    # plot the maximum posterior trajectory
    geom_line(aes(x = epoch, y = freq, color = snp_label, alpha = as.numeric(significant)), cex = 1, na.rm = TRUE) +

    # display as a grid
    facet_grid(description ~ ancestry, labeller = labeller(description = label_wrap_gen())) +

    # set the SNP colors
    scale_color_manual(values = snp_colors) +

    # print the labels
    geom_dl(aes(x = epoch, y = freq, label = snp_label, color = snp_label, alpha = as.numeric(significant)), method = list(dl.trans(x = x + 0.1), "last.qp", cex = 0.8), na.rm = TRUE) +

    # plot non-significant trajectories as transparent
    scale_alpha(range = c(0.3, 1)) +

    # set the axis breaks
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), position = "left", expand = expansion(add = c(0.03, 0.03))) +
    scale_x_continuous(limits = c(xmin, xmax), breaks = xbreaks, labels = xlabels, expand = expansion(add = c(0, 230))) +
    ylab("DAF") +
    xlab("kyr BP") +

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
ggsave(output_png, width = 16, height = 4)
