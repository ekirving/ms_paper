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

# load the helper functions
source("scripts/clues_utils.R")

# get the command line arguments
p <- arg_parser("Plot the density of the sampling and ancestry painting over time")
p <- add_argument(p, "--tsv", help = "PALM report", default = "results/palm/ancestral_paths_new-ALL-ms-palm_report.tsv")
p <- add_argument(p, "--json", help = "PALM json file", default = "results/palm/ancestral_paths_new-ALL-ms-palm.json")
p <- add_argument(p, "--trait", help = "The complex trait name", default = "ms")
p <- add_argument(p, "--dataset", help = "The dataset", default = "ancestral_paths_new")
p <- add_argument(p, "--gen-time", help = "Generation time", default = 28)
p <- add_argument(p, "--output", help = "PALM trajectory", default = "results/palm/ancestral_paths_new-ALL-ms-palm_trajectory.png")

argv <- parse_args(p)

# load the PALM results
results <- fromJSON(argv$json)

# load the list of SNPs and their effect sizes
snps <- read_tsv(argv$tsv, col_types = cols())

snps <- snps %>%
    # compose the model prefixes from the PALM metadata
    mutate(prefix = paste0("results/clues/", rsid, "/", argv$dataset, "-", chrom, ":", pos, ":", ancestral_allele, ":", derived_allele))

snps <- head(snps, n = 5)


models <- list()
for (i in 1:nrow(snps)) {
    # load the diploid pan-ancestry models
    model <- read_table(paste0(snps[i, ]$prefix, "-ALL.ancient"), col_names = c("generations", "hom_ref", "het", "hom_alt"), col_types = cols()) %>%
        mutate(ancestry = "ALL", rsid = snps[i, ]$rsid) %>%
        select(generations, ancestry, rsid)

    # double the count of ALL ancestries as they are diploid and the others are haploid
    model <- bind_rows(model, model)

    models <- append(models, list(model))
}

for (i in 1:nrow(snps)) {
    for (ancestry in c("ANA", "CHG", "EHG", "WHG")) {
        # load the CLUES genotype files
        model <- read_table(paste0(snps[i, ]$prefix, "-", ancestry, ".ancient"), col_names = c("generations", "hom_ref", "hom_alt"), col_types = cols()) %>%
            mutate(ancestry = ancestry, rsid = snps[i, ]$rsid) %>%
            select(generations, ancestry, rsid)

        models <- append(models, list(model))
    }
}

# merge all the data together
samples <- bind_rows(models) %>%
    mutate(age = round(generations * argv$gen_time))

plot_title <- traits[[argv$trait]]

# plt <-
ggplot(samples, aes(x = ancestry, y = age)) +
    # add a red line at the zero mark
    geom_hline(yintercept = 0, color = "red") +

    # density plot of the distribution
    ggdist::stat_halfeye(aes(fill = ancestry), adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA) +

    # boxplot of the distribution
    geom_boxplot(width = .25, outlier.shape = NA) +

    # jittered scatter plot of the distribution
    # geom_point(aes(colour=ancestry), size = 1.3, alpha = .3, position = position_jitter(seed = 1, width = .03)) +
    geom_point(aes(colour = ancestry), shape = 95, size = 10, alpha = .2) +

    # flip the plot from vertical to horizontal
    coord_flip() +

    # ancestry_colors

    labs(
        title = plot_title,
        fill = "Ancestry"
    ) +
    ylab("Δ PRS per SNP") +
    xlab("Ancestry painting") +

    # set the colour scales
    scale_fill_manual(values = ancestry_colors) +
    scale_color_manual(values = ancestry_colors) +

    # hide these legends
    guides(color = "none") +

    # basic styling
    theme_bw() +
    theme(
        # legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
    )

plt

# save the plot
ggsave(filename = argv$output, plt, width = 12, height = 8)
