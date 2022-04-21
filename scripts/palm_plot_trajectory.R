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
p <- arg_parser("Convert the GWAS metadata into PALM input format")
p <- add_argument(p, "--palm", help = "PALM metadata", default = "data/targets/all_clumped_annotated_ms_ancestral_paths_new_palm.tsv")
p <- add_argument(p, "--trait", help = "The complex trait name", default = "Multiple sclerosis")
p <- add_argument(p, "--datasource", help = "The datasource", default = "ancestral_paths_new")
p <- add_argument(p, "--ancestry", help = "The ancestry path", default = "ALL")
p <- add_argument(p, "--gen-time", help = "Generation time", default = 28)
p <- add_argument(p, "--output", help = "PALM trajectory", default = "data/targets/all_clumped_annotated_ms_ancestral_paths_new_palm.png")

argv <- parse_args(p)

# load the list of SNPs and their effect sizes
palm <- read_tsv(argv$palm, col_types = cols())

# compose the model prefixes from the PALM metadata
prefixes <- palm %>%
    mutate(prefix = paste0("results/clues/", rsid, "/", argv$datasource, "-", chrom, ":", pos, ":", ancestral_allele, ":", derived_allele, "-", argv$ancestry)) %>%
    pull(prefix)

# load all the models
models <- lapply(prefixes, clues_load_data)

df_prob <- models[[1]]

for (model in models[-1]) {
    # get the joint probability of this model and the previous models
    df_prob <- inner_join(df_prob, model, by = "epoch") %>%
        mutate(
            freq = freq.x + freq.y,
            density = density.x * density.y
        ) %>%
        group_by(epoch, freq) %>%
        summarise(density = sum(density), .groups = "drop") %>% 
        # prevent combinatorial explosion by dropping extremely low probability points
        filter(density > 1e-20)

    # force garbage collection so we don't cause an out-of-memory crash
    gc()
}

# get the original frequency bins
bins <- models[[1]] %>%
    filter(epoch == 0) %>%
    mutate(bin = row_number()) %>%
    select(bin, freq, height)

# group the summed marginal frequencies into the original frequency bins
df_prob <- df_prob %>%
    mutate(freq = freq / length(models)) %>%
    mutate(bin = cut(freq, breaks = c(-1, bins$freq), labels = FALSE)) %>%
    group_by(epoch, bin) %>%
    summarise(density = sum(density), .groups = "drop") %>%
    left_join(bins, by = "bin")

# constrain the extent of the plotting
xmin <- min(-df_prob$epoch)
xmax <- max(-df_prob$epoch)
xbreaks <- seq(xmin, xmax + 1, round(1000 / argv$gen_time))
xlabels <- round(xbreaks * argv$gen_time / 1000)

# get the maximum likelihood trajectory
max_traj <- df_prob %>%
    group_by(epoch) %>%
    top_n(1, density) %>%
    ungroup() %>%
    arrange(epoch)

plt <- df_prob %>%
    # plot the heatmap
    ggplot(aes(x = epoch, y = freq, height = height, fill = density)) +
    geom_tile() +

    # plot the maximum likelihood trajectory
    # geom_line(mapping = aes(x = epoch, y = freq), data = max_traj, color = "red") +

    # set the axis breaks
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), expand = c(0, 0), position = "right") +
    scale_x_continuous(limits = c(-xmax, xmin), breaks = -xbreaks, labels = xlabels, expand = c(0, 0)) +
    scale_fill_viridis_c(limits = c(0, 0.5), option = "plasma") +
    labs(title = argv$trait, fill = "Density") +
    ylab("Scaled PRS") +
    xlab("kyr BP") +

    # basic styling
    theme_minimal() +
    theme(
        # legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        # fill in all the blanks with the zero density colour
        panel.background = element_rect(fill="#0D1687", color=NA)
    )

png(file = argv$output, width = 16, height = 9, units = "in", res = 300)
plt
dev <- dev.off()

plt