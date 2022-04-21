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
quiet(library(RcppCNPy))  # v0.2.11
quiet(library(RcppRoll))  # v0.3.0

# load the helper functions
source("scripts/clues_utils.R")

# get the command line arguments
p <- arg_parser("Convert the GWAS metadata into PALM input format")
p <- add_argument(p, "--palm", help = "PALM metadata", default = "data/targets/all_clumped_annotated_ms_ancestral_paths_new_palm.tsv")
p <- add_argument(p, "--gen-time", help = "Generation time", default = 28)
p <- add_argument(p, "--output", help = "PALM trajectory", default = "data/targets/all_clumped_annotated_ms_ancestral_paths_new_palm.png")

argv <- parse_args(p)

# load the list of SNPs and their effect sizes
palm <- read_tsv(argv$palm, col_types = cols())

prefixes <- c(
    "results/clues/rs10175798/ancestral_paths_new-2:30449594:G:A-ALL",
    "results/clues/rs12527959/ancestral_paths_new-6:29537426:C:T-ALL"
)

models <- list(
    clues_load_data(prefixes[[1]], "rs10175798"),
    clues_load_data(prefixes[[2]], "rs12527959")
)

df_prob <- models[[1]]

for (model in models[-1]) {
    # model <- models[-1][[1]]
    df_prob <- inner_join(df_prob, model, by = 'epoch') %>%
        mutate(
            freq=freq.x+freq.y,
            # TODO float underflow issue
            density=density.x*density.y
        ) %>%
        group_by(epoch, freq) %>%
        summarise(density=sum(density), .groups="drop")
}

# get the original frequency bins
bins <- models[[1]] %>% 
    filter(epoch==0) %>% 
    mutate(bin=row_number()) %>%
    select(bin, freq, height)

# convert the summed marginal frequencies back to the original frequency bins
df_prob <- df_prob %>%
    mutate(freq=freq/length(models)) %>%
    mutate(bin=cut(df_prob$freq, breaks=c(-1, bins$freq), labels=FALSE)) %>%
    group_by(epoch, bin) %>%
    summarise(density=sum(density), .groups="drop") %>%
    left_join(bins, by="bin")

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
    arrange(epoch) %>%
    mutate(freq = roll_mean(freq, 5, align = "left", fill = max(freq)))

df_prob %>%
    # plot the heatmap
    ggplot(aes(x = epoch, y = freq, height = height, fill = density)) +
    geom_tile() +
    
    # plot the maximum posterior trajectory
    geom_line(aes(x = epoch, y = freq), max_traj) +
    
    # set the axis breaks
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), expand = c(0, 0), position = "right") +
    scale_x_continuous(limits = c(-xmax, xmin), breaks = -xbreaks, labels = xlabels, expand = c(0, 0)) +
    scale_fill_viridis_c(limits = c(0, 0.5)) +
    # scale_fill_gradient(low = "white", high = "black", limits = c(0, 0.5)) +
    
    labs(title = argv$prefix) +
    ylab("Derived Allele Frequency") +
    xlab("kyr BP") +
    
    # basic styling
    theme_minimal() +
    theme(
        # legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
    )

