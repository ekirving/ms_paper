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
quiet(library(RcppCNPy)) # v0.2.11
quiet(library(RcppRoll)) # v0.3.0

# get the command line arguments
p <- arg_parser("Convert the GWAS metadata into PALM input format")
p <- add_argument(p, "--palm", help = "PALM metadata", default = "data/targets/all_clumped_annotated_ms_ancestral_paths_new_palm.tsv")
p <- add_argument(p, "--gen-time", help = "Generation time", default = 28)
p <- add_argument(p, "--output", help = "PALM trajectory", default = "data/targets/all_clumped_annotated_ms_ancestral_paths_new_palm.png")

argv <- parse_args(p)

# load the list of SNPs and their effect sizes
palm <- read_tsv(argv$palm, col_types = cols())

clues_trajectory <- function(prefix) {

    # load the CLUES model data
    epochs <- npyLoad(paste0(prefix, ".epochs.npy"))
    freqs <- npyLoad(paste0(prefix, ".freqs.npy"))
    logpost <- npyLoad(paste0(prefix, ".post.npy"), dotranspose = F)

    # .name_repair = make.names
    df <- as_tibble(logpost) %>%
        # convert posterior densities from log-likelihoods
        exp() %>%
        # add the frequency labels
        add_column(freq = freqs, .before = 1) %>%
        # add the title heights (and a little padding to make sure there are no gaps)
        # tile heights are not equal because there is higher sampling density near 0 and 1
        add_column(height = diff(c(0, freqs)) + 1e-4, .before = 2) %>%
        # pivot the columns into long format
        pivot_longer(-c(freq, height), names_to = "epoch", values_to = "density") %>%
        # convert the column names into epochs (and switch the direction of time)
        mutate(epoch = -epochs[as.numeric(gsub("V", "", epoch))])

    # extract the maximum likelihood trajectory
    traj <- df %>%
        group_by(epoch) %>%
        top_n(1, density) %>%
        ungroup() %>%
        arrange(epoch) %>%
        mutate(freq = roll_mean(freq, 5, align = "left", fill = max(freq)))

    traj
}

# TODO remove when done testing
models <- c(
    # "results/clues/rs10175798/ancestral_paths_new-2:30449594:G:A-ALL",
    "results/clues/rs12527959/ancestral_paths_new-6:29537426:C:T-ALL"
)

# load the trajectories
traj <- bind_rows(
    lapply(models, function(model) {
        clues_trajectory(model) %>% mutate(label = model)
    })
)

# constrain the extent of the plotting
xmin <- min(-traj$epoch)
xmax <- max(-traj$epoch)
xbreaks <- seq(xmin, xmax + 1, round(1000 / argv$gen_time))
xlabels <- round(xbreaks * argv$gen_time / 1000)

plt <- ggplot(traj) +

    # plot the density raster
    # geom_tile(aes(x = epoch, y = freq, height = height, fill = density)) +

    # plot the maximum posterior trajectory
    geom_line(aes(x = epoch, y = freq, color = label)) +

    # set the axis breaks
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), expand = c(0, 0), position = "right") +
    scale_x_continuous(limits = c(-xmax, xmin), breaks = -xbreaks, labels = xlabels, expand = c(0, 0)) +
    scale_fill_viridis_c(limits = c(0, 0.5)) +
    labs(title = "Trait name") +
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

plt
