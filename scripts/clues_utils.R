#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(tidyverse)) # v1.3.1
quiet(library(jsonlite)) # v1.8.0
quiet(library(RcppCNPy)) # v0.2.11
quiet(library(zoo)) # v1.8_10

clues_load_data <- function(prefix) {

    # load the model data
    epochs <- npyLoad(paste0(prefix, ".epochs.npy"))
    freqs <- npyLoad(paste0(prefix, ".freqs.npy"))
    logpost <- npyLoad(paste0(prefix, ".post.npy"), dotranspose = F)

    # add column names
    colnames(logpost) <- paste0("V", seq(ncol(logpost)))

    model <- as_tibble(logpost) %>%
        # convert posterior densities from log-likelihoods
        exp() %>%
        # add the frequency labels
        add_column(freq = freqs, .before = 2) %>%
        # add the title heights (and a little padding to make sure there are no gaps)
        # tile heights are not equal because there is higher sampling density near 0 and 1
        add_column(height = diff(c(0, freqs)) + 1e-4, .before = 3) %>%
        # pivot the columns into long format
        pivot_longer(-c(freq, height), names_to = "epoch", values_to = "density") %>%
        # convert the column names into epochs (and switch the direction of time)
        mutate(epoch = -epochs[as.numeric(gsub("V", "", epoch))])

    model
}

clues_trajectory <- function(rsid, ancestry, prefix, smooth = 10, ancestral = NULL, threshold = NA) {

    # load the model
    model <- clues_load_data(prefix) %>% mutate(rsid = rsid, ancestry = ancestry)

    # extract the maximum posterior trajectory
    traj <- model %>%
        group_by(epoch) %>%
        top_n(1, density) %>%
        ungroup() %>%
        arrange(epoch)

    # trim the portion of the trajectory with low posterior density
    if (!is.na(threshold)) {
        # find the most recent epoch which is below the threshold
        cutoff <- traj %>%
            filter(density < threshold) %>%
            group_by() %>%
            slice_max(epoch)

        if (length(cutoff$epoch)) {
            # drop everything before that epoch
            traj <- traj %>%
                filter(epoch > cutoff$epoch)
        }
    }

    # apply a little smoothing to the jagged steps in the trajectory
    if (smooth) {
        traj$freq <- rollapply(c(traj$freq, rep_len(NA, smooth - 1)), width = smooth, by = 1, FUN = mean, na.rm = TRUE, align = "left")
    }

    if (!is.null(ancestral)) {
        metadata <- fromJSON(file = paste0("variants/metadata/GRCh37/", rsid, ".json"))

        # re-polarise, if necessary
        if (metadata$ancestral != ancestral) {
            traj$freq <- 1 - traj$freq
        }
    }

    traj
}

traits <- list(
    "cd-all" = "Celiac disease (all genome-wide significant)",
    "cd-r0.05-kb250" = "Celiac disease (r^2 < 0.05; window 250 kb)",
    "ms-all" = "Multiple sclerosis (all genome-wide significant)",
    "ms-r0.05-kb250" = "Multiple sclerosis (r^2 < 0.05; window 250 kb)",
    "ra-all" = "Rheumatoid arthritis (all genome-wide significant)",
    "ra-r0.05-kb250" = "Rheumatoid arthritis (r^2 < 0.05; window 250 kb)"
)

# decode the ancestry and trait names
ancestries <- list(
    "ALL" = "All ancestries",
    "ANA" = "Anatolian Farmers",
    "CHG" = "Caucasus Hunter-gatherers",
    "WHG" = "Western Hunter-gatherers",
    "EHG" = "Eastern Hunter-gatherers"
)

ancestry_colors <- c(
    "ALL" = "#66c2a5",
    "WHG" = "#fc8d62",
    "EHG" = "#8da0cb",
    "CHG" = "#e78ac3",
    "ANA" = "#a6d854"
)

# color brewer https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=10
snp_colors <- c(
    "#1f78b4",
    "#33a02c",
    "#e31a1c",
    "#ff7f00",
    "#6a3d9a",
    "#a6cee3",
    "#b2df8a",
    "#fb9a99",
    "#fdbf6f",
    "#cab2d6"
)
