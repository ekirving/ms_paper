#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(tidyverse)) # v1.3.1
quiet(library(RcppCNPy)) # v0.2.11

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

traits <- list(
    "celiac-UKB-r0.05-kb250" = "Celiac disease (R2 < 0.05; window 250 kb)",
    "celiac-UKB-r0.1-kb1000" = "Celiac disease (R2 < 0.1; window 1000 kb)",
    "ms-r0.05-kb250" = "Multiple sclerosis (R2 < 0.05; window 250 kb)",
    "ms-r0.1-kb1000" = "Multiple sclerosis (R2 < 0.1; window 1000 kb)",
    "ra-r0.05-kb250" = "Rheumatoid arthritis (R2 < 0.05; window 250 kb)",
    "ra-r0.1-kb1000" = "Rheumatoid arthritis (R2 < 0.1; window 1000 kb)",
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
