#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(tidyverse)) # v1.3.1

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
