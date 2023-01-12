#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(tidyverse))

pathogen <- "brucella"
# pathogen <- "listeria"

# get the sample metadata
samples <- read_tsv("data/Ancestral_paths_new/ancestral_paths_merged_filtered_age.sampleInfo.tsv", col_types = cols()) %>%
    select(sampleId, popId, groupAge, ageAverage)

BIN_SIZE <- 1000

# get the pathogen presence/absence data
hits <- read_tsv(paste0("data/pathogens/presence_absence_", pathogen, ".tsv"), col_types = cols()) %>%
    rename(sampleId = Original_sample) %>%
    inner_join(samples, by = "sampleId")

binned <- hits %>%
    pivot_longer(cols = -c("sampleId", "popId", "groupAge", "ageAverage"), names_to = "species", values_to = "present") %>%
    # add a new age bin
    mutate(bin = -ceiling(ageAverage / BIN_SIZE) * BIN_SIZE) %>%
    group_by(species, bin) %>%
    summarise(sum = sum(present), count = n(), freq = sum / count, .groups = "drop_last")

# constrain the extent of the plotting
xmin <- min(binned$bin)
xmax <- max(binned$bin)
xbreaks <- seq(xmin, xmax, BIN_SIZE)
xlabels <- round(xbreaks / BIN_SIZE)

plt <- ggplot(binned, aes(x = bin, y = freq)) +
    geom_point() +
    geom_smooth(method = "loess", formula = "y ~ x", fullrange = TRUE, se = FALSE) +
    facet_wrap(~species, ncol = 3) +
    scale_x_continuous(breaks = xbreaks, labels = xlabels) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    ylab("Frequency") +
    xlab("kyr BP") +

    # basic styling
    theme_minimal() +
    theme(
        # legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white", color = "white"),
        panel.spacing = unit(1, "lines")
    )

height <- (binned %>% pull(species) %>% unique() %>% length()) / 3 * 1.5

ggsave(paste0("figs/time-series-", pathogen, ".png"), plot = plt, width = 8, height = height, units = "in")
