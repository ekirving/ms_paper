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
p <- arg_parser("Plot the trajectory of the polygenic risk score as a posterior density raster")
p <- add_argument(p, "--tsv", help = "PALM report", default = "results/palm/ancestral_paths_new-ALL-ms-palm_report.tsv")
p <- add_argument(p, "--json", help = "PALM json file", default = "results/palm/ancestral_paths_new-ALL-ms-palm.json")
p <- add_argument(p, "--trait", help = "The complex trait name", default = "ms")
p <- add_argument(p, "--dataset", help = "The dataset", default = "ancestral_paths_new")
p <- add_argument(p, "--ancestry", help = "The ancestry path", default = "ALL")
p <- add_argument(p, "--gen-time", help = "Generation time", default = 28)
p <- add_argument(p, "--min-density", help = "Minimum posterior density", default = "1e-20")
p <- add_argument(p, "--output", help = "PALM trajectory", default = "results/palm/ancestral_paths_new-ALL-ms-palm_trajectory.png")

argv <- parse_args(p)

# reduce the memory overhead by filtering out frequency bins below this probability threshold
MIN_POSTERIOR_DENSITY <- as.numeric(argv$min_density)

# reduce the memory overhead by rounding to this maximum precision
MAX_PRECISION <- 5

# load the PALM results
results <- fromJSON(argv$json)

# load the list of SNPs and their effect sizes
snps <- read_tsv(argv$tsv, col_types = cols())

# get the maximum possible PRS
prs_max <- sum(abs(snps$beta))

snps <- snps %>%
    # PALM assumes that betas measure the effect of the ALT allele, and CLUES models the frequency of the derived allele
    # here we want to polarize the trajectories by the positive effect allele
    mutate(flip = (alt == derived_allele & beta < 0) | (alt == ancestral_allele & beta > 0)) %>%
    # normalise the beta scores by the maximum PRS
    mutate(prs = abs(beta) / prs_max) %>%
    # compose the model prefixes from the PALM metadata
    mutate(prefix = paste0("results/clues/", rsid, "/", argv$dataset, "-", chrom, ":", pos, ":", ancestral_allele, ":", derived_allele, "-", argv$ancestry))

models <- list()
for (i in 1:nrow(snps)) {
    # load the CLUES model
    model <- clues_load_data(snps[i, ]$prefix)

    if (snps[i, ]$flip) {
        # polarize the model (if necessary)
        model <- model %>%
            mutate(freq = 1.0 - freq) %>%
            arrange(freq)
    }

    model <- model %>%
        # apply the scaled PRS
        mutate(prs_freq = snps[i, ]$prs * freq)

    models <- append(models, list(model))
}

df_prob <- models[[1]] %>%
    filter(density > 0)

for (model in models[-1]) {
    # get the joint probability of this model and all the previous models
    df_prob <- model %>%
        filter(density > 0) %>%
        inner_join(df_prob, by = "epoch") %>%
        mutate(
            prs_freq = round(prs_freq.x + prs_freq.y, MAX_PRECISION),
            density = density.x * density.y
        ) %>%
        group_by(epoch, prs_freq) %>%
        summarise(density = sum(density), .groups = "drop") %>%
        # prevent combinatorial explosion by dropping extremely low probability points
        filter(density > MIN_POSTERIOR_DENSITY)

    # force garbage collection so we don't cause an out-of-memory crash
    gc()
}

# get the original frequency bins
bins <- models[[1]] %>%
    filter(epoch == 0) %>%
    arrange(freq) %>%
    mutate(bin = row_number()) %>%
    select(bin, freq, height)

# group the summed marginal frequencies into the original frequency bins
df_prob <- df_prob %>%
    mutate(bin = cut(prs_freq, breaks = c(-1, bins$freq), labels = FALSE)) %>%
    group_by(epoch, bin) %>%
    summarise(density = sum(density), .groups = "drop") %>%
    left_join(bins, by = "bin")

# decode the ancestry and trait names
ancestries <- list(
    "ALL" = "All ancestries",
    "ANA" = "Anatolian Farmers",
    "CHG" = "Caucasus Hunter-gatherers",
    "WHG" = "Western Hunter-gatherers",
    "EHG" = "Eastern Hunter-gatherers"
)

traits <- list(
    "ms" = "Multiple sclerosis",
    "ra" = "Rheumatoid arthritis",
    "ibd" = "Inflammatory bowel disease",
    "celiac" = "Celiac disease"
)

# add the model results to the plot title
plot_title <- paste0(
    traits[[argv$trait]],
    " (n = ", results$num_loci, ")",
    " | ", ancestries[[results$ancestry]],
    " | ω = ", results$sel,
    " | se = ", results$se,
    " | z = ", results$z,
    " | p = ", signif(pnorm(q = abs(as.numeric(results$z)), lower.tail = FALSE) * 2, 3)
)

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
    labs(
        title = plot_title,
        fill = "Density"
    ) +
    ylab("Scaled PRS") +
    xlab("kyr BP") +

    # set the colour scales
    scale_fill_viridis_c(option = "plasma") +

    # basic styling
    theme_bw() +
    theme(
        # legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        # fill in all the blanks with the zero density colour
        panel.background = element_rect(fill = "#0D1687", color = NA)
    )

# save the plot
ggsave(filename = argv$output, plt, width = 12, height = 8)
