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
p <- arg_parser("Plot the trajectory of the polygenic risk score as a stacked line")
p <- add_argument(p, "--palm", help = "PALM metadata", default = "data/targets/all_clumped_annotated_ms_ancestral_paths_new_palm.tsv")
p <- add_argument(p, "--json", help = "PALM json file", default = "results/palm/ancestral_paths_new/ALL/ms/ms_palm.json")
p <- add_argument(p, "--trait", help = "The complex trait name", default = "ms")
p <- add_argument(p, "--dataset", help = "The dataset", default = "ancestral_paths_new")
p <- add_argument(p, "--ancestry", help = "The ancestry path", default = "ALL")
p <- add_argument(p, "--gen-time", help = "Generation time", default = 28)
p <- add_argument(p, "--output", help = "PALM trajectory", default = "results/palm/ancestral_paths_new/ALL/ms/ms_palm.png")

argv <- parse_args(p)

# load the PALM results
results <- fromJSON(argv$json)

# load the list of SNPs and their effect sizes
palm <- read_tsv(argv$palm, col_types = cols())

# TODO remove when done testing
palm <- head(palm, n = 10)

# get the maximum possible PRS
prs_max <- sum(abs(palm$beta))

palm <- palm %>%
    # PALM assumes that betas measure the effect of the ALT allele, and CLUES models the frequency of the derived allele
    # here we want to polarize the trajectories by the positive effect allele
    mutate(flip = (alt == derived_allele & beta < 0) | (alt == ancestral_allele & beta > 0)) %>%
    # normalise the beta scores by the maximum PRS
    mutate(prs = abs(beta) / prs_max) %>%
    # compose the model prefixes from the PALM metadata
    mutate(prefix = paste0("results/clues/", rsid, "/", argv$dataset, "-", chrom, ":", pos, ":", ancestral_allele, ":", derived_allele, "-", argv$ancestry))

models <- list()
for (i in 1:nrow(palm)) {
    # load the CLUES model, and extract the maximum likelihood path
    model <- clues_load_data(palm[i, ]$prefix) %>%
        group_by(epoch) %>%
        top_n(1, density) %>%
        ungroup() %>%
        arrange(epoch) %>%
        mutate(prefix=palm[i, ]$prefix)

    if (palm[i, ]$flip) {
        # polarize the model (if necessary)
        model <- model %>%
            mutate(freq = 1.0 - freq)
    }

    model <- model %>%
        # apply the scaled PRS
        mutate(prs_freq = palm[i, ]$prs * freq)

    models <- append(models, list(model))
}

df_ml <- bind_rows(models)

# constrain the extent of the plotting
xmin <- min(-df_ml$epoch)
xmax <- max(-df_ml$epoch)
xbreaks <- seq(xmin, xmax + 1, round(1000 / argv$gen_time))
xlabels <- round(xbreaks * argv$gen_time / 1000)

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

df_ml %>%
    # plot the heatmap
    ggplot(aes(x = epoch, y = prs_freq, color = prefix, fill = prefix)) +

    # plot the maximum likelihood trajectories as stacked lines
    geom_line(position = "stack", size = 2) +
    geom_area(position="stack", stat="identity",alpha = 0.5) +

    # set the axis breaks
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), expand = c(0, 0), position = "right") +
    scale_x_continuous(limits = c(-xmax, xmin), breaks = -xbreaks, labels = xlabels, expand = c(0, 0)) +
    labs(title = plot_title) +
    ylab("Scaled PRS") +
    xlab("kyr BP") +

    # basic styling
    theme_minimal() +
    theme(
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
    )

# save the plot
ggsave(filename = argv$output, plt, width=12, heigh=8)
