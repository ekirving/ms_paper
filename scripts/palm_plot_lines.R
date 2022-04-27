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
quiet(library(directlabels)) # v2021.1.13

# load the helper functions
source("scripts/clues_utils.R")

# get the command line arguments
p <- arg_parser("Plot the trajectory of the polygenic risk score as stacked lines")
p <- add_argument(p, "--tsv", help = "PALM report", default = "results/palm/ancestral_paths_new/ALL/ms/ms_palm_report.tsv")
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
    # load the CLUES model, and extract the maximum likelihood path
    model <- clues_load_data(snps[i, ]$prefix) %>%
        group_by(epoch) %>%
        top_n(1, density) %>%
        ungroup() %>%
        arrange(epoch) %>%
        mutate(
            rsid = snps[i, ]$rsid, 
            significant = snps[i, ]$significant
        )

    if (snps[i, ]$flip) {
        # polarize the model (if necessary)
        model <- model %>%
            mutate(freq = 1.0 - freq)
    }

    model <- model %>%
        # apply the scaled PRS
        mutate(prs_freq = snps[i, ]$prs * freq)

    models <- append(models, list(model))
}

df_ml <- bind_rows(models) %>%
    # only label the significant SNPs
    mutate(label = ifelse(significant, rsid, NA))

# sort the SNPs by the increase in their scaled PRS
snp_order <- df_ml %>%
    filter(epoch %in% c(0, -528)) %>%
    select(rsid, epoch, prs_freq) %>%
    pivot_wider(rsid, names_from="epoch", values_from="prs_freq") %>%
    mutate(effect=`0` - `-528`) %>%
    arrange(desc(effect)) %>%
    pull(rsid)

# apply the sort ordering
df_ml$rsid <- factor(df_ml$rsid, levels = snp_order)

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

plt <- df_ml %>%
    # plot the heatmap
    ggplot(aes(x = epoch, y = prs_freq, color = rsid)) +

    # plot the maximum likelihood trajectories as stacked lines
    geom_line(position = "stack") +
    geom_area(aes(fill = rsid), position = "stack", stat = "identity", alpha = 0.5) +

    # label the ends of each line
    geom_dl(aes(label = label), method = list(dl.trans(x = x + 0.1), "last.qp", cex = 0.8), position = "stack", na.rm = TRUE) +

    # set the axis breaks
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .1), expand = c(0, 0), position = "right") +
    scale_x_continuous(limits = c(-xmax, xmin), breaks = -xbreaks, labels = xlabels, expand = expansion(add = c(0, 55))) +
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
