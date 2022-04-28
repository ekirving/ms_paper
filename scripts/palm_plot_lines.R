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
            logp = -log10(snps[i, ]$p.value),
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

# sort the SNPs by their -log10(p) and the direction of the effect (positive on top)
snp_order <- bind_rows(
    df_ml %>% group_by(rsid) %>% slice_min(epoch) %>% mutate(type = "min_epoch"),
    df_ml %>% group_by(rsid) %>% slice_max(epoch) %>% mutate(type = "max_epoch"),
) %>%
    select(rsid, logp, type, prs_freq) %>%
    pivot_wider(id_cols = c(rsid, logp), names_from = "type", values_from = "prs_freq") %>%
    mutate(direction = ifelse(max_epoch > min_epoch, 1, -1)) %>%
    mutate(signed_logp = direction * logp) %>%
    arrange(desc(signed_logp)) %>%
    select(rsid, signed_logp)

# apply the sort ordering, by setting the factor levels
df_ml$rsid <- factor(df_ml$rsid, levels = snp_order$rsid)

# get the Bonferroni corrected significance threshold
bonferroni <- round(-log10(0.05 / nrow(snps)), 1)

# get the maximum -log10(p.value)
max_logp <- max(df_ml$logp)

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
xmin <- min(-df_ml$epoch)
xmax <- max(-df_ml$epoch)
xbreaks <- seq(xmin, xmax + 1, round(1000 / argv$gen_time))
xlabels <- round(xbreaks * argv$gen_time / 1000)

# set the colour bar breaks, ensuring that the first break is the Bonferroni threshold
bar_breaks <- seq(0, max_logp, bonferroni)
bar_labels <- sprintf("%.1f", bar_breaks)
bar_labels[2] <- paste0(bar_labels[2], "*")

plt <- df_ml %>%
    # plot the heatmap
    ggplot(aes(x = epoch, y = prs_freq, group = rsid)) +

    # plot the maximum likelihood trajectories as stacked lines
    geom_area(aes(fill = logp, color = logp), position = "stack", stat = "identity", outline.type = "full") +
    # geom_line(aes(color = logp), position = "stack") +

    # label the ends of each line
    geom_dl(aes(label = label, color = logp), method = list(dl.trans(x = x + 0.1), "last.qp", cex = 0.8), position = "stack", na.rm = TRUE) +

    # set the axis breaks
    scale_y_continuous(limits = c(0, 0.55), breaks = seq(0, 1, .05), expand = c(0, 0), position = "right") +
    scale_x_continuous(limits = c(-xmax, xmin), breaks = -xbreaks, labels = xlabels, expand = expansion(add = c(1, 55))) +
    labs(
        title = plot_title,
        fill = "-log10(p.value)"
    ) +
    ylab("Scaled PRS") +
    xlab("kyr BP") +

    # set the colour scales
    scale_fill_viridis_c(option = "plasma", breaks = bar_breaks, labels = bar_labels, end = 0.85, alpha = 0.9) +
    scale_color_viridis_c(option = "plasma", end = 0.85) +

    # hide these legends
    guides(color = "none") +

    # basic styling
    theme_bw() +
    theme(
        # legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()
    )

# save the plot
ggsave(filename = argv$output, plt, width = 12, heigh = 8)
