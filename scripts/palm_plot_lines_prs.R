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
p <- add_argument(p, "--tsv", help = "PALM report", default = "results/palm/ancestral_paths_v3-ALL-ms-r0.05-kb250-palm_report.tsv")
p <- add_argument(p, "--json", help = "PALM json file", default = "results/palm/ancestral_paths_v3-ALL-ms-r0.05-kb250-palm.json")
p <- add_argument(p, "--trait", help = "The complex trait name", default = "ms-r0.05-kb250")
p <- add_argument(p, "--dataset", help = "The dataset", default = "ancestral_paths_v3")
p <- add_argument(p, "--ancestry", help = "The ancestry path", default = "ALL")
p <- add_argument(p, "--gen-time", help = "Generation time", default = 28)
p <- add_argument(p, "--output", help = "PALM trajectory", default = "results/palm/ancestral_paths_v3-ALL-ms-palm_lines-prs.png")

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
    # record the focal allele (to avoid any ambiguity)
    mutate(focal_allele = ifelse(beta < 0, ref, alt), .after = "derived_allele") %>%
    # normalise the beta scores by the maximum PRS
    mutate(scaled_beta = abs(beta) / prs_max) %>%
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
            significant = snps[i, ]$significant,
            chrom = snps[i, ]$chrom,
            pos = snps[i, ]$pos
        )

    if (snps[i, ]$flip) {
        # polarize the model (if necessary)
        model <- model %>%
            mutate(freq = 1.0 - freq)
    }

    model <- model %>%
        # apply the scaled PRS
        mutate(prs_freq = snps[i, ]$scaled_beta * freq)

    models <- append(models, list(model))
}

df_ml <- bind_rows(models)

# sort the SNPs by the magnitude of their change in PRS
snp_order <- bind_rows(
    df_ml %>% group_by(rsid) %>% slice_min(epoch) %>% mutate(name = "prs_start"),
    df_ml %>% group_by(rsid) %>% slice_max(epoch) %>% mutate(name = "prs_end"),
) %>%
    select(rsid, name, logp, prs_freq) %>%
    pivot_wider(names_from = "name", values_from = "prs_freq") %>%
    mutate(delta_prs = prs_end - prs_start) %>%
    arrange(desc(delta_prs)) %>%
    select(rsid, delta_prs)

# define a delta PRS threshold for labeling SNPs
delta_prs_threshold <- 0.0025

# join the delta_prs column so we can plot it as a fill colour
df_ml <- df_ml %>%
    inner_join(snp_order, by = "rsid") %>%
    # only label the SNPs above the threshold
    mutate(label = ifelse(abs(delta_prs) > delta_prs_threshold, ifelse(chrom == 6 & pos >= 28477797 & pos <= 33448354, paste0(rsid, "†"), rsid), NA))

# apply the sort ordering, by setting the factor levels
df_ml$rsid <- factor(df_ml$rsid, levels = snp_order$rsid)

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
# some SNPs have trajectories that don't go back the as far as others
limits <- df_ml %>%
    group_by(rsid) %>%
    summarise(xmin = min(-epoch), xmax = max(-epoch)) %>%
    group_by() %>%
    summarise(xmin = min(xmin), xmax = min(xmax))

xbreaks <- seq(limits$xmin, limits$xmax + 1, round(1000 / argv$gen_time))
xlabels <- round(xbreaks * argv$gen_time / 1000)

# get the maximum scaled PRS
max_prs <- df_ml %>%
    group_by(epoch) %>%
    summarise(prs = sum(prs_freq)) %>%
    pull(prs) %>%
    max()

# determine a sensible y-axis limit
ylimit <- ceiling(max_prs * 10) / 10

# set the range of the colorbar breaks
min_break <- round(min(df_ml$delta_prs) / delta_prs_threshold) * delta_prs_threshold
max_break <- round(max(df_ml$delta_prs) / delta_prs_threshold) * delta_prs_threshold

# set the colour bar breaks, ensuring that the first break is the delta_prs_threshold threshold
bar_breaks <- seq(min(min_break, 0), max(max_break, delta_prs_threshold), delta_prs_threshold)
bar_labels <- sprintf("%.4f", bar_breaks)

# define a rescaling function which caps the upper range of the colorbar at the delta_prs threshold
show_significant <- function(x, to = c(0, 1), from = NULL) {
    ifelse(x < delta_prs_threshold, scales::rescale(x, to = to, from = c(min(x, na.rm = TRUE), delta_prs_threshold)), 1)
}

plt <- df_ml %>%
    # plot the heatmap
    ggplot(aes(x = epoch, y = prs_freq, group = rsid)) +

    # plot the maximum likelihood trajectories as stacked lines
    geom_area(aes(fill = delta_prs, color = delta_prs), position = "stack", stat = "identity", outline.type = "full") +

    # label the ends of each line
    geom_dl(aes(label = label, color = delta_prs), method = list(dl.trans(x = x + 0.1), "last.qp", cex = 0.8), position = "stack", na.rm = TRUE) +

    # set the axis breaks
    scale_y_continuous(limits = c(0, ylimit), breaks = seq(0, 1, .05), expand = c(0, 0), position = "right") +
    scale_x_continuous(limits = c(-limits$xmax, limits$xmin), breaks = -xbreaks, labels = xlabels, expand = expansion(add = c(1, 80))) +
    labs(
        title = plot_title,
        fill = "Δ PRS"
    ) +
    ylab("Scaled PRS") +
    xlab("kyr BP") +

    # set the colour scales
    scale_fill_viridis_c(option = "plasma", rescaler = show_significant, breaks = bar_breaks, labels = bar_labels, end = 0.85, alpha = 0.9) +
    scale_color_viridis_c(option = "plasma", rescaler = show_significant, end = 0.85) +

    # hide these legends
    guides(color = "none") +

    # basic styling
    theme_bw() +
    theme(
        # legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size = 11)
    )

# save the plot
ggsave(filename = argv$output, plt, width = 12, height = 8)
