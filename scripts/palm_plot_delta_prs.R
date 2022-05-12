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

# get the command line arguments
p <- arg_parser("Plot the density of the delta PRS, stratified by ancestry")
p <- add_argument(p, "--trait", help = "The complex trait name", default = "ms")
p <- add_argument(p, "--dataset", help = "The dataset", default = "ancestral_paths_new")
p <- add_argument(p, "--all-json", help = "PALM json file for ALL", default = "results/palm/ancestral_paths_new-ALL-ms-palm.json")
p <- add_argument(p, "--ana-json", help = "PALM json file for ANA", default = "results/palm/ancestral_paths_new-ANA-ms-palm.json")
p <- add_argument(p, "--chg-json", help = "PALM json file for CHG", default = "results/palm/ancestral_paths_new-CHG-ms-palm.json")
p <- add_argument(p, "--ehg-json", help = "PALM json file for EHG", default = "results/palm/ancestral_paths_new-EHG-ms-palm.json")
p <- add_argument(p, "--whg-json", help = "PALM json file for WHG", default = "results/palm/ancestral_paths_new-WHG-ms-palm.json")
p <- add_argument(p, "--all-tsv", help = "PALM report for ALL", default = "results/palm/ancestral_paths_new-ALL-ms-palm_report_prs.tsv")
p <- add_argument(p, "--ana-tsv", help = "PALM report for ANA", default = "results/palm/ancestral_paths_new-ANA-ms-palm_report_prs.tsv")
p <- add_argument(p, "--chg-tsv", help = "PALM report for CHG", default = "results/palm/ancestral_paths_new-CHG-ms-palm_report_prs.tsv")
p <- add_argument(p, "--ehg-tsv", help = "PALM report for EHG", default = "results/palm/ancestral_paths_new-EHG-ms-palm_report_prs.tsv")
p <- add_argument(p, "--whg-tsv", help = "PALM report for WHG", default = "results/palm/ancestral_paths_new-WHG-ms-palm_report_prs.tsv")
p <- add_argument(p, "--output", help = "PALM trajectory", default = "results/palm/ancestral_paths_new-ms-delta_prs.png")

argv <- parse_args(p)

# load all the PALM results
results <- bind_rows(
    fromJSON(argv$all_json),
    fromJSON(argv$ana_json),
    fromJSON(argv$chg_json),
    fromJSON(argv$ehg_json),
    fromJSON(argv$whg_json),
)

# load the list of SNPs and their effect sizes
snps <- bind_rows(
    read_tsv(argv$all_tsv, col_types = cols()),
    read_tsv(argv$ana_tsv, col_types = cols()),
    read_tsv(argv$chg_tsv, col_types = cols()),
    read_tsv(argv$ehg_tsv, col_types = cols()),
    read_tsv(argv$whg_tsv, col_types = cols())
)

# set the sort order of the ancestries based on their Z-score
snps$ancestry <- factor(snps$ancestry, levels = results %>% arrange(z) %>% pull(ancestry))

traits <- list(
    "ms" = "Multiple sclerosis",
    "ra" = "Rheumatoid arthritis",
    "ibd" = "Inflammatory bowel disease",
    "celiac" = "Celiac disease"
)

# add the model results to the plot title
plot_title <- traits[[argv$trait]]

ancestry_colors <- c(
    "ALL" = "#66c2a5",
    "WHG" = "#fc8d62",
    "EHG" = "#8da0cb",
    "CHG" = "#e78ac3",
    "ANA" = "#a6d854"
)

plt <- ggplot(snps, aes(x = ancestry, y = delta_prs)) + 
    # add a red line at the zero mark
    geom_hline(yintercept=0, color="red") +
    
    # density plot of the distribution
    ggdist::stat_halfeye(aes(fill=ancestry), adjust = .5, width = .6, .width = 0, justification = -.3, point_colour = NA) + 
    
    # boxplot of the distribution 
    geom_boxplot(width = .25, outlier.shape = NA) +
    
    # jittered scatter plot of the distribution
    geom_point(aes(colour=ancestry), size = 1.3, alpha = .3, position = position_jitter(seed = 1, width = .03)) +
    
    # flip the plot from vertical to horizontal
    coord_flip() +

    # ancestry_colors
    
    labs(
        title = plot_title,
        fill = "Ancestry"
    ) +
    ylab("Δ PRS per SNP") +
    xlab("Ancestry painting") +
    
    # set the colour scales
    scale_fill_manual(values=ancestry_colors) +
    scale_color_manual(values=ancestry_colors) +
    
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
ggsave(filename = argv$output, plt, width = 12, height = 8)
