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
quiet(library(ggrepel)) # v0.9.1

# load the helper functions
source("scripts/clues_utils.R")

# get the command line arguments
p <- arg_parser("Plot the joint PALM models")
p <- add_argument(p, "--biobank", help = "The GWAS biobank", default = "ukbb")
p <- add_argument(p, "--report", help = "PALM multi-trait report", default = "results/palm/ancestral_paths_v3-ms-r0.05-kb250-palm_report_multi_trait.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "results/palm/ancestral_paths_v3-ms-r0.05-kb250-palm_report_multi_trait-ukbb.png")

argv <- parse_args(p)

report <- read_tsv(argv$report, col_types = cols())

data <- report %>%
    filter(trait2_gwas == argv$biobank)

top_ancestry <- report %>% group_by() %>% slice_max(abs(trait1_r)) %>% pull(ancestry)

# set the sort order of the marginal traits based on their R-score in the ancestry with the largest score
data$trait2 <- factor(data$trait2, levels = data %>% filter(ancestry == top_ancestry) %>% arrange(trait1_r) %>% pull(trait2) %>% unique())

# apply a Bonferroni correction for the number of tests
num_traits <- length(unique(data$trait2))
p.bonferroni <- 0.05 / num_traits
r.significant <- qnorm(p.bonferroni / 2, lower.tail = FALSE)

biobanks <- list("ukbb" = "UKBB", "finngen" = "FinnGen")

plt <- data %>%
    ggplot(aes(x = trait2, y = trait1_r, color = ancestry)) +

    # add a solid line at the zero mark
    geom_hline(yintercept = 0, color = "darkgrey") +

    # add dotted lines at the significance threshold
    geom_hline(yintercept = -r.significant, linetype = "dashed", color = "red") +
    geom_hline(yintercept = r.significant, linetype = "dashed", color = "red") +
    geom_point() +
    facet_wrap(~ancestry, ncol = 1) +
    ylab("R-score") +
    xlab(paste0(biobanks[[argv$biobank]], " marginal trait (n=", num_traits, ")")) +

    # set the colour scales
    scale_color_manual(values = ancestry_colors) +

    # basic styling
    theme_bw() +
    theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white", color = "white"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.5, "cm"),
        panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_text(angle = 45, hjust = 1)
    )

# save the plot
ggsave(argv$output, plot = plt, heigh = 9, width = num_traits / 6)
