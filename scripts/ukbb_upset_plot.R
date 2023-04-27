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
quiet(library(directlabels)) # v2021.1.13
quiet(library(R.utils)) # v2.12.0
quiet(library(ggupset)) # v0.3.0

# load the helper functions
source("scripts/clues_utils.R")

# get the command line arguments
p <- arg_parser("Compare associations between trait associated SNPs and all other traits in UKBB")
p <- add_argument(p, "--trait", help = "Focal trait name", default = "ms-r0.05-kb250")
p <- add_argument(p, "--ukbb", help = "UKBB associations", default = "results/compare/ukbb/ancestral_paths_v3-ms-r0.05-kb250.both_sexes.significant.tsv.gz")
p <- add_argument(p, "--pheno", help = "UKBB phenotypes", default = "data/ukbb/nealelab/phenotypes.both_sexes.tsv.bgz")
p <- add_argument(p, "--palm", help = "PALM report for all ancestries", default = "results/palm/ancestral_paths_v3-ms-r0.05-kb250-palm_report_prs.tsv")
p <- add_argument(p, "--infectious", help = "Ony show infectious diesase phenotypes", flag = TRUE)
p <- add_argument(p, "--num-traits", help = "Number of traits to show", default = 20)
p <- add_argument(p, "--out-png", help = "Output file", default = "results/compare/ancestral_paths_v3-ms-r0.05-kb250-ukbb-upset-top.png")
p <- add_argument(p, "--out-tsv", help = "Output file", default = "results/compare/ancestral_paths_v3-ms-r0.05-kb250-ukbb-upset-top.tsv")

argv <- parse_args(p)

ukbb <- read_tsv(argv$ukbb, col_types = cols())
pheno <- read_tsv(argv$pheno, col_types = cols())
palm <- read_tsv(argv$palm, col_types = cols())

# get the list of SNPs with marginally significant p-values in CLUES in any ancestry
selected_snps <- palm %>%
    filter(significant == TRUE) %>%
    # UKBB uses `{chr}:{pos}:{ref}:{alt}` as the variant ID
    mutate(variant = paste(chrom, pos, ref, alt, sep = ":")) %>%
    pull(variant) %>%
    unique()

# only retain UKBB associations for SNPs used in the PALM analysis
ukbb_selected <- ukbb %>%
    # join the phenotype description to the code
    inner_join(
        pheno %>% select(phenotype, description),
        by = "phenotype"
    ) %>%
    # only keep selected SNPs
    filter(variant %in% selected_snps) %>%
    # remove long and unnecessary prefixes
    mutate(description = str_replace(description, "Diagnoses - main ICD10: ", "")) %>%
    mutate(description = str_replace(description, "Non-cancer illness code, self-reported: ", "")) %>%
    # add the phenotype code as a suffix
    mutate(description = paste0(str_replace(description, phenotype, ""), " [", phenotype, "]")) %>%
    # capitalize first word and strip whitespace
    mutate(description = str_squish(capitalize(description)))


if (argv$infectious) {
    # infectious disease category
    infect_keywords <- c("infection", "infectious", "virus", "viral", "bacterial", "bacteria")

    ukbb_selected <- ukbb_selected %>%
        # only keep infectious disease phenotypes
        filter(grepl(paste(infect_keywords, collapse = "|"), description, ignore.case = TRUE))
} else {
    # only retain the top overlapping phenotypes
    top_traits <- ukbb_selected %>%
        # for some measures, UKBB has both a `raw` and an `irnt` (inverse rank-normal transformed) phenotype
        filter(!grepl("_raw$", phenotype)) %>%
        # drop Celiac duplicates
        filter(!grepl("K11_COELIAC", phenotype)) %>%
        group_by(phenotype) %>%
        tally() %>%
        slice_max(order_by = n, n = as.numeric(argv$num_traits), with_ties = FALSE) %>%
        pull(phenotype)

    ukbb_selected <- ukbb_selected %>%
        filter(phenotype %in% top_traits)
}

# save the data
write_tsv(ukbb_selected, argv$out_tsv)

# append the count of SNPs to the phenotype name
ukbb_selected <- ukbb_selected %>%
    group_by(phenotype) %>%
    add_count() %>%
    mutate(phenotype = paste0(description, " (n=", n, ")")) %>%
    ungroup()

total_selected <- length(selected_snps)
total_plotted <- length(unique(ukbb_selected$variant))

# determine the max height of the bar chart
bar_height <- ukbb_selected %>%
    group_by(variant) %>%
    summarise(traits = paste0(sort(phenotype), collapse = " ")) %>%
    group_by(traits) %>%
    tally() %>%
    pull(n) %>%
    max()

# determine the size of the upset plot
upset_height <- ukbb_selected %>%
    pull(phenotype) %>%
    unique() %>%
    length()

upset_width <- ukbb_selected %>%
    group_by(variant) %>%
    summarise(traits = paste0(sort(phenotype), collapse = " ")) %>%
    pull(traits) %>%
    unique() %>%
    length()

data <- ukbb_selected %>%
    # reformat the data for plotting
    group_by(variant) %>%
    summarize(phenotypes = list(phenotype))

plot_title <- paste0(traits[[argv$trait]], " selected SNPs")

if (argv$infectious) {
    x_label <- paste0("UKBB infectious disease phenotypes and symptoms (", total_plotted, " of ", total_selected, " selected SNPs)")
} else {
    x_label <- paste0("UKBB top ", argv$num_traits, " overlapping phenotypes (", total_plotted, " of ", total_selected, " selected SNPs)")
}

plt <- data %>%
    ggplot(aes(x = phenotypes)) +
    geom_bar() +
    # stat_count(geom = "text", size = 3.5, aes(label = ..count..), position = "stack", nudge_y = 0.15) +
    # geom_text(aes(label = after_stat(paste0(round(count / total_selected * 100, 1), "%")), group = 1), stat = "count", size = 3.5, nudge_y = 0.15) +
    scale_x_upset() +
    scale_y_continuous(breaks = seq(0, bar_height, ifelse(bar_height > 10, 5, 1))) +
    ggtitle(plot_title) +
    ylab("Count of SNPs") +
    xlab(x_label)

# save the plot
ggsave(argv$out_png, plt, width = min(upset_width / 4 + 4, 30), height = min((bar_height + upset_height) / 3, 15), limitsize = FALSE)
