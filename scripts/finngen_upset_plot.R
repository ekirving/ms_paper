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
p <- arg_parser("Compare associations between trait associated SNPs and all other traits in FinnGen")
p <- add_argument(p, "--trait", help = "Focal trait name", default = "ms-r0.05-kb250")
p <- add_argument(p, "--finngen", help = "FinnGen associations", default = "results/compare/finngen/ancestral_paths_new-ms-r0.05-kb250.significant.tsv.gz")
p <- add_argument(p, "--pheno", help = "FinnGen phenotypes", default = "data/finngen/finngen_R8_manifest.tsv")
p <- add_argument(p, "--palm", help = "PALM report for all ancestries", default = "results/palm/ancestral_paths_new-ms-r0.05-kb250-palm_report_prs.tsv")
p <- add_argument(p, "--infectious", help = "Ony show infectious diesase phenotypes", flag = TRUE)
p <- add_argument(p, "--num-traits", help = "Number of traits to show", default = 20)
p <- add_argument(p, "--out-png", help = "Output file", default = "results/compare/ancestral_paths_new-ms-r0.05-kb250-finngen-upset-top.png")
p <- add_argument(p, "--out-tsv", help = "Output file", default = "results/compare/ancestral_paths_new-ms-r0.05-kb250-finngen-upset-top.tsv")

argv <- parse_args(p)

finngen <- read_tsv(argv$finngen, col_types = cols())
pheno <- read_tsv(argv$pheno, col_types = cols()) %>% rename(phenotype = phenocode, description = name)
palm <- read_tsv(argv$palm, col_types = cols())

# get the list of SNPs with marginally significant p-values in CLUES in any ancestry
selected_snps <- palm %>%
    filter(significant == TRUE) %>%
    pull(rsid) %>%
    unique()

finngen_selected <- finngen %>%
    # join the phenotype description to the code
    inner_join(
        pheno %>% select(phenotype, description, category),
        by = "phenotype"
    ) %>%
    # only keep selected SNPs
    filter(rsid %in% selected_snps)

if (argv$infectious) {
    # infectious disease category
    infect_category <- "I Certain infectious and parasitic diseases (AB1_)"

    # infectious disease symptoms
    infect_symptoms <- c(
        "ASTHMA_INFECTIONS",
        "C3_CERVIX_UTERI_EXALLC",
        "H7_IRIDOACUTE",
        "H7_OPTNEURITIS",
        "H8_OTHEREAR",
        "J10_CHRONTONSADEN",
        "J10_LOWCHRON",
        "J10_SINUSITIS",
        "L12_BULLOUS",
        "N14_CERVICAL_DYSPLASIA_ALL",
        "N14_CERVICAL_HSIL",
        "N14_FEMALE_GENITAL_DYSPLASIA_ALL",
        "N14_FEMALE_GENITAL_HSIL",
        "N14_FEMALE_GENITAL_LSIL",
        "N14_VULVAL_DYSPLASIA_ALL",
        "N14_VULVAL_LSIL",
        "O15_PREG_OTHER_MAT_DISORD",
        "RHEUMA_ARHTROPAT_REACTIVE"
    )

    finngen_selected <- finngen_selected %>%
        # only keep infectious disease phenotypes
        filter(category == infect_category | phenotype %in% infect_symptoms)
} else {
    # only retain the top overlapping phenotypes
    top_traits <- finngen_selected %>%
        # drop T1D duplicates
        filter(!grepl("T1D_|E4_DM1", phenotype)) %>%
        # drop "strict" phenotypes
        filter(!grepl("_strict", phenotype, ignore.case = TRUE)) %>%
        group_by(phenotype) %>%
        tally() %>%
        slice_max(order_by = n, n = as.numeric(argv$num_traits), with_ties = FALSE) %>%
        pull(phenotype)

    finngen_selected <- finngen_selected %>%
        filter(phenotype %in% top_traits)
}

# save the data
write_tsv(finngen_selected, argv$out_tsv)

# append the count of SNPs to the phenotype name
finngen_selected <- finngen_selected %>%
    group_by(phenotype) %>%
    add_count() %>%
    mutate(phenotype = paste0(phenotype, " (n=", n, ")")) %>%
    ungroup()

total_selected <- length(selected_snps)
total_plotted <- length(unique(finngen_selected$rsid))

# determine the max height of the bar chart
bar_height <- finngen_selected %>%
    group_by(rsid) %>%
    summarise(traits = paste0(sort(phenotype), collapse = " ")) %>%
    group_by(traits) %>%
    tally() %>%
    pull(n) %>%
    max()

# determine the size of the upset plot
upset_height <- finngen_selected %>%
    pull(phenotype) %>%
    unique() %>%
    length()

upset_width <- finngen_selected %>%
    group_by(rsid) %>%
    summarise(traits = paste0(sort(phenotype), collapse = " ")) %>%
    pull(traits) %>%
    unique() %>%
    length()

data <- finngen_selected %>%
    # reformat the data for plotting
    group_by(rsid) %>%
    summarize(phenotypes = list(phenotype))

plot_title <- paste0(traits[[argv$trait]], " selected SNPs")

if (argv$infectious) {
    x_label <- paste0("FinnGen infectious disease phenotypes and symptoms (", total_plotted, " of ", total_selected, " selected SNPs)")
} else {
    x_label <- paste0("FinnGen top ", argv$num_traits, " overlapping phenotypes (", total_plotted, " of ", total_selected, " selected SNPs)")
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
ggsave(argv$out_png, plt, width = min(upset_width / 3 + 4, 30), height = min((bar_height + upset_height) / 3, 15), limitsize = FALSE)
