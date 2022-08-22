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

# load the helper functions
source("scripts/clues_utils.R")

# get the command line arguments
p <- arg_parser("Compare associations between trait associated SNPs and all other traits in the GWAS catalog")
p <- add_argument(p, "--catalog", help = "GWAS catalog", default = "data/gwascat/gwas_catalog_significant_ontology.tsv.gz")
p <- add_argument(p, "--palm", help = "PALM report for all ancestries", default = "results/palm/ancestral_paths_new-ms-all-palm_report_prs.tsv")
# p <- add_argument(p, "--output", help = "Output file", default = "data/targets/gwas_ms.tsv")

argv <- parse_args(p)

gwas_catalog <- read_tsv(argv$catalog, col_types = cols(), guess_max = 40000)
palm <- read_tsv(argv$palm, col_types = cols())

# filter the GWAS catalog
assocations <- gwas_catalog %>%
    # select only the necessary columns
    select(snp = `STRONGEST SNP-RISK ALLELE`, or_beta = `OR or BETA`, trait = `MAPPED_TRAIT`, ontology = `ONTOLOGY_LABELS`) %>%
    # split records with multiple SNPs
    separate_rows(snp, sep = "[;x]") %>%
    # trim whitespace
    mutate(snp = str_trim(snp)) %>%
    # only retain associations with an rsID and risk allele
    filter(grepl("^rs[0-9]+ ?- ?[ACGT?]+$", snp, ignore.case = TRUE)) %>%
    # split rsID and risk allele
    separate(snp, into = c("rsid", "risk_allele")) %>%
    # standardise case
    mutate(risk_allele = toupper(risk_allele))

# get the list of SNPs in the GWAS Catalog
gwas_snps <- assocations %>%
    pull(rsid) %>%
    unique()

# get the list of SNPs with marginally significant p-values in CLUES
selected_snps <- palm %>%
    filter(significant == TRUE) %>%
    pull(rsid) %>%
    unique()

# get the models to plot
snps <- palm %>%
    # only retain significant SNPs
    filter(rsid %in% selected_snps) %>%
    # with a GWAS Catalog entry
    filter(rsid %in% gwas_snps) %>%
    # compose the model prefixes from the PALM metadata
    mutate(prefix = paste0("results/clues/", rsid, "/ancestral_paths_new-", chrom, ":", pos, ":", ancestral_allele, ":", derived_allele, "-", ancestry))

models <- list()
for (i in 1:nrow(snps)) {
    # load the CLUES trajectory
    model <- clues_trajectory(snps[i, ]$rsid, snps[i, ]$ancestry, snps[i, ]$prefix)
    models <- append(models, list(model))
}

# load all the trajectories
traj <- bind_rows(models) %>%
    # join the GWAS catalog
    inner_join(assocations, by = "rsid")

gen_time <- 28
max_age <- 13665

# constrain the extent of the plotting
xmin <- round(-max_age / gen_time)
xmax <- max(traj$epoch)
xbreaks <- -seq(-xmax, -xmin, round(2000 / gen_time))
xlabels <- round(xbreaks * gen_time / 1000)

plt <- traj %>%
    # plot the heatmap
    ggplot(aes(x = epoch, y = freq)) +
    
    # plot the maximum posterior trajectory
    geom_line(aes(color = rsid), cex = 1, na.rm = TRUE) +
    
    # diplay as a grid
    facet_grid(trait ~ ancestry) +
    
    # print the labels
    geom_dl(aes(label = rsid, color = rsid), method = list(dl.trans(x = x + 0.1), "last.qp", cex = 0.8), na.rm = TRUE) +

    # set the axis breaks
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), position = "left", expand = expansion(add = c(0.03, 0.03))) +
    scale_x_continuous(limits = c(xmin, xmax), breaks = xbreaks, labels = xlabels, expand = expansion(add = c(0, 230))) +

    ylab("DAF") +
    xlab("kyr BP") +

    # basic styling
    theme_minimal() +
    theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white", color = "white"),
        panel.spacing = unit(1, "lines")
    )

# size the plot based on the facet dimensions
num_traits <- assocations %>%
    filter(rsid %in% selected_snps) %>%
    pull(trait) %>%
    unique() %>%
    length()
num_ancestries <- palm %>%
    pull(ancestry) %>%
    unique() %>%
    length()

# save the plot
ggsave("pleiotropy.png", plt, width = num_ancestries * 3, height = num_traits * 2, limitsize = FALSE)
