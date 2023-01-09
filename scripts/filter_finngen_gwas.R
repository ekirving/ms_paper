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

# get the command line arguments
p <- arg_parser("Filter the FinnGen GWAS file to retain only the SNPs with a significant GWAS assocation in the focal trait")
p <- add_argument(p, "--pheno", help = "FinnGen phenotype", default = "R8_AB1_OTHER_BACTERIAL")
p <- add_argument(p, "--finngen", help = "FinnGen associations", default = "data/finngen/gwas/finngen_R8_AB1_OTHER_BACTERIAL.significant.tsv.bgz")
p <- add_argument(p, "--palm", help = "PALM report for all ancestries", default = "results/palm/ancestral_paths_new-ms-all-palm_report_prs.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "results/compare/finngen/ancestral_paths_new-ms-r0.05-kb250.AB1_OTHER_BACTERIAL.significant.tsv")

argv <- parse_args(p)

finngen <- read_tsv(argv$finngen, col_types = cols()) %>% separate_rows(rsids) %>% rename(rsid=rsids, chrom=`#chrom`)
palm <- read_tsv(argv$palm, col_types = cols())

# get the list of SNPs used in the PALM analysis
snps <- palm %>%
    pull(rsid) %>%
    unique()

# only retain FinnGen associations for SNPs used in the PALM analysis
finngen <- finngen %>%
    filter(rsid %in% snps) %>%
    mutate(phenotype = argv$pheno, .before = 1)

write_tsv(finngen, argv$output)
