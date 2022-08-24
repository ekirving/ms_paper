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
p <- arg_parser("Filter the UKBB GWAS file to retain only the SNPs with a significant GWAS assocation in the focal trait")
p <- add_argument(p, "--pheno", help = "UKBB phenotype", default = "C21")
p <- add_argument(p, "--ukbb", help = "UKBB associations", default = "data/ukbb/nealelab/gwas/C21.gwas.imputed_v3.both_sexes.significant.tsv.bgz")
p <- add_argument(p, "--palm", help = "PALM report for all ancestries", default = "results/palm/ancestral_paths_new-ms-all-palm_report_prs.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "results/palm/ukbb/ancestral_paths_new-ms-all.C21.gwas.imputed_v3.both_sexes.significant.tsv")

argv <- parse_args(p)

ukbb <- read_tsv(argv$ukbb, col_types = cols())
palm <- read_tsv(argv$palm, col_types = cols())

# get the list of SNPs used in the PALM analysis
snps <- palm %>%
    # UKBB uses `{chr}:{pos}:{ref}:{alt}` as the variant ID
    mutate(variant = paste(chrom, pos, ref, alt, sep = ":")) %>%
    pull(variant) %>%
    unique()

# only retain UKBB associations for SNPs used in the PALM analysis
ukbb <- ukbb %>%
    filter(variant %in% snps) %>%
    mutate(phenotype = argv$pheno, .before = 1) %>%
    # drop this column if it exists (so all sheets have the same columns)
    select_if(!names(.) == "expected_case_minor_AC")

write_tsv(ukbb, argv$output)
