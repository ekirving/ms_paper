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
p <- arg_parser("Convert the metadata from the GWAS Catalog format")
p <- add_argument(p, "--gwas", help = "GWAS associations", default = "data/targets/gwas-association-accessionId_GCST90013445.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "data/targets/gwas_GCST90013445.tsv")

argv <- parse_args(p)

data <- read_tsv(argv$gwas, col_types = cols()) %>%
    # rename the columns
    select(SNP = `STRONGEST SNP-RISK ALLELE`, P = `P-VALUE`, OR = `OR or BETA`, OR_CI = `95% CI (TEXT)`) %>%
    # split the rsID and effect allele
    separate(SNP, into = c("SNP", "effect_allele"), sep = "-") %>%
    # split the OR confidence intervals into their own columns
    mutate(OR_CI = str_remove_all(OR_CI, "[\\[\\]]")) %>%
    separate(OR_CI, into = c("OR_lower", "OR_upper"), sep = "-") %>%
    # convert OR to beta score and CI intervals to SE
    mutate(beta = log(as.numeric(OR)), se = (log(as.numeric(OR_upper)) - log(as.numeric(OR_lower))) / (2 * 1.96))

data <- data %>%
    # for each row, fetch the GRCh37 coordinates by calling the Ensembl web API
    rowwise() %>%
    mutate(location = fromJSON(paste0("https://grch37.rest.ensembl.org/variation/human/", SNP))$mappings$location[1]) %>%
    separate(location, into = c("CHR", "BP", "END")) %>%
    # sort the SNPs
    arrange(as.numeric(CHR), as.numeric(BP)) %>%
    # fetch all the expected columns
    select(CHR, BP, SNP, effect_allele, P, beta, se)

write_tsv(data, argv$output)
