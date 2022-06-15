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
p <- arg_parser("Convert the metadata from Robertson et al., 2021 (doi: 10.1038/s41588-021-00880-5)")
p <- add_argument(p, "--gwas", help = "GWAS associations", default = "data/targets/Robertson_et_al_2021_ST7.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "data/targets/gwas_t1d.tsv")

argv <- parse_args(p)

data <- read_tsv(argv$gwas, col_types = cols()) %>%
    # N.B. `alleleB` is the effect allele; we use the effect sizes from the Phase I meta-analysis
    select(SNP = ID, effect_allele = alleleB, other_allele = alleleA, P = pphaseI, beta = EffectphaseI, se = SEphaseI) %>%
    # for each row, fetch the GRCh37 coordinates by calling the Ensembl web API
    rowwise() %>%
    mutate(location = fromJSON(paste0("https://grch37.rest.ensembl.org/variation/human/", SNP))$mappings$location[1]) %>%
    separate(location, into = c("CHR", "BP", "END")) %>%
    # sort the SNPs
    arrange(as.numeric(CHR), as.numeric(BP)) %>%
    # fetch all the expected columns
    select(CHR, BP, SNP, effect_allele, other_allele, P, beta, se)

write_tsv(data, argv$output)
