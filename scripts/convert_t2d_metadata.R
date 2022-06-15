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
p <- arg_parser("Convert the metadata from Vujkovic et al., 2020 (doi: 10.1038/s41588-020-0637-y)")
p <- add_argument(p, "--gwas", help = "GWAS associations", default = "data/targets/Vujkovic_et_al_2020_ST6.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "data/targets/gwas_t2d.tsv")

argv <- parse_args(p)

data <- read_tsv(argv$gwas, col_types = cols()) %>%
    # sort the SNPs
    arrange(as.numeric(CHR), as.numeric(BP)) %>%
    # fetch all the expected columns (N,B. CHR and BP are already in GRCh37 coordinates)
    select(CHR, BP, SNP = rsid, effect_allele = EA, other_allele = NEA, P, beta = Beta, se = SE)

write_tsv(data, argv$output)
