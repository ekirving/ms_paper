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
p <- arg_parser("Convert the metadata from Ju and Mathieson, 2021 (doi: 10.1073/pnas.2009227118)")
p <- add_argument(p, "--gwas", help = "GWAS associations", default = "data/targets/Ju_and_Mathieson_2021_S1b.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "data/targets/gwas_skinp.tsv")

argv <- parse_args(p)

data <- read_tsv(argv$gwas, col_types = cols()) %>%
    # sort the SNPs
    arrange(as.numeric(CHROM), as.numeric(POS)) %>%
    # fetch all the expected columns (N,B. CHR and BP are already in GRCh37 coordinates)
    select(CHR = CHROM, BP = POS, SNP = RSID, effect_allele = Dark, other_allele = Light, P = p, beta = Beta, se = SE)

write_tsv(data, argv$output)
