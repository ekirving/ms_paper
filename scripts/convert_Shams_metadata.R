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
p <- arg_parser("Convert the metadata from Shams et al., 2022 (doi: 10.1093/brain/awac092)")
p <- add_argument(p, "--gwas", help = "GWAS associations", default = "data/targets/Shams_et_al_2022_S15.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "data/targets/gwas_ms-mr.tsv")

argv <- parse_args(p)

data <- read_tsv(argv$gwas, col_types = cols()) %>%
    # fetch all the expected columns
    select(CHR, BP, SNP, effect_allele = EA, other_allele = OA, P, beta = `log(OR)`, SE)

write_tsv(data, argv$output)
