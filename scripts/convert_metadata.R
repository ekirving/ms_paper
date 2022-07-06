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
p <- arg_parser("Convert the metadata for Multiple sclerosis")
p <- add_argument(p, "--gwas", help = "GWAS associations", default = "data/targets/ms_snps_final_discovery_0.7_combined.csv")
p <- add_argument(p, "--output", help = "Output file", default = "data/targets/gwas_ms.tsv")

argv <- parse_args(p)

data <- read_csv(argv$gwas, col_types = cols()) %>%
    # fetch all the expected columns
    select(
        CHR = proxy_chr,
        BP = proxy_position,
        SNP = proxy_rs,
        effect_allele = proxy_effect_allele,
        other_allele = proxy_other_allele,
        P = finemapped_p_value,
        beta = finemapped_est,
        se = finemapped_se
    ) %>%
    arrange(CHR, BP)

write_tsv(data, argv$output)
