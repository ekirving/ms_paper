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
p <- arg_parser("Filter one GWAS based on the positions in another")
p <- add_argument(p, "--gwas", help = "The GWAS to be filtered", default = "data/ukbb/nealelab/clump/20002_1456.gwas.imputed_v3.both_sexes.ms-full.tsv")
p <- add_argument(p, "--clump", help = "The secondary GWAS", default = "data/ukbb/nealelab/clump/20002_1456.gwas.imputed_v3.both_sexes.ms-full.clumped")
p <- add_argument(p, "--output1", help = "Output file", default = "data/targets/gwas_20002_1456-ukbb-r0.05-kb250.tsv")
p <- add_argument(p, "--output2", help = "Output file", default = "data/targets/gwas_20002_1456-ukbb-full.tsv")

argv <- parse_args(p)

gwas <- read_tsv(argv$gwas, col_types = cols())
clump <- read_table(argv$clump, col_types = cols())

data1 <- gwas %>%
    filter(rsid %in% clump$SNP) %>%
    mutate(OR = exp(beta)) %>%
    # fetch all the expected columns
    select(SNP = rsid, CHR = chr, BP = pos, effect_allele = alt, other_allele = ref, OR, beta, se, P = pval)

write_tsv(data1, argv$output1)


data2 <- gwas %>%
    mutate(OR = exp(beta)) %>%
    # fetch all the expected columns
    select(SNP = rsid, CHR = chr, BP = pos, effect_allele = alt, other_allele = ref, OR, beta, se, P = pval)

write_tsv(data2, argv$output2)
