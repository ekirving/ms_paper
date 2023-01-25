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
p <- add_argument(p, "--gwas1", help = "The primary GWAS to be filtered", default = "data/ukbb/nealelab/gwas/20002_1456.gwas.imputed_v3.both_sexes.tsv.bgz")
p <- add_argument(p, "--gwas2", help = "The secondary GWAS", default = "data/targets/gwas_ms-full_ancestral_paths_new_palm.tsv")
p <- add_argument(p, "--vars", help = "Variant details", default = "data/ukbb/nealelab/variants.tsv.bgz")
p <- add_argument(p, "--output", help = "Output file", default = "data/ukbb/nealelab/clump/20002_1456.gwas.imputed_v3.both_sexes.ms-full.tsv")

argv <- parse_args(p)

gwas1 <- read_tsv(argv$gwas1, col_types = cols(), na = c("", "NA", "-"))
gwas2 <- read_tsv(argv$gwas2, col_types = cols())
vars <- read_tsv(argv$vars, col_types = cols()) %>% select(variant, chr, pos, ref, alt, rsid)

intersect <- vars %>%
    inner_join(gwas1, by = "variant") %>%
    inner_join(gwas2 %>% select(rsid), by = "rsid")

write_tsv(intersect, argv$output)
