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
p <- add_argument(p, "--gwas1", help = "The primary GWAS to be filtered", default = "data/finngen/gwas/finngen_R8_T1D.significant.tsv.bgz")
p <- add_argument(p, "--gwas2", help = "The secondary GWAS", default = "data/targets/gwas_ms-full_ancestral_paths_v3_palm.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "data/finngen/clump/finngen_R8_T1D.significant.ms-full.tsv")

argv <- parse_args(p)

gwas1 <- read_tsv(argv$gwas1, col_types = cols())
gwas2 <- read_tsv(argv$gwas2, col_types = cols())

intersect <- gwas1 %>%
    separate_rows(rsids) %>%
    rename(rsid = rsids, chrom = `#chrom`) %>%
    # FinnGen uses GRCh38 not GRCh37, so we need to switch the positions
    select(-c(chrom, pos)) %>%
    inner_join(gwas2 %>% select(rsid, chrom, pos), by = "rsid")

write_tsv(intersect, argv$output)
