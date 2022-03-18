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
p <- arg_parser("Convert the GWAS metadata into PALM input format")
p <- add_argument(p, "--gwas", help = "GWAS associations", default = "data/targets/all_clumped_annotated_ra.tsv")
p <- add_argument(p, "--variants", help = "Variant metadata", default = "data/targets/all_clumped_annotated_ra_ancestral_paths_new_variants.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "data/targets/all_clumped_annotated_ra_ancestral_paths_new_palm.tsv")

argv <- parse_args(p)

gwas <- read_tsv(argv$gwas, col_types = cols())
vars <- read_tsv(argv$variants, col_names = c("CHR", "BP", "REF", "ALT", "ANC"), col_types = cols())

data <- gwas %>%
    # join the variant metadata
    inner_join(vars, by = c("CHR", "BP")) %>%
    mutate(
        # our SNPs are already LP-pruned
        ld_block = row_number(),

        # compose the variant ID
        variant = paste(CHR, BP, REF, ALT, sep = ":"),

        # default to ALT if the ancestral allele is unknown
        derived_allele = ifelse(ANC == ALT, REF, ALT),

        # PALM assumes that betas measure the effect of the ALT allele
        beta = ifelse(effect_allele == ALT, beta, -beta)
    ) %>%
    rename(pvalue = P, rsid = SNP) %>%
    select(ld_block, variant, rsid, derived_allele, beta, se, pvalue)

write_tsv(data, argv$output)
