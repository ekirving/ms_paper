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
p <- add_argument(p, "--gwas", help = "GWAS associations", default = "data/targets/gwas_ms.tsv")
p <- add_argument(p, "--variants", help = "Variant metadata", default = "data/targets/gwas_ms_ancestral_paths_new_variants.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "data/targets/gwas_ms_ancestral_paths_new_palm.tsv")

argv <- parse_args(p)

gwas <- read_tsv(argv$gwas, col_types = cols(), na = c("", "NA", "-"))
vars <- read_tsv(argv$variants, col_names = c("CHR", "BP", "ID", "REF", "ALT", "ANC"), col_types = cols())

data <- gwas %>%
    # join the variant metadata (this drops any SNPs that are not in the callset)
    inner_join(vars, by = c("CHR", "BP")) %>%
    # fill any missing rsIDs
    separate(col = ID, into = c("rsid", "chr_pos"), sep = ";", fill = "left") %>%
    mutate(SNP = coalesce(SNP, rsid)) %>%
    select(-rsid) %>%
    # add the extra columns
    mutate(
        # our SNPs are already LP-pruned
        ld_block = row_number(),

        # compose the variant ID
        variant = paste(CHR, BP, REF, ALT, sep = ":"),

        # default to ALT if the ancestral allele is unknown
        derived_allele = ifelse(ANC == ALT, REF, ALT),
        ancestral_allele = ifelse(derived_allele == ALT, REF, ALT),

        # PALM assumes that betas measure the effect of the ALT allele
        beta = ifelse(effect_allele == ALT, beta, -beta)
    ) %>%
    # remove any SNPs without a valid p-value or standard error for the association
    filter(P != 0 & se != 0) %>%
    rename(pval = P, rsid = SNP, chrom = CHR, pos = BP, ref = REF, alt = ALT) %>%
    select(ld_block, variant, chrom, pos, rsid, ref, alt, ancestral_allele, derived_allele, beta, se, pval)

write_tsv(data, argv$output)
