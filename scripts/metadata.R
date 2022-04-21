#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

library(tidyverse) # v1.3.1

setwd("/Users/evan/Dropbox/Code/ms")

ms <- read_tsv("data/targets/all_clumped_annotated_ms.tsv", col_types = cols(.default = "c"))
ra <- read_tsv("data/targets/all_clumped_annotated_ra.tsv", col_types = cols(.default = "c"))

# table_cols <- colnames(ms)
#
# # convert the `ra` columns to the same format as `ms`
# ra <- ra %>%
#     rename(effect_allele = A1, other_allele = A2, OR = `OR(A1)`) %>%
#     separate(`OR_95%CIup-OR_95%CIlow`, into = c("OR_upper", "OR_lower"), sep = "-") %>%
#     mutate(beta = log(as.numeric(OR)), se = NA) %>%
#     select(all_of(table_cols))
#
# write_tsv(ra, "data/targets/all_clumped_annotated_ra.tsv")

# -------------------------------------------------
# multiple sclerosis
# -------------------------------------------------

# check that all the SNPs are present and that the reported alleles match
variants_ms <- read_tsv("variants_ms.list", col_types = cols(.default = "c"), col_names = c("chrom", "pos", "ref", "alt", "id")) %>%
    separate(id, into = c("rsid", "other"), sep = ";")

joined_ms <- ms %>%
    left_join(variants_ms, by = c("CHR" = "chrom", "BP" = "pos"))

joined_ms %>%
    filter(is.na(ref)) %>%
    pull(SNP)
# none!

joined_ms %>%
    filter((effect_allele != ref & effect_allele != alt) | (other_allele != ref & other_allele != alt)) %>%
    nrow() # 91
# none!

# -------------------------------------------------
# rheumatoid arthritis
# -------------------------------------------------

variants_ra <- read_tsv("variants_ra.list", col_types = cols(.default = "c"), col_names = c("chrom", "pos", "ref", "alt", "id")) %>%
    separate(id, into = c("rsid", "other"), sep = ";")

joined_ra <- ra %>%
    left_join(variants_ra, by = c("CHR" = "chrom", "BP" = "pos"))

joined_ra %>%
    filter(is.na(ref)) %>%
    View()
# "rs9275601" "rs3873444"

joined_ra %>%
    filter((effect_allele != ref & effect_allele != alt) | (other_allele != ref & other_allele != alt)) %>%
    nrow()
# none!
