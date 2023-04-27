#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(tidyverse)) # v1.3.1


ms_all <- read_tsv("results/palm/ancestral_paths_v3-ms-all-palm_report_prs.tsv", col_types = cols()) %>% mutate(trait = "ms")
ra_all <- read_tsv("results/palm/ancestral_paths_v3-ra-all-palm_report_prs.tsv", col_types = cols()) %>% mutate(trait = "ra")
cd_all <- read_tsv("results/palm/ancestral_paths_v3-cd-all-palm_report_prs.tsv", col_types = cols()) %>% mutate(trait = "cd")

traits <- bind_rows(ms_all, ra_all, cd_all)

overlap <- traits %>%
    filter(significant) %>%
    group_by(rsid) %>%
    summarise(num_traits = n_distinct(trait)) %>%
    filter(num_traits > 1) %>%
    pull(rsid)

cross <- traits %>%
    filter(rsid %in% overlap) %>%
    select(variant, chrom, pos, rsid, ancestry, snp_effect, trait) %>%
    pivot_wider(names_from = trait, values_from = snp_effect)

write_tsv(cross, "results/compare/compare-ms-rs-cd-all.tsv")

cross.n <- cross %>%
    group_by(ancestry, ms, ra, cd) %>%
    tally() %>%
    arrange(ancestry, desc(n))

write_tsv(cross.n, "results/compare/compare-ms-rs-cd-all-tally.tsv")
