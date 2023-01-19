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
p <- add_argument(p, "--trait1", help = "The first complex trait name", default = "cd")
p <- add_argument(p, "--trait2", help = "The second complex trait name", default = "ra")
p <- add_argument(p, "--gwas1-ind", help = "Independent GWAS associations for trait 1", default = "data/targets/gwas_cd-r0.05-kb250_ancestral_paths_new_palm.tsv")
p <- add_argument(p, "--gwas2-ind", help = "Independent GWAS associations for trait 2", default = "data/targets/gwas_ra-r0.05-kb250_ancestral_paths_new_palm.tsv")
p <- add_argument(p, "--gwas1-all", help = "All GWAS associations for trait 1", default = "data/targets/gwas_cd-full_ancestral_paths_new_palm.tsv")
p <- add_argument(p, "--gwas2-all", help = "All GWAS associations for trait 2", default = "data/targets/gwas_ra-full_ancestral_paths_new_palm.tsv")
p <- add_argument(p, "--ld", help = "Pairwise LD for finding proxy SNPs", default = "data/targets/gwas_cd-r0.05-kb250_ld.tsv.gz")
p <- add_argument(p, "--min-ld", help = "Minimum LD threshold for pairing variants", default = 0.3)
p <- add_argument(p, "--output", help = "Output file", default = "data/targets/gwas_cd~ra_ancestral_paths_new_palm.tsv")

argv <- parse_args(p)

# load all the GWAS sum stats
gwas1_all <- read_tsv(argv$gwas1_all, col_types = cols(), na = c("", "NA", "-")) %>% select(-ld_block)
gwas2_all <- read_tsv(argv$gwas2_all, col_types = cols(), na = c("", "NA", "-")) %>% select(-ld_block)
gwas1_ind <- read_tsv(argv$gwas1_ind, col_types = cols(), na = c("", "NA", "-"))
gwas2_ind <- read_tsv(argv$gwas2_ind, col_types = cols(), na = c("", "NA", "-"))

# load the pairwise LD
ld <- read_table(argv$ld, col_types = cols(), col_names = c("GWAS_CHR", "GWAS_BP", "GWAS_SNP", "PROXY_CHR", "PROXY_BP", "PROXY_SNP", "PHASE", "R2", "blank"), skip = 1) %>% select(-blank)

# PALM requires that all included SNPs have betas, p-values and standard errors for both traits
gwas_all <- gwas1_all %>%
    inner_join(gwas2_all,
        by = c("variant", "chrom", "pos", "rsid", "ref", "alt", "ancestral_allele", "derived_allele"),
        suffix = paste0("@", c(argv$trait1, argv$trait2))
    )

# drop any SNPs that are not called in both GWAS sets
trait1 <- gwas_all %>% filter(variant %in% gwas1_ind$variant)
trait2 <- gwas_all %>% filter(variant %in% gwas2_ind$variant)

# determine which of the cross-GWAS SNPs are in LD with each other
paired <- ld %>%
    # we only care about SNPs inside the 250 Kb wide LD window
    filter(abs(GWAS_BP - PROXY_BP) <= 250000) %>%
    # enforce a minimum R^2 threshold
    filter(R2 >= argv$min_ld) %>%
    # joint the two sets of independent markers
    inner_join(trait1, by = c("GWAS_CHR" = "chrom", "GWAS_BP" = "pos")) %>%
    inner_join(trait2, by = c("PROXY_CHR" = "chrom", "PROXY_BP" = "pos")) %>%
    # drop exact matches
    filter(GWAS_SNP != PROXY_SNP) %>%
    select(GWAS_SNP, PROXY_SNP, R2)

joined <- gwas_all %>%
    # get the combined set of SNPs
    filter(variant %in% c(gwas1_ind$variant, gwas2_ind$variant)) %>%
    # exclude the paired SNPs from the second GWAS
    filter(!(rsid %in% paired$PROXY_SNP)) %>%
    # put the SNPs in order
    arrange(chrom, pos) %>%
    # add the LD block numbers
    mutate(ld_block = row_number(), .before = "variant")

# now assign the paired SNPs the LD block code from the main list
others <- gwas_all %>%
    inner_join(paired %>% select(GWAS_SNP, PROXY_SNP), by = c("rsid" = "PROXY_SNP")) %>%
    inner_join(joined %>% select(ld_block, rsid), by = c("GWAS_SNP" = "rsid")) %>%
    select(-GWAS_SNP)

# and finally, put everything back together
final <- bind_rows(joined, others) %>% arrange(ld_block, pos)

write_tsv(final, argv$output)
