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
p <- add_argument(p, "--gwas", help = "GWAS associations", default = "data/targets/gwas_ms-auto.tsv")
p <- add_argument(p, "--ld", help = "Pairwise LD for finging proxy SNPs", default = "data/targets/gwas_ms-auto_ld.tsv.gz")
p <- add_argument(p, "--sites", help = "List of callable sites in the current dataset", default = "data/sites/ancestral_paths_new_sites.tsv.gz")
p <- add_argument(p, "--min-ld", help = "Minimum LD threhold for proxy variants", default = 0.7)
p <- add_argument(p, "--output", help = "Output file", default = "data/targets/gwas_ms-auto_ancestral_paths_new_palm.tsv")

argv <- parse_args(p)

gwas <- read_tsv(argv$gwas, col_types = cols(), na = c("", "NA", "-"))
ld <- read_table(argv$ld, col_types = cols(), col_names = c("GWAS_CHR", "GWAS_BP", "GWAS_SNP", "PROXY_CHR", "PROXY_BP", "PROXY_SNP", "PHASE", "R2", "blank"), skip = 1) %>% select(-blank)
sites <- read_tsv(argv$sites, col_types = cols(), col_names = c("chrom", "pos", "id", "REF", "ALT", "ancestral_allele"), skip = 1) %>% select(-id)

# find the callable SNP in strongest LD with each GWAS SNP (this will retain the original SNP where possible)
proxy <- ld %>%
    # filter the LD SNPs by sites that are present in the current dataset
    inner_join(sites, by = c("PROXY_CHR" = "chrom", "PROXY_BP" = "pos")) %>%
    # pick the LD SNP with the highest R^2
    group_by(GWAS_SNP) %>%
    slice_max(R2, with_ties = TRUE) %>%
    # break R^2 ties by choosing the closest SNP
    mutate(distance = abs(GWAS_BP - PROXY_BP)) %>%
    group_by(GWAS_SNP) %>%
    slice_min(distance, with_ties = FALSE) %>%
    ungroup() %>%
    # enforce a minimum R^2 threshold
    filter(R2 >= argv$min_ld) %>%
    # split the phase block from plink into individual alleles
    separate(col = PHASE, into = c("PHASE_A1", "PHASE_A2"), sep = "/") %>%
    separate(col = PHASE_A1, into = c("GWAS_A1", "PROXY_A1"), sep = c(1)) %>%
    separate(col = PHASE_A2, into = c("GWAS_A2", "PROXY_A2"), sep = c(1))

data <- proxy %>%
    inner_join(gwas, by = c("GWAS_CHR" = "CHR", "GWAS_BP" = "BP", "GWAS_SNP" = "SNP")) %>%
    # rename the old effect and other alleles
    rename(GWAS_EA = effect_allele, GWAS_OA = other_allele) %>%
    # map the GWAS alleles onto the proxy alleles
    mutate(effect_allele = ifelse(GWAS_EA == GWAS_A1, PROXY_A1, PROXY_A2)) %>%
    mutate(other_allele = ifelse(GWAS_OA == GWAS_A1, PROXY_A1, PROXY_A2)) %>%
    # rename the proxy columns
    rename_with(~ str_remove(., "PROXY_")) %>%
    # sort the SNPs
    arrange(CHR, BP) %>%
    # add the extra columns
    mutate(
        # our SNPs are already LP-pruned
        ld_block = row_number(),

        # compose the variant ID
        variant = paste(CHR, BP, REF, ALT, sep = ":"),

        # default to ALT if the ancestral allele is unknown
        derived_allele = ifelse(ancestral_allele == ALT, REF, ALT),
        ancestral_allele = ifelse(derived_allele == ALT, REF, ALT),

        # PALM assumes that betas measure the effect of the ALT allele
        beta = ifelse(effect_allele == ALT, beta, -beta)
    ) %>%
    # remove any SNPs without a valid p-value or standard error for the association
    filter(P != 0 & se != 0) %>%
    rename(pval = P, rsid = SNP, chrom = CHR, pos = BP, ref = REF, alt = ALT) %>%
    select(ld_block, variant, chrom, pos, rsid, ref, alt, ancestral_allele, derived_allele, beta, se, pval)

write_tsv(data, argv$output)
