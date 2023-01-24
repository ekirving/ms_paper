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
p <- add_argument(p, "--gwas", help = "GWAS associations", default = "data/targets/gwas_ms-r0.05-kb250.tsv")
p <- add_argument(p, "--ld", help = "Pairwise LD for finding proxy SNPs", default = "data/targets/gwas_ms-r0.05-kb250_ld.tsv.gz")
p <- add_argument(p, "--sites", help = "List of callable sites in the current dataset", default = "data/sites/ancestral_paths_new_sites.tsv.gz")
p <- add_argument(p, "--proxy", help = "Should we replace missing GWAS SNPs with proxies", flag = TRUE)
p <- add_argument(p, "--min-ld", help = "Minimum LD threshold for proxy variants", default = 0.7)
p <- add_argument(p, "--output", help = "Output file", default = "data/targets/gwas_ms-r0.05-kb250_ancestral_paths_new_palm.tsv")

argv <- parse_args(p)

gwas <- read_tsv(argv$gwas, col_types = cols(CHR = "c"), na = c("", "NA", "-"))
sites <- read_tsv(argv$sites, col_types = cols(chrom = "c"), col_names = c("chrom", "pos", "id", "REF", "ALT", "ancestral_allele"), skip = 1) %>% select(-id)

if (argv$proxy) {
    print("INFO: Checking all GWAS SNPs are callable, and if not, searching for proxy SNPs...")

    ld <- read_table(argv$ld, col_types = cols(PROXY_CHR = "c"), col_names = c("GWAS_CHR", "GWAS_BP", "GWAS_SNP", "PROXY_CHR", "PROXY_BP", "PROXY_SNP", "PHASE", "R2", "blank"), skip = 1) %>% select(-blank)

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

    # get the list of SNPs that have been replaced
    replaced <- proxy %>%
        filter(GWAS_SNP != PROXY_SNP)

    # replace the GWAS SNPs with the callable proxies
    data <- proxy %>%
        inner_join(gwas, by = c("GWAS_CHR" = "CHR", "GWAS_BP" = "BP", "GWAS_SNP" = "SNP")) %>%
        # rename the old effect and other alleles
        rename(GWAS_EA = effect_allele, GWAS_OA = other_allele) %>%
        # map the GWAS alleles onto the proxy alleles
        mutate(effect_allele = ifelse(GWAS_EA == GWAS_A1, PROXY_A1, PROXY_A2)) %>%
        mutate(other_allele = ifelse(GWAS_OA == GWAS_A1, PROXY_A1, PROXY_A2)) %>%
        # rename the proxy columns
        rename_with(~ str_remove(., "PROXY_"))

    if (nrow(replaced) > 0) {
        print(paste0("WARN: Replaced ", nrow(replaced), " SNPs with proxy SNPs"))
    }

    if (nrow(gwas) != nrow(proxy)) {
        print(paste0("WARN: Removed ", nrow(gwas) - nrow(proxy), " SNPs with no good proxies (R^2 >= ", argv$min_ld, ")"))
    }

    if (nrow(replaced) == 0 && nrow(gwas) == nrow(proxy)) {
        print("INFO: All GWAS SNPs are present in the current dataset")
    }
} else {
    print("INFO: Checking all GWAS SNPs are callable...")

    # filter the GWAS SNPs by sites that are present in the current dataset
    data <- gwas %>%
        inner_join(sites, by = c("CHR" = "chrom", "BP" = "pos"))

    if (nrow(gwas) != nrow(data)) {
        print(paste0("WARN: Removed ", nrow(gwas) - nrow(data), " SNPs not present in the current dataset"))
    } else {
        print("INFO: All GWAS SNPs are present in the current dataset")
    }
}

data <- data %>%
    # sort the SNPs
    arrange(CHR, BP) %>%
    # fix any cases of zero p-values, due to very large negative exponents
    mutate(P = ifelse(P != 0, P, .Machine$double.xmin)) %>%
    mutate(se = ifelse(se > 0, se, abs(beta / qnorm(P / 2)))) %>%
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
    rename(pval = P, rsid = SNP, chrom = CHR, pos = BP, ref = REF, alt = ALT) %>%
    select(ld_block, variant, chrom, pos, rsid, ref, alt, ancestral_allele, derived_allele, beta, se, pval)

write_tsv(data, argv$output)
