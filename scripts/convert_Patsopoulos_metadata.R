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
p <- arg_parser("Convert the metadata for Multiple sclerosis")
p <- add_argument(p, "--auto", help = "Genome-wide significant autosomal SNPs", default = "data/targets/Patsopoulos_et_al_2019_ST7.tsv")
p <- add_argument(p, "--mhc", help = "MHC SNPs, classical model", default = "data/targets/Patsopoulos_et_al_2019_ST11.tsv")
p <- add_argument(p, "--sS", help = "Strongly suggestive SNPs", default = "data/targets/Patsopoulos_et_al_2019_ST14_sS.tsv")
p <- add_argument(p, "--wS", help = "Weakly suggestive SNPs", default = "data/targets/Patsopoulos_et_al_2019_ST14_wS.tsv")
p <- add_argument(p, "--out-auto", help = "Output autosomal SNPs", default = "data/targets/gwas_ms_auto.tsv")
p <- add_argument(p, "--out-mhc", help = "Output MHC SNPs", default = "data/targets/gwas_ms_mhc.tsv")
p <- add_argument(p, "--out-sig", help = "Output significant SNPs", default = "data/targets/gwas_ms_auto_mhc.tsv")
p <- add_argument(p, "--out-sug", help = "Output significant + suggestive SNPs", default = "data/targets/gwas_ms_auto_mhc_sS.tsv")
p <- add_argument(p, "--out-all", help = "Output all SNPs", default = "data/targets/gwas_ms_auto_mhc_sS_wS.tsv")

argv <- parse_args(p)

# load all the fine-mapped data
fine <- bind_rows(
    read_tsv(argv$auto, col_types = cols(Step = col_character())) %>% mutate(type = "auto"),
    read_tsv(argv$sS, col_types = cols(Step = col_character())) %>% mutate(type = "sS"),
    read_tsv(argv$wS, col_types = cols(Step = col_character())) %>% mutate(type = "wS"),
) %>%
    # fetch all the expected columns
    select(
        CHR = Chromosome,
        BP = `Position (hg19)`,
        SNP = `Effect`,
        effect_allele = `A1 effect`,
        other_allele = `A2 effect`,
        P = `P effect conditional`, # The p-value of the "Effect" variant in the conditional model, i.e. the one that includes all variants from previous steps.
        OR = `OR effect conditional`, # The OR of the A1 allele of the "Effect" variant in the conditional model, i.e. the one that includes all variants from previous steps.
        type,
    )

# handle the MHC SNPs from the classical model
mhc <- read_tsv(argv$mhc, col_types = cols()) %>% mutate(CHR = 6)

# TODO how best to handle the 12 HLA allele calls (e.g., `HLA-DRB1*15:01`, `AA B position 45 TK`)
# mhc %>% filter(!grepl("rs[0-9]+", Variant))

mhc_snps <- mhc %>%
    # get the associations with rsIDs
    filter(grepl("rs[0-9]+", Variant)) %>%
    mutate(type = "mhc_snps") %>%
    # SNPs with alleles in the variant column represent multiallelic SNPs
    separate(Variant, into = c("Variant", "allele"), fill = "right") %>%
    mutate(A1 = coalesce(allele, A1), A2 = ifelse(is.na(allele), A2, NA)) %>%
    # fetch all the expected columns
    select(
        CHR,
        BP = `Position (hg19)`,
        SNP = `Variant`,
        effect_allele = `A1`,
        other_allele = `A2`,
        P = `P-value`,
        OR,
        type,
    )

# merge all the associations
assoc <- bind_rows(fine, mhc_snps) %>%
    # convert OR to beta score
    mutate(beta = log(as.numeric(OR)), .keep = "unused", .before = type) %>%
    # calculate the SE from the beta and p-value
    mutate(se = abs(beta / qnorm(P / 2)), .before = type) %>%
    # sort the SNPs
    arrange(CHR, BP)


# TODO what about the non-rsID SNP names?

# save the genome-wide significant autosomal SNPs
assoc %>%
    filter(type == "auto") %>%
    write_tsv(argv$out_auto)

# save the MHC SNPs
assoc %>%
    filter(type == "mhc_snps") %>%
    write_tsv(argv$out_mhc)

# save all the genome-wide significant autosomal SNPs (auto + MHC)
assoc %>%
    filter(type %in% c("auto", "mhc_snps")) %>%
    write_tsv(argv$out_sig)

# save significant and strongly suggestive SNPs
assoc %>%
    filter(type %in% c("auto", "mhc_snps", "sS")) %>%
    write_tsv(argv$out_sug)

# save all the associations
assoc %>% write_tsv(argv$out_all)
