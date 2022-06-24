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
p <- add_argument(p, "--out-auto", help = "Output autosomal SNPs", default = "data/targets/gwas_ms-auto.tsv")
p <- add_argument(p, "--out-mhc", help = "Output MHC SNPs", default = "data/targets/gwas_ms-mhc.tsv")
p <- add_argument(p, "--out-sig", help = "Output significant SNPs", default = "data/targets/gwas_ms-auto-mhc.tsv")
p <- add_argument(p, "--out-sug", help = "Output significant + suggestive SNPs", default = "data/targets/gwas_ms-auto-mhc-sS.tsv")
p <- add_argument(p, "--out-all", help = "Output all SNPs", default = "data/targets/gwas_ms-auto-mhc-sS-wS.tsv")

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
        P = `P (joined)`, # The p-value of the joint analysis (discovery + MS Chip + ImmunoChip) for the "SNP discovery" variant
        OR = `OR (joined)`, # The OR of the joint analysis (discovery + MS Chip + ImmunoChip) for the A1 allele of the "SNP discovery" variant
        # P = `P effect conditional`, # The p-value of the "Effect" variant in the conditional model, i.e. the one that includes all variants from previous steps.
        # OR = `OR effect conditional`, # The OR of the A1 allele of the "Effect" variant in the conditional model, i.e. the one that includes all variants from previous steps.
        type,
    )

get_rsid <- function(chr, pos, a1, a2) {
    # use the Variant Effect Predictor to fetch the rsID of the SNP (which requires the non-ref allele, which we don't know)
    json_text <- tryCatch(
        {
            readLines(paste0("https://grch37.rest.ensembl.org/vep/human/region/", chr, ":", pos, "/", a1), warn = FALSE)
        },
        warning = function(e) {
            readLines(paste0("https://grch37.rest.ensembl.org/vep/human/region/", chr, ":", pos, "/", a2), warn = FALSE)
        }
    )

    rsid <- str_extract(paste(json_text, collapse = " "), "rs[0-9]+")[1]

    return(rsid)
}

# fill in any missing rsIDs, as we need them for the rest of the pipeline
fixed <- fine %>%
    filter(!grepl("rs[0-9]+", SNP)) %>%
    # for each row, fetch the rsID by calling the Ensembl web API
    rowwise() %>%
    mutate(SNP = get_rsid(CHR, BP, effect_allele, other_allele))

# merge the data back together
fine <- bind_rows(fine %>% filter(grepl("rs[0-9]+", SNP)), fixed)

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
