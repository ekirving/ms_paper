#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

library(tidyverse) # v1.3.1

# -------------------------------------------------
# Multiple sclerosis
# -------------------------------------------------

meta_tsv <- "data/targets/variants_for_evan_ms.tsv"
vars_tsv <- "variants_ms.list"
gwas_tsv <- "data/targets/gwas_ms.tsv"

meta <- read_tsv(meta_tsv, col_types = cols(.default = "c")) %>%
    select(-c(ref, alt, mismatch))

suppressWarnings(
    variants <- read_tsv(vars_tsv, col_types = cols(.default = "c"), col_names = c("chrom", "pos", "ref", "alt", "id")) %>%
        separate(id, into = c("rsid", "other"), sep = ";")
)

# check that all the SNPs are present and that the reported alleles match
meta_joined <- meta %>%
    left_join(variants, by = c("CHR" = "chrom", "BP" = "pos"))

meta_joined %>%
    filter(is.na(ref)) %>%
    pull(SNP)
# None!

meta_joined %>%
    filter((effect_allele != ref & effect_allele != alt) || (other_allele != ref & other_allele != alt)) %>%
    nrow()
# none!

meta_joined %>%
    # drop missing sites
    filter(!is.na(ref)) %>%
    # fix the missing other_allele
    mutate(other_allele = ifelse(effect_allele == ref, alt, ref)) %>%
    # standardise the column names
    select(CHR, BP, SNP, effect_allele, other_allele, P, beta, se) %>%
    write_tsv(gwas_tsv)

# -------------------------------------------------
# Celiac disease
# -------------------------------------------------

meta_tsv <- "data/targets/variants_for_evan_celiac.tsv"
vars_tsv <- "variants_celiac.list"
gwas_tsv <- "data/targets/gwas_celiac.tsv"

meta <- read_tsv(meta_tsv, col_types = cols(.default = "c"))

suppressWarnings(
    variants <- read_tsv(vars_tsv, col_types = cols(.default = "c"), col_names = c("chrom", "pos", "ref", "alt", "id")) %>%
        separate(id, into = c("rsid", "other"), sep = ";")
)

# check that all the SNPs are present and that the reported alleles match
meta_joined <- meta %>%
    left_join(variants, by = c("CHR" = "chrom", "BP" = "pos"))

meta_joined %>%
    filter(is.na(ref)) %>%
    pull(SNP)
# [1] "rs9275601"  "rs56044559" "rs3873444"

meta_joined %>%
    filter(effect_allele != ref & effect_allele != alt) %>%
    nrow()
# none!

meta_joined %>%
    # drop missing sites
    filter(!is.na(ref)) %>%
    # fix the missing other_allele
    mutate(other_allele = ifelse(effect_allele == ref, alt, ref)) %>%
    # standardise the column names
    select(CHR, BP, SNP, effect_allele, other_allele, P, beta, se = se_beta) %>%
    write_tsv(gwas_tsv)

# -------------------------------------------------
# IBD disease
# -------------------------------------------------

meta_tsv <- "data/targets/variants_for_evan_ibd.tsv"
vars_tsv <- "variants_ibd.list"
gwas_tsv <- "data/targets/gwas_ibd.tsv"

meta <- read_tsv(meta_tsv, col_types = cols(.default = "c"))

suppressWarnings(
    variants <- read_tsv(vars_tsv, col_types = cols(.default = "c"), col_names = c("chrom", "pos", "ref", "alt", "id")) %>%
        separate(id, into = c("rsid", "other"), sep = ";")
)

# check that all the SNPs are present and that the reported alleles match
meta_joined <- meta %>%
    left_join(variants, by = c("CHR" = "chrom", "BP" = "pos"))

meta_joined %>%
    filter(is.na(ref)) %>%
    pull(SNP)
# [1] "rs7517810"  "rs3197999"  "rs7714584"  "rs6556412"  "rs1799964"  "rs1456896"  "rs10758669" "rs3810936"  "rs10761659" "rs4409764"  "rs2076756"
# [12] "rs740495"   "rs1998598"  "rs2058660"  "rs212388"   "rs4077515"  "rs12722489" "rs1819658"  "rs102275"   "rs694739"   "rs4809330"  "rs181359"
# [23] "rs713875"

meta_joined %>%
    filter((effect_allele != ref & effect_allele != alt) || (other_allele != ref & other_allele != alt)) %>%
    pull(SNP) #
# [1] "rs2797685"  "rs3024505"  "rs1847472"  "rs1250550"  "rs151181"   "rs12720356"

meta_joined %>%
    # drop missing sites
    filter(!is.na(ref)) %>%
    # drop misspecified sites
    filter(effect_allele == ref || effect_allele == alt) %>%
    # fix the missing other_allele
    mutate(other_allele = ifelse(effect_allele == ref, alt, ref)) %>%
    # split the OR confidence intervals into their own columns
    mutate(OR_CI = str_remove_all(str_extract(`OR (95% CI)`, "\\(.+\\)"), "[()]")) %>%
    separate(OR_CI, into = c("OR_lower", "OR_upper"), sep = "–") %>%
    # convert OR to beta score and CI intervals to SE
    mutate(beta = log(as.numeric(OR)), se = (log(as.numeric(OR_upper)) - log(as.numeric(OR_lower))) / (2 * 1.96)) %>%
    # standardise the column names
    select(CHR, BP, SNP, effect_allele, other_allele, P, beta, se) %>%
    write_tsv(gwas_tsv)


# -------------------------------------------------
# Rheumatoid arthritis
# -------------------------------------------------

meta_tsv <- "data/targets/variants_for_evan_ra.tsv"
vars_tsv <- "variants_ra.list"
gwas_tsv <- "data/targets/gwas_ra.tsv"

meta <- read_tsv(meta_tsv, col_types = cols(.default = "c"))

suppressWarnings(
    variants <- read_tsv(vars_tsv, col_types = cols(.default = "c"), col_names = c("chrom", "pos", "ref", "alt", "id")) %>%
        separate(id, into = c("rsid", "other"), sep = ";")
)

# check that all the SNPs are present and that the reported alleles match
meta_joined <- meta %>%
    left_join(variants, by = c("CHR" = "chrom", "BP" = "pos"))

meta_joined %>%
    filter(is.na(ref)) %>%
    pull(SNP)
# [1] "rs9275601" "rs3873444"

meta_joined %>%
    filter((effect_allele != ref & effect_allele != alt) || (other_allele != ref & other_allele != alt)) %>%
    nrow()
# none!

meta_joined %>%
    # drop missing sites
    filter(!is.na(ref)) %>%
    # split the OR confidence intervals into their own columns
    separate(`OR_95%CIup-OR_95%CIlow`, into = c("OR_upper", "OR_lower"), sep = "-") %>%
    # convert OR to beta score and CI intervals to SE
    mutate(beta = log(as.numeric(OR)), se = (log(as.numeric(OR_upper)) - log(as.numeric(OR_lower))) / (2 * 1.96)) %>%
    # standardise the column names
    select(CHR, BP, SNP, effect_allele, other_allele, P, beta, se) %>%
    write_tsv(gwas_tsv)


# -------------------------------------------------
# Schizophrenia
# -------------------------------------------------

meta_tsv <- "data/targets/variants_for_evan_schizophrenia.tsv"
vars_tsv <- "variants_schizophrenia.list"
gwas_tsv <- "data/targets/gwas_schizophrenia.tsv"

meta <- read_tsv(meta_tsv, col_types = cols(.default = "c"))

suppressWarnings(
    variants <- read_tsv(vars_tsv, col_types = cols(.default = "c"), col_names = c("chrom", "pos", "ref", "alt", "id")) %>%
        separate(id, into = c("rsid", "other"), sep = ";")
)

# check that all the SNPs are present and that the reported alleles match
meta_joined <- meta %>%
    left_join(variants, by = c("CHR" = "chrom", "BP" = "pos"))

meta_joined %>%
    filter(is.na(ref)) %>%
    pull(SNP) # 110
# [1] "rs2802535"   "rs6673880"   "rs11121172"  "rs4915203"   "rs59519965"  "rs61786047"  "rs28366567"  "rs9803993"   "rs10127983"  "rs2717003"
# [11] "rs1658810"   "rs2914983"   "rs12991836"  "rs13016542"  "rs140001745" "rs3770752"   "rs72974269"  "rs1347692"   "rs12712510"  "rs9636429"
# [21] "rs999494"    "rs10178713"  "rs13032111"  "rs11692435"  "rs12713008"  "rs2909457"   "rs35026989"  "rs141216273" "rs2710323"   "rs13080668"
# [31] "rs75968099"  "rs6577597"   "rs308697"    "rs13107325"  "rs41533650"  "rs61405217"  "rs35734242"  "rs215483"    "rs6848123"   "rs6839635"
# [41] "rs10035564"  "rs10515678"  "rs177001"    "rs16867576"  "rs10052736"  "rs34879738"  "rs9267575"   "rs3131295"   "rs6922815"   "rs2261033"
# [51] "rs3130924"   "rs113113059" "rs1856507"   "rs12661552"  "rs1264369"   "rs2206956"   "rs2517524"   "rs6919146"   "rs4947288"   "rs6925079"
# [61] "rs62407633"  "rs9259545"   "rs13219424"  "rs58120505"  "rs6943762"   "rs6969410"   "rs2968532"   "rs10243922"  "rs1589726"   "rs4129585"
# [71] "rs73229090"  "rs10086619"  "rs10503253"  "rs7816998"   "rs1915019"   "rs11136325"  "rs12285419"  "rs1440480"   "rs7938083"   "rs6484367"
# [81] "rs12883788"  "rs17571951"  "rs12432904"  "rs12431743"  "rs34611983"  "rs17104975"  "rs2456020"   "rs11638554"  "rs2929278"   "rs4779050"
# [91] "rs56205728"  "rs3814883"   "rs9926049"   "rs8047438"   "rs154433"    "rs77463171"  "rs352935"    "rs72980087"  "rs9636107"   "rs7238071"
# [101] "rs4632195"   "rs56040937"  "rs2905432"   "rs322128"    "rs11696755"  "rs4328688"   "rs9975024"   "rs5751191"   "rs5995756"   "rs9623320"

meta_joined %>%
    filter((effect_allele != ref & effect_allele != alt) || (other_allele != ref & other_allele != alt)) %>%
    nrow()
# none!

meta_joined %>%
    # drop missing sites
    filter(!is.na(ref)) %>%
    # convert OR to beta score and CI intervals to SE
    mutate(beta = log(as.numeric(OR))) %>%
    # standardise the column names
    select(CHR, BP, SNP, effect_allele, other_allele, P, beta, se = SE) %>%
    write_tsv(gwas_tsv)
