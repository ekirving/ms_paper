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

# load the helper functions
source("scripts/clues_utils.R")

# get the command line arguments
p <- arg_parser("Add an extra column to the PALM report that contains the delta PRS for each SNP")
p <- add_argument(p, "--data", help = "PALM report", default = "results/palm/ancestral_paths_new-ALL-ms-palm_report.tsv")
p <- add_argument(p, "--dataset", help = "The dataset", default = "ancestral_paths_new")
p <- add_argument(p, "--ancestry", help = "The ancestry path", default = "ALL")
p <- add_argument(p, "--gen-time", help = "Generation time", default = 28)
p <- add_argument(p, "--output", help = "PALM trajectory", default = "results/palm/ancestral_paths_new-ALL-ms-palm_report_prs.tsv")

argv <- parse_args(p)

# TODO this only applies to MS!
# get the nearest epochs to the beginning of the visible selection peak
PEAK_START <- -round(6000 / argv$gen_time)
PEAK_END <- -round(2000 / argv$gen_time)

# load the PALM report
snps <- read_tsv(argv$data, col_types = cols())

# get the maximum possible PRS
prs_max <- sum(abs(snps$beta))

snps <- snps %>%
    # PALM assumes that betas measure the effect of the ALT allele, and CLUES models the frequency of the derived allele
    # here we want to polarize the trajectories by the positive effect allele
    mutate(flip = (alt == derived_allele & beta < 0) | (alt == ancestral_allele & beta > 0)) %>%
    # record the focal allele (to avoid any ambiguity)
    mutate(focal_allele = ifelse(beta < 0, ref, alt), .after = "derived_allele") %>%
    # normalise the beta scores by the maximum PRS
    mutate(scaled_beta = abs(beta) / prs_max) %>%
    # compose the model prefixes from the PALM metadata
    mutate(prefix = paste0("results/clues/", rsid, "/", argv$dataset, "-", chrom, ":", pos, ":", ancestral_allele, ":", derived_allele, "-", argv$ancestry))

models <- list()
for (i in 1:nrow(snps)) {
    # load the CLUES model, and extract the maximum likelihood path
    model <- clues_load_data(snps[i, ]$prefix) %>%
        group_by(epoch) %>%
        top_n(1, density) %>%
        ungroup() %>%
        arrange(epoch) %>%
        mutate(
            rsid = snps[i, ]$rsid,
            logp = -log10(snps[i, ]$p.value),
            significant = snps[i, ]$significant
        )

    if (snps[i, ]$flip) {
        # polarize the model (if necessary)
        model <- model %>%
            mutate(freq = 1.0 - freq)
    }

    model <- model %>%
        # apply the scaled PRS
        mutate(prs_freq = snps[i, ]$scaled_beta * freq)

    models <- append(models, list(model))
}

df_ml <- bind_rows(models)

# calculate the delta PRS for each SNP (i.e., the difference between the staring and ending frequency, weighted by the scaled effect size)
snp_order <- bind_rows(
    df_ml %>% group_by(rsid) %>% slice_min(epoch) %>% mutate(name = "prs_start"),
    df_ml %>% group_by(rsid) %>% slice_max(epoch) %>% mutate(name = "prs_end"),
    df_ml %>% filter(epoch == PEAK_START) %>% mutate(name = "peak_start"),
    df_ml %>% filter(epoch == PEAK_END) %>% mutate(name = "peak_end"),
) %>%
    select(rsid, name, logp, prs_freq) %>%
    pivot_wider(names_from = "name", values_from = "prs_freq") %>%
    mutate(delta_prs = prs_end - prs_start) %>%
    mutate(delta_prs_6k_2k = peak_end - peak_start) %>%
    select(rsid, delta_prs)

snp_order <- bind_rows(
    df_ml %>% group_by(rsid) %>% slice_min(epoch) %>% mutate(name = "prs_start"),
    df_ml %>% group_by(rsid) %>% slice_max(epoch) %>% mutate(name = "prs_end"),
) %>%
    select(rsid, name, logp, prs_freq) %>%
    pivot_wider(names_from = "name", values_from = "prs_freq") %>%
    mutate(delta_prs = prs_end - prs_start) %>%
    select(rsid, delta_prs)

# join the delta_prs column to the main sheet
snps <- snps %>%
    inner_join(snp_order, by = "rsid") %>%
    select(-prefix)

write_tsv(snps, argv$output)
