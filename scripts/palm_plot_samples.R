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
quiet(library(jsonlite)) # v1.8.0

# get the command line arguments
p <- arg_parser("Plot the density of the sampling and ancestry painting over time")
p <- add_argument(p, "--tsv", help = "PALM report", default = "results/palm/ancestral_paths_new-ALL-ms-palm_report.tsv")
p <- add_argument(p, "--json", help = "PALM json file", default = "results/palm/ancestral_paths_new-ALL-ms-palm.json")
p <- add_argument(p, "--trait", help = "The complex trait name", default = "ms")
p <- add_argument(p, "--dataset", help = "The dataset", default = "ancestral_paths_new")
p <- add_argument(p, "--gen-time", help = "Generation time", default = 28)
p <- add_argument(p, "--output", help = "PALM trajectory", default = "results/palm/ancestral_paths_new-ALL-ms-palm_trajectory.png")

argv <- parse_args(p)

# load the PALM results
results <- fromJSON(argv$json)

# load the list of SNPs and their effect sizes
snps <- read_tsv(argv$tsv, col_types = cols())

snps <- snps %>%
    # compose the model prefixes from the PALM metadata
    mutate(prefix = paste0("results/clues/", rsid, "/", argv$dataset, "-", chrom, ":", pos, ":", ancestral_allele, ":", derived_allele))

snps <- head(snps, n=5)


models <- list()
for (i in 1:nrow(snps)) {
    # load the diploid pan-ancestry models
    read_table(paste0(snps[i, ]$prefix, "-ALL.ancient"), col_names = c("generations", "hom_ref", "het", "hom_alt"), col_types = cols()) %>%
        mutate(ancestry="ALL", rsid=snps[i, ]$rsid)
}

for (i in 1:nrow(snps)) {
    for (ancestry in c("ANA", "CHG", "EHG", "WHG")) {
        # load the CLUES genotype files
        read_table(paste0(snps[i, ]$prefix, "-", ancestry, ".ancient"), col_names = c("generations", "hom_ref", "hom_alt"), col_types = cols()) %>%
            mutate(ancestry=ancestry, rsid=snps[i, ]$rsid)
        
    }
}





