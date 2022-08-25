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
quiet(library(directlabels)) # v2021.1.13
quiet(library(R.utils))

# load the helper functions
source("scripts/clues_utils.R")

# get the command line arguments
p <- arg_parser("Compare associations between trait associated SNPs and all other traits in the UKBB")
p <- add_argument(p, "--ukbb", help = "UKBB assocations", default = "results/compare/ukbb/ancestral_paths_new-ms-r0.05-kb250.both_sexes.significant.tsv.gz")
p <- add_argument(p, "--pheno", help = "UKBB phenotypes", default = "data/ukbb/nealelab/phenotypes.both_sexes.tsv.bgz")
p <- add_argument(p, "--palm", help = "PALM report for all ancestries", default = "results/palm/ancestral_paths_new-ms-r0.05-kb250-palm_report_prs.tsv")
p <- add_argument(p, "--polarize", help = "How should we polarize the trajectories", default = "focal")
p <- add_argument(p, "--output", help = "Output file", default = "results/compare/ancestral_paths_new-ms-r0.05-kb250-ukbb-%03d.png")

argv <- parse_args(p)

ukbb <- read_tsv(argv$ukbb, col_types = cols())
pheno <- read_tsv(argv$pheno, col_types = cols())
palm <- read_tsv(argv$palm, col_types = cols())

# get the list of all genome-wide significant SNPs from the UKBB that intersect this trait
gwas_snps <- ukbb %>%
    pull(variant) %>%
    unique()

# join the phenotype description to the code
ukbb <- ukbb %>%
    inner_join(
        pheno %>% select(phenotype, description),
        by = "phenotype"
    ) %>%
    # remove long and unnecessary prefixes
    mutate(description = str_replace(description, "Diagnoses - main ICD10: ", "")) %>%
    mutate(description = str_replace(description, "Non-cancer illness code, self-reported: ", "")) %>%
    # add the phenotype code as a suffix
    mutate(description = paste0(str_replace(description, phenotype, ""), " (", phenotype, ")")) %>%
    # capitalize first word and strip whitespace
    mutate(description = str_squish(capitalize(description))) %>%
    # for some measures, UKBB has both a `raw` and an `irnt` (inverse rank-normal transformed) phenotype
    filter(!grepl("_raw$", phenotype))

# get the list of SNPs with marginally significant p-values in CLUES
selected_snps <- palm %>%
    filter(significant == TRUE) %>%
    pull(rsid) %>%
    unique()

# get the models to plot
snps <- palm %>%
    # only retain significant SNPs
    filter(rsid %in% selected_snps) %>%
    # UKBB uses `{chr}:{pos}:{ref}:{alt}` as the variant ID
    mutate(variant = paste(chrom, pos, ref, alt, sep = ":")) %>%
    # filter on SNPs with a UKBB association
    filter(variant %in% gwas_snps) %>%
    # compose the model prefixes from the PALM metadata
    mutate(prefix = paste0("results/clues/", rsid, "/ancestral_paths_new-", chrom, ":", pos, ":", ancestral_allele, ":", derived_allele, "-", ancestry))

# count the number of SNPs with significant CLUES p-values that intersect with each UKBB phenotype
snp_count <- ukbb %>%
    filter(variant %in% snps$variant) %>%
    group_by(phenotype, description) %>%
    tally(name = "num_snps") %>%
    arrange(desc(num_snps))

# drop phenotypes with less than the minimum number of intersecting SNPs
pheno_list <- snp_count %>%
    filter(num_snps > 1) %>%
    pull(phenotype)
ukbb <- ukbb %>% filter(phenotype %in% pheno_list)

models <- list()
for (i in 1:nrow(snps)) {
    # load the CLUES trajectory
    model <- clues_trajectory(snps[i, ]$rsid, snps[i, ]$ancestry, snps[i, ]$prefix)

    # handle the polarization
    if (argv$polarize == "focal" && snps[i, ]$flip) {
        # use the `flip` flag in the PALM metadata to polarize by risk in the focal trait
        model <- model %>% mutate(freq = 1.0 - freq, snp_label = paste0(rsid, ":", snps[i, ]$focal_allele))
    } else {
        # otherwise we're using ancestral state polarization
        model <- model %>% mutate(snp_label = paste0(rsid, ":", snps[i, ]$derived_allele))
    }

    models <- append(models, list(model))
}

# load all the trajectories
traj <- bind_rows(models) %>%
    inner_join(snps, by = c("rsid", "ancestry"), suffix = c("", ".focal")) %>%
    inner_join(ukbb, by = "variant", suffix = c(".focal", ".marginal"))


# display the phenotypes and ancestries in custom sorted order
traj$description <- factor(traj$description, levels = snp_count$description)
traj$ancestry <- factor(traj$ancestry, levels = c("ALL", "WHG", "EHG", "CHG", "ANA"))

gen_time <- 28
max_age <- 13665

# constrain the extent of the plotting
xmin <- round(-max_age / gen_time)
xmax <- max(traj$epoch)
xbreaks <- -seq(-xmax, -xmin, round(2000 / gen_time))
xlabels <- round(xbreaks * gen_time / 1000)

# how many rows do we need to print
num_rows <- traj %>%
    pull(description) %>%
    unique() %>%
    length()
per_page <- 8
num_pages <- ceiling(num_rows / per_page)

# how many columns
num_cols <- 5

for (page in 1:num_pages) {
    # subset by UKBB phenotype, so we can paginate the output (or else the plot is too big to create)
    traj_subset <- traj %>% filter(phenotype %in% pheno_list[((page - 1) * per_page):(page * per_page)])

    plt <- traj_subset %>%
        # plot the trajectories
        ggplot(aes(x = epoch, y = freq)) +

        # plot the maximum posterior trajectory
        geom_line(aes(color = snp_label, alpha = as.numeric(significant)), cex = 1, na.rm = TRUE) +

        # display as a grid
        facet_grid(description ~ ancestry, labeller = labeller(description = label_wrap_gen())) +

        # print the labels
        geom_dl(aes(label = snp_label, color = snp_label, alpha = as.numeric(significant)), method = list(dl.trans(x = x + 0.1), "last.qp", cex = 0.8), na.rm = TRUE) +

        # plot non-significant trajectories as transparent
        scale_alpha(range = c(0.3, 1)) +

        # set the axis breaks
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), position = "left", expand = expansion(add = c(0.03, 0.03))) +
        scale_x_continuous(limits = c(xmin, xmax), breaks = xbreaks, labels = xlabels, expand = expansion(add = c(0, 230))) +
        ylab("DAF") +
        xlab("kyr BP") +

        # basic styling
        theme_minimal() +
        theme(
            legend.position = "none",
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.background = element_rect(fill = "white", color = "white"),
            panel.spacing = unit(1, "lines")
        )

    # save the plot
    ggsave(sprintf(argv$output, page), plt, width = num_cols * 3, height = per_page * 2.3, limitsize = FALSE)
}
