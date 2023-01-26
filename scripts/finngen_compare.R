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
quiet(library(R.utils)) # v2.12.0

# load the helper functions
source("scripts/clues_utils.R")

# get the command line arguments
p <- arg_parser("Compare associations between trait associated SNPs and all other traits in FinnGen")
p <- add_argument(p, "--finngen", help = "FinnGen associations", default = "results/compare/finngen/ancestral_paths_new-ms-r0.05-kb250.significant.tsv.gz")
p <- add_argument(p, "--pheno", help = "FinnGen phenotypes", default = "data/finngen/finngen_R8_manifest.tsv")
p <- add_argument(p, "--palm", help = "PALM report for all ancestries", default = "results/palm/ancestral_paths_new-ms-r0.05-kb250-palm_report_prs.tsv")
p <- add_argument(p, "--polarize", help = "How should we polarize the trajectories", default = "marginal")
p <- add_argument(p, "--out-png", help = "Output file", default = "results/compare/ancestral_paths_new-ms-r0.05-kb250-finngen-marginal-%03d.png")
p <- add_argument(p, "--out-tsv", help = "Output file", default = "results/compare/ancestral_paths_new-ms-r0.05-kb250-finngen-marginal.tsv")

argv <- parse_args(p)

finngen <- read_tsv(argv$finngen, col_types = cols())
pheno <- read_tsv(argv$pheno, col_types = cols()) %>% rename(phenotype = phenocode, description = name)
palm <- read_tsv(argv$palm, col_types = cols())

# get the list of all genome-wide significant SNPs from the FinnGen that intersect this trait
gwas_snps <- finngen %>%
    pull(rsid) %>%
    unique()

# get the list of SNPs with marginally significant p-values in CLUES
selected_snps <- palm %>%
    filter(significant == TRUE) %>%
    pull(rsid) %>%
    unique()

# get the models to plot
snps <- palm %>%
    # only retain significant SNPs
    filter(rsid %in% selected_snps) %>%
    # filter on SNPs with a FinnGen association
    filter(rsid %in% gwas_snps) %>%
    # compose the model prefixes from the PALM metadata
    mutate(prefix = paste0("results/clues/", rsid, "/ancestral_paths_new-", chrom, ":", pos, ":", ancestral_allele, ":", derived_allele, "-", ancestry))

# count the number of SNPs with significant CLUES p-values that intersect with each FinnGen phenotype
snp_count <- finngen %>%
    filter(rsid %in% snps$rsid) %>%
    group_by(phenotype) %>%
    tally(name = "num_snps") %>%
    mutate(frac_snps = num_snps / length(selected_snps)) %>%
    # join the phenotype description to the code
    inner_join(
        pheno %>% select(phenotype, description),
        by = "phenotype"
    ) %>%
    # add the phenotype code as a suffix
    mutate(description = paste0(str_replace(description, phenotype, ""), " [", phenotype, "]")) %>%
    # capitalize first word and strip whitespace
    mutate(description = str_squish(capitalize(description))) %>%
    # add the SNP count to the phenotype description
    mutate(description = paste0(description, " (n=", num_snps, ")")) %>%
    # sort by the count
    arrange(desc(num_snps))

write_tsv(snp_count, argv$out_tsv)

# get the list of phenotypes with more than 1 intersecting SNPs
pheno_list <- snp_count %>%
    filter(num_snps > 1) %>%
    pull(phenotype)

finngen <- finngen %>%
    inner_join(snp_count, by = "phenotype") %>%
    # drop phenotypes with less than the minimum number of intersecting SNPs
    filter(phenotype %in% pheno_list)

models <- list()
for (i in 1:nrow(snps)) {
    # load the CLUES trajectory
    model <- clues_trajectory(snps[i, ]$rsid, snps[i, ]$ancestry, snps[i, ]$prefix, threshold = 0.08)
    models <- append(models, list(model))
}

# load all the trajectories
traj <- bind_rows(models) %>%
    # and join all the metadata
    inner_join(snps, by = c("rsid", "ancestry"), suffix = c("", ".focal")) %>%
    inner_join(finngen, by = "rsid", suffix = c(".focal", ".marginal"))

# both PALM and FinnGen report betas for the ALT allele, and CLUES models the frequency of the derived allele
if (argv$polarize == "focal") {
    # polarize by the positive effect allele in the focal trait (e.g., MS or RA)
    traj <- traj %>% mutate(flip = (alt.focal == derived_allele & beta.focal < 0) | (alt.focal == ancestral_allele & beta.focal > 0))
} else if (argv$polarize == "marginal") {
    # polarize by the positive effect allele in the FinnGen marginal trait
    traj <- traj %>% mutate(flip = (alt.focal == derived_allele & beta.marginal < 0) | (alt.focal == ancestral_allele & beta.marginal > 0))
} else {
    # the default polarization from CLUES is by ancestral/derived state
    traj <- traj %>% mutate(flip = FALSE)
}

traj <- traj %>%
    mutate(
        # use the `flip` flag to polarize the trajectories
        freq = ifelse(flip, 1.0 - freq, freq),
        # add the focal allele to the SNP label
        snp_label = paste0(rsid, ":", ifelse(flip, ancestral_allele, derived_allele))
    )

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

# the approximate time (in generations) during which these ancestries existed as discrete populations
ancestry_epochs <- tibble(
    ancestry = c("ALL", "WHG", "EHG", "CHG", "ANA"),
    start = c(-150, -500, -500, -800, -800),
    finish = c(0, -250, -200, -200, -250)
) %>%
    # don't exceed the range of the plot
    rowwise() %>%
    mutate(start = max(start, xmin), finish = min(finish, xmax))

# how many rows do we need to print
num_rows <- traj %>%
    pull(description) %>%
    unique() %>%
    length()
per_page <- 5
num_pages <- ceiling(num_rows / per_page)

# how many columns
num_cols <- 5

for (page in 1:num_pages) {
    # subset by FinnGen phenotype, so we can paginate the output (or else the plot is too big to create)
    pheno_subset <- pheno_list[((page - 1) * per_page):(page * per_page)]
    traj_subset <- traj %>% filter(phenotype %in% pheno_subset)
    num_pheno <- length(pheno_subset[!is.na(pheno_subset)])

    plt <- ggplot(traj_subset) +

        # shade the ancestry epoch
        geom_rect(data = ancestry_epochs, aes(xmin = start, xmax = finish, ymin = 0, ymax = 1), alpha = 0.5, fill = "#F4F4F4") +

        # plot the maximum posterior trajectory
        geom_line(aes(x = epoch, y = freq, color = snp_label, alpha = as.numeric(significant)), cex = 1, na.rm = TRUE) +

        # display as a grid
        facet_grid(description ~ ancestry, labeller = labeller(description = label_wrap_gen())) +

        # print the labels
        geom_dl(aes(x = epoch, y = freq, label = snp_label, color = snp_label, alpha = as.numeric(significant)), method = list(dl.trans(x = x + 0.1), "last.qp", cex = 0.8), na.rm = TRUE) +

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
    ggsave(sprintf(argv$out_png, page), plt, width = num_cols * 3.3, height = num_pheno * 2.3, limitsize = FALSE)
}
