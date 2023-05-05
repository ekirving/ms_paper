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
quiet(library(directlabels)) # v2021.1.13
quiet(library(zoo)) # v1.8_10

# load the helper functions
source("scripts/clues_utils.R")

# DRB1*15:01 tag SNPs
description <- "DRB1*15:01"
report_tsv <- "results/palm/ancestral_paths_v3-ra-all-palm_report_prs.tsv"
output_png <- "DRB1-15-01.png"
rsid_list <- c("rs3129934", "rs3135391", "rs3135388", "rs3129889")

# # DRB1*04:01 tag SNPs
# description <- "DRB1*04:01"
# report_tsv <- "results/palm/ancestral_paths_v3-ra-all-palm_report_prs.tsv"
# output_png <- "DRB1-04:01.png"
# rsid_list <- c("rs660895", "rs6910071", "rs3817964")

# # TB SNPs
# description <- "Tuberculosis (TB)"
# report_tsv <- "results/palm/ancestral_paths_v3-ms-all-palm_report_prs.tsv"
# output_png <- "TB-SNPs.png"
# rsid_list <- c("rs422951", "rs3135363", "rs3806156", "rs2395162", "rs9268530", "rs2596472", "rs2844494")

# TB Meta-GWAS SNPs
# description <- "Tuberculosis (TB) - MetaGWAS SNPs"
# report_tsv <- "results/palm/ancestral_paths_v3-ra-all-palm_report_prs.tsv"
# output_png <- "TB-metaGWAS-SNPs.png"
# rsid_list <- c("rs2858331", "rs2857106", "rs2856993")
# rsid_list <- c("rs10245298", "rs10497744", "rs10515787", "rs10738171", "rs10956514", "rs11031728", "rs11031731", "rs11041435", "rs11041441", "rs11041442", "rs11041443", "rs11041444", "rs11041451", "rs111299137", "rs111718686", "rs111930050", "rs112848051", "rs113399573", "rs113551974", "rs1136744", "rs114204347", "rs114472823", "rs114538034", "rs115752743", "rs12362545", "rs12437118", "rs141823983", "rs144861045", "rs146049519", "rs150158763", "rs17061034", "rs17205212", "rs17394081", "rs17590261", "rs181393949", "rs188872", "rs191177369", "rs191199007", "rs191695775", "rs200554917", "rs2008560", "rs200953013", "rs2057178", "rs2269497", "rs2273061", "rs2307058", "rs2326344", "rs2621322", "rs2647071", "rs28383206", "rs28383310", "rs28383311", "rs28383323", "rs2856993", "rs2857106", "rs28578990", "rs2858331", "rs28680981", "rs358793", "rs369917557", "rs373685708", "rs376144936", "rs41541115", "rs41553512", "rs4240897", "rs4331426", "rs4348560", "rs4461087", "rs4576509", "rs4696852", "rs4733781", "rs57294695", "rs57640649", "rs58809711", "rs59979761", "rs6053639", "rs6071980", "rs6114027", "rs61321957", "rs61890220", "rs6588110", "rs6786408", "rs6913309", "rs6985962", "rs7029867", "rs71510656", "rs73226617", "rs73401332", "rs73401358", "rs73403123", "rs73409538", "rs73423701", "rs73642961", "rs76827747", "rs77924639", "rs79165841", "rs79298537", "rs79456496", "rs8006139", "rs80154958", "rs80234155", "rs80282780", "rs916943", "rs9365798", "rs9897835")

# load the palm report, with all the results
palm <- read_tsv(report_tsv, col_types = cols())

# get the models to plot
snps <- palm %>%
    filter(rsid %in% rsid_list) %>%
    # compose the model prefixes from the PALM metadata
    mutate(prefix = paste0("results/clues/", rsid, "/ancestral_paths_v3-", chrom, ":", pos, ":", ancestral_allele, ":", derived_allele, "-", ancestry))

models <- list()
for (i in 1:nrow(snps)) {
    # load the CLUES trajectory
    model <- clues_trajectory(snps[i, ]$rsid, snps[i, ]$ancestry, snps[i, ]$prefix, threshold = 0.08)
    models <- append(models, list(model))
}

# load all the trajectories
traj <- bind_rows(models) %>%
    # and join all the metadata
    inner_join(snps, by = c("rsid", "ancestry"), suffix = c("", ".focal"))

# flip `rs3135391`
traj$flip <- traj$rsid %in% c("rs3135391")

traj <- traj %>%
    mutate(
        # use the `flip` flag to polarize the trajectories
        freq = ifelse(flip, 1.0 - freq, freq),
        # add the focal allele to the SNP label
        snp_label = paste0(rsid, ":", ifelse(flip, ancestral_allele, derived_allele)),
        # add the description
        description = description
    )

# display the ancestries in custom sorted order
traj$ancestry <- factor(traj$ancestry, levels = c("ALL", "WHG", "EHG", "CHG", "ANA"))

gen_time <- 28
max_age <- 13665

# constrain the extent of the plotting
xmin <- round(-max_age / gen_time)
xmax <- max(traj$epoch)
xbreaks <- -seq(-xmax, -xmin, round(4000 / gen_time))
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

plt <- ggplot(traj) +

    # shade the ancestry epoch
    geom_rect(data = ancestry_epochs, aes(xmin = start, xmax = finish, ymin = 0, ymax = Inf), alpha = 0.5, fill = "#F4F4F4") +

    # plot the maximum posterior trajectory
    geom_line(aes(x = epoch, y = freq, color = snp_label), cex = 1, na.rm = TRUE) +

    # display as a grid
    facet_grid(~ancestry, labeller = labeller(description = label_wrap_gen())) +

    # set the SNP colors
    scale_color_manual(values = snp_colors) +

    # # print the labels
    # geom_dl(aes(x = epoch, y = freq, label = snp_label, color = snp_label, alpha = as.numeric(significant)),
    #         method = list(dl.trans(x = x + 0.1), "last.qp"), na.rm = TRUE) +

    # plot non-significant trajectories as transparent
    scale_alpha(range = c(0.3, 1), guide = "none") +
    labs(color = description) +

    # set the axis breaks
    scale_y_continuous(limits = c(0, .45), breaks = seq(0, 1, .2), position = "left") + # , expand = expansion(add = c(0.05, 0.03))) +
    scale_x_continuous(limits = c(xmin, xmax), breaks = xbreaks, labels = xlabels) + # expand = expansion(add = c(0, 300))) +
    ylab("DAF") +
    xlab("kyr BP") +

    # basic styling
    theme_minimal() +
    theme(
        # legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white", color = "white"),
        panel.spacing = unit(1, "lines"),
        text = element_text(size = 10)
    )

# save the plot
ggsave(output_png, width = 16, height = 4 * .9)

# save figure 5b for the main text
if (description == "DRB1*15:01") {
    ggsave(filename = "figure/fig_5b.png", plt, width = 806, height = 190, units = "px", dpi = 90)
}
