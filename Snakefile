#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"


import pandas as pd


configfile: "config.yaml"


wildcard_constraints:
    ancestry="[A-Z]{3}",
    dataset="ancestral_paths_new",
    trait="\w+-[^_]+",


include: "rules/dbsnp.smk"
include: "rules/proxy-snps.smk"
include: "rules/gwascat.smk"
include: "rules/ukbb.smk"
include: "rules/finngen.smk"
include: "rules/variants.smk"
include: "rules/clues.smk"
include: "rules/palm.smk"
include: "rules/compare.smk"


ANCESTRIES = ["ALL", "ANA", "CHG", "WHG", "EHG"]
TRAITS = [
    "ms-r0.05-kb250",
    "ra-r0.05-kb250",
    "cd-r0.05-kb250",
]


def all_clues_plots(_):
    """
    Run CLUES for all the SNPs in all the traits
    """
    files = []

    trait = config.get("trait", False)
    traits = [trait] if trait else TRAITS

    for trait in traits:
        # noinspection PyUnresolvedReferences
        meta_tsv = checkpoints.palm_metadata_single_trait.get(dataset="ancestral_paths_new", trait=trait).output.tsv

        snp = pd.read_table(meta_tsv)

        for row in snp.itertuples():
            rsid = row.rsid
            variant = f"{row.chrom}:{row.pos}:{row.ancestral_allele}:{row.derived_allele}"

            files += [f"results/clues/{rsid}/ancestral_paths_new-{variant}-{ancestry}.png" for ancestry in ANCESTRIES]

    return files


rule all:
    input:
        # plot all the CLUES models
        all_clues_plots,
        # run PALM and plot all the traits
        expand(
            [
                "results/palm/ancestral_paths_new-{ancestry}-{trait}-palm_lines-pval.png",
                "results/palm/ancestral_paths_new-{ancestry}-{trait}-palm_lines-prs.png",
                "results/palm/ancestral_paths_new-{trait}-delta_prs.png",
                "results/palm/ancestral_paths_new-{trait}-palm_report_prs.tsv",
                "results/palm/ancestral_paths_new-{trait}-scatter.png",
            ],
            ancestry=ANCESTRIES,
            trait=config.get("trait", TRAITS),
        ),
        # make a report for each trait using the "all genome-wide SNPs" ascertainments
        expand("results/palm/ancestral_paths_new-{trait}-palm_report_prs.tsv", trait=["ms-all", "ra-all", "cd-all"]),
        # plot the UKBB and FinnGen comparisons
        expand(
            "results/compare/ancestral_paths_new-{trait}-{biobank}-marginal-001.png",
            trait=TRAITS,
            biobank=["ukbb", "finngen"],
        ),
        # make a PALM multi-trait report for all overlapping traits in UKBB and FinnGEN
        expand("results/palm/ancestral_paths_new-{trait}-palm_report_multi_trait.tsv", trait=["ms-r0.05-kb250"]),
