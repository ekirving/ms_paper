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
    dataset="[^-]+",
    trait="[^_]+",


include: "rules/dbsnp.smk"
include: "rules/proxy-snps.smk"
include: "rules/gwascat.smk"
include: "rules/variants.smk"
include: "rules/clues.smk"
include: "rules/palm.smk"


ANCESTRIES = ["ALL", "ANA", "CHG", "WHG", "EHG"]
TRAITS = [
    "cd-all",
    "ms-all",
    "ra-all",
]


def all_clues_plots(_):
    """
    Run CLUES for all the SNPs in all the traits
    """
    files = []

    for trait in TRAITS:
        # noinspection PyUnresolvedReferences
        meta_tsv = checkpoints.palm_metadata.get(dataset="ancestral_paths_new", trait=trait).output.tsv

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
            ],
            ancestry=ANCESTRIES,
            trait=config["trait"],
        ),
