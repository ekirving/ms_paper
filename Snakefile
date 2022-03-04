#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"


import pandas as pd


configfile: "config.yaml"


include: "rules/gwascat.smk"
include: "rules/variants.smk"
include: "rules/clues.smk"


ANCESTRIES = ["ALL", "ANA", "CHG", "WHG", "EHG"]


def all_clumped_annotated(trait):
    """
    Run all the SNPs in the `all_clumped_annotated` ascertainments
    """
    df = pd.read_table(f"data/targets/all_clumped_annotated_{trait}.tsv")

    files = []

    for row in df.itertuples():
        rsid = row.SNP
        variant = f"{row.CHR}:{row.BP}:{row.other_allele}:{row.effect_allele}"

        files += [f"results/clues/{rsid}/ancestral_paths_new-{variant}-{ancestry}.png" for ancestry in ANCESTRIES]

    return files


rule all:
    input:
        all_clumped_annotated("ms"),
        all_clumped_annotated("ra"),
