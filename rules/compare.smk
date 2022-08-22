#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd


"""
Rules for comparing pleiotropic between traits  
"""


rule compare_gwas_catalog:
    """
    Compare trajectories between trait associated SNPs and all other traits in the GWAS catalog
    """
    input:
        catalog="data/gwascat/gwas_catalog_significant_ontology.tsv.gz",
        palm="results/palm/{dataset}-{trait}-palm_report_prs.tsv",
    output:
        png="results/compare_gwas_catalog-{dataset}-{trait}.png",
    shell:
        "Rscript scripts/compare_gwas_catalog.R"
        " --catalog {input.catalog}"
        " --palm {input.palm}"
        " --output {output.png}"
