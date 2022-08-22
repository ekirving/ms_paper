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


rule filter_ukbb_gwas:
    """
    Filter the UKBB GWAS file to retain only the SNPs with a significant selection test in the CLUES analysis
    """
    input:
        gwas="data/ukbb/nealelab/gwas/{pheno}.gwas.imputed_v3.{sex}.significant.tsv.bgz",
        palm="results/palm/{dataset}-{trait}-palm_report_prs.tsv",
    output:
        tsv=temp("results/palm/ukbb/{dataset}-{trait}.{pheno}.gwas.imputed_v3.{sex}.significant.tsv"),
    shell:
        "Rscript scripts/filter_ukbb_phenotypes.R"
        " --pheno {wildcards.pheno}"
        " --gwas {input.gwas}"
        " --palm {input.palm}"
        " --output {output.tsv}"
