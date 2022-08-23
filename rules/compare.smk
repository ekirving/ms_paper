#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd


"""
Rules for comparing pleiotropic associations between the focal trait (e.g., MS, RA) and a database of associations (e.g., GWAS Catalog)
"""


rule compare_gwas_catalog:
    """
    Compare trajectories between trait associated SNPs and all other traits in the GWAS catalog
    """
    input:
        catalog="data/gwascat/gwas_catalog_significant_ontology.tsv.gz",
        palm="results/palm/{dataset}-{trait}-palm_report_prs.tsv",
    output:
        png="results/compare/{dataset}-{trait}-gwas_catalog.png",
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
        ukbb="data/ukbb/nealelab/gwas/{pheno}.gwas.imputed_v3.{sex}.significant.tsv.bgz",
        palm="results/palm/{dataset}-{trait}-palm_report_prs.tsv",
    output:
        tsv=temp("results/compare/ukbb/{dataset}-{trait}.{pheno}.gwas.imputed_v3.{sex}.significant.tsv"),
    shell:
        "Rscript scripts/filter_ukbb_gwas.R"
        " --pheno {wildcards.pheno}"
        " --ukbb {input.ukbb}"
        " --palm {input.palm}"
        " --output {output.tsv}"


def merge_all_filtered_phenotypes_input(wildcards):
    """
    Return a list of all the phenotype association files for the NealeLab UKBB GWAS
    """
    # noinspection PyUnresolvedReferences
    pheno_tsv = checkpoints.ukbb_nealelab_phenotypes.get(sex=wildcards.sex).output.bgz

    phenotypes = pd.read_table(pheno_tsv, compression="gzip")["phenotype"].unique()

    return expand(
        "results/compare/ukbb/{dataset}-{trait}.{pheno}.gwas.imputed_v3.{sex}.significant.tsv",
        pheno=phenotypes,
        **wildcards
    )


rule merge_all_filtered_phenotypes:
    """
    Merge all the filtered UKBB association files
    """
    input:
        merge_all_filtered_phenotypes_input,
    output:
        tsv="results/compare/ukbb/{dataset}-{trait}.{sex}.significant.tsv.gz",
    shell:
        "head -n1 {input[0]} | gzip -c > {output.tsv} && "
        "tail -n +2 -q {input} | gzip -c >> {output.tsv}"


rule compare_ukbb:
    """
    Compare trajectories between trait associated SNPs and all other traits in the UKBB
    """
    input:
        ukbb="results/palm/ukbb/{dataset}-{trait}.both_sexes.significant.tsv.gz",
        pheno="data/ukbb/nealelab/phenotypes.both_sexes.tsv.bgz",
        palm="results/palm/{dataset}-{trait}-palm_report_prs.tsv",
    output:
        png="results/compare/{dataset}-{trait}-ukbb.png",
    shell:
        "Rscript scripts/compare_ukbb.R"
        " --ukbb {input.ukbb}"
        " --pheno {input.pheno}"
        " --palm {input.palm}"
        " --output {output.png}"
