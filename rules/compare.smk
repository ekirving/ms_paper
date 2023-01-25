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


rule gwas_catalog_compare:
    """
    Compare trajectories between trait associated SNPs and all other traits in the GWAS catalog
    """
    input:
        catalog="data/gwascat/gwas_catalog_significant_ontology.tsv.gz",
        palm="results/palm/{dataset}-{trait}-palm_report_prs.tsv",
    output:
        png="results/compare/{dataset}-{trait}-gwas_catalog.png",
    shell:
        "Rscript scripts/gwas_catalog_compare.R"
        " --catalog {input.catalog}"
        " --palm {input.palm}"
        " --output {output.png}"


rule ukbb_filter_gwas:
    """
    Filter the UKBB GWAS file to retain only the SNPs with a significant selection test in the CLUES analysis
    """
    input:
        ukbb="data/ukbb/nealelab/gwas/{pheno}.gwas.imputed_v3.{sex}.significant.tsv.bgz",
        palm="results/palm/{dataset}-{trait}-palm_report_prs.tsv",
    output:
        tsv=temp("results/compare/ukbb/{dataset}-{trait}.{pheno}.gwas.imputed_v3.{sex}.significant.tsv"),
    shell:
        "Rscript scripts/ukbb_filter_gwas.R"
        " --pheno {wildcards.pheno}"
        " --ukbb {input.ukbb}"
        " --palm {input.palm}"
        " --output {output.tsv}"


def ukbb_merge_all_filtered_phenotypes_input(wildcards):
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


rule ukbb_merge_all_filtered_phenotypes:
    """
    Merge all the filtered UKBB association files
    """
    input:
        ukbb_merge_all_filtered_phenotypes_input,
    output:
        tsv="results/compare/ukbb/{dataset}-{trait}.{sex}.significant.tsv.gz",
    shell:
        "head -n1 {input[0]} | gzip -c > {output.tsv} && "
        "tail -n +2 -q {input} | gzip -c >> {output.tsv}"


rule ukbb_compare:
    """
    Compare trajectories between trait associated SNPs and all other traits in the UKBB
    """
    input:
        ukbb="results/compare/ukbb/{dataset}-{trait}.both_sexes.significant.tsv.gz",
        pheno="data/ukbb/nealelab/phenotypes.both_sexes.tsv.bgz",
        palm="results/palm/{dataset}-{trait}-palm_report_prs.tsv",
    output:
        png="results/compare/{dataset}-{trait}-ukbb-{polarize}-001.png",
    params:
        # plotting script produces multiple PNG files (based on number of UKBB traits)
        png="results/compare/{dataset}-{trait}-ukbb-{polarize}-%03d.png",
    wildcard_constraints:
        polarize="ancestral|focal|marginal",
    shell:
        "Rscript scripts/ukbb_compare.R"
        " --ukbb {input.ukbb}"
        " --pheno {input.pheno}"
        " --palm {input.palm}"
        " --polarize {wildcards.polarize}"
        " --output {params.png}"


rule finngen_filter_gwas:
    """
    Filter the FinnGen GWAS file to retain only the SNPs with a significant selection test in the CLUES analysis
    """
    input:
        finngen="data/finngen/gwas/finngen_R8_{pheno}.significant.tsv.bgz",
        palm="results/palm/{dataset}-{trait}-palm_report_prs.tsv",
    output:
        tsv=temp("results/compare/finngen/{dataset}-{trait}.{pheno}.significant.tsv"),
    shell:
        "Rscript scripts/finngen_filter_gwas.R"
        " --pheno {wildcards.pheno}"
        " --finngen {input.finngen}"
        " --palm {input.palm}"
        " --output {output.tsv}"


def finngen_merge_all_filtered_phenotypes_input(wildcards):
    """
    Return a list of all the phenotype association files for the FinnGen GWAS
    """
    # noinspection PyUnresolvedReferences
    pheno_tsv = checkpoints.finngen_phenotypes.get().output.tsv

    phenotypes = pd.read_table(pheno_tsv)["phenocode"].unique()

    return expand("results/compare/finngen/{dataset}-{trait}.{pheno}.significant.tsv", pheno=phenotypes, **wildcards)


rule finngen_merge_all_filtered_phenotypes:
    """
    Merge all the filtered FinnGen association files
    """
    input:
        finngen_merge_all_filtered_phenotypes_input,
    output:
        tsv="results/compare/finngen/{dataset}-{trait}.significant.tsv.gz",
    shell:
        "head -n1 {input[0]} | gzip -c > {output.tsv} && "
        "tail -n +2 -q {input} | gzip -c >> {output.tsv}"


rule finngen_compare:
    """
    Compare trajectories between trait associated SNPs and all other traits in the FinnGen
    """
    input:
        finngen="results/compare/finngen/{dataset}-{trait}.significant.tsv.gz",
        pheno="data/finngen/finngen_R8_manifest.tsv",
        palm="results/palm/{dataset}-{trait}-palm_report_prs.tsv",
    output:
        png="results/compare/{dataset}-{trait}-finngen-{polarize}-001.png",
    params:
        # plotting script produces multiple PNG files (based on number of FinnGen traits)
        png="results/compare/{dataset}-{trait}-finngen-{polarize}-%03d.png",
    wildcard_constraints:
        polarize="ancestral|focal|marginal",
    shell:
        "Rscript scripts/finngen_compare.R"
        " --finngen {input.finngen}"
        " --pheno {input.pheno}"
        " --palm {input.palm}"
        " --polarize {wildcards.polarize}"
        " --output {params.png}"
