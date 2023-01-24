#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import gzip

"""
FinnGen GWAS

Rules for downloading and processing the FinnGen GWAS summary statistics.

DF8 - Dec 1 2022
Total sample size: 342,499 (190,879 females and 151,620 males)

Total number of variants analyzed: 20,175,454 variants

Number of disease endpoints (phenotypes) available: 2,202 endpoints*

Data released to the partners: Q4 2021

Public release: Dec 1, 2022

See https://www.finngen.fi/en/access_results
"""

import pandas as pd

# genome-wide significance threshold (see https://www.nature.com/articles/ejhg2015269)
FINNGEN_PVALUE = "5e-8"

# parameters for PLINK clumping
PLINK_CLUMP_PVAL = 5e-8
PLINK_CLUMP_R2 = 0.05
PLINK_CLUMP_KB = 250


wildcard_constraints:
    pheno="T1D",


checkpoint finngen_phenotypes:
    """
    Download the FinnGen phenotype manifest file

    see https://console.cloud.google.com/storage/browser/finngen-public-data-r8/summary_stats
    """
    output:
        tsv="data/finngen/finngen_R8_manifest.tsv",
    shell:
        "wget --quiet -O {output.tsv} https://storage.googleapis.com/finngen-public-data-r8/summary_stats/R8_manifest.tsv"


rule finngen_gwas:
    """
    Download the FinnGen GWAS summary statistics for a specific phenotype

    see https://console.cloud.google.com/storage/browser/finngen-public-data-r8/summary_stats
    """
    input:
        tsv="data/finngen/finngen_R8_manifest.tsv",
    output:
        bgz="data/finngen/gwas/finngen_R8_{pheno}.tsv.bgz",
    shell:
        "awk -F'\\t' '$1==\"{wildcards.pheno}\" {{print $7}}' {input.tsv} | "
        "xargs wget --quiet -O {output.bgz} && gunzip --test {output.bgz}"


rule finngen_gwas_significant:
    """
    Extract the genome-wide significant GWAS variants from the summary stats (`pval` < 5e-8)
    """
    input:
        bgz="data/finngen/gwas/finngen_R8_{pheno}.tsv.bgz",
    output:
        bgz="data/finngen/gwas/finngen_R8_{pheno}.significant.tsv.bgz",
    params:
        col=lambda wildcards, input: gzip.open(input[0], "r").readline().decode().strip().split("\t").index("pval") + 1,
    threads: 4
    shell:
        "bgzip --decompress --stdout --threads {threads} {input.bgz} | "
        "awk -F'\\t' 'NR==1 || ${params.col} < {FINNGEN_PVALUE} {{print $0}}' | "
        "bgzip --stdout --threads {threads} > {output.bgz}"


def finngen_all_phenotypes(wildcards):
    """
    Return a list of all the phenotype association files for the FinnGen GWAS
    """
    # noinspection PyUnresolvedReferences
    pheno_tsv = checkpoints.finngen_phenotypes.get().output.tsv

    phenotypes = pd.read_table(pheno_tsv)["phenocode"].unique()

    return expand("data/finngen/gwas/finngen_R8_{pheno}.significant.tsv.bgz", pheno=phenotypes)


rule finngen_download_all:
    """
    Download all the phenotype association files for the FinnGen GWAS
    """
    input:
        finngen_all_phenotypes,
    output:
        touch("data/finngen/download.done"),


rule intersect_finngen_gwas:
    """
    Find the intersections between the the FinnGen GWAS and MS
    """
    input:
        gwas1="data/finngen/gwas/finngen_R8_{pheno}.significant.tsv.bgz",
        gwas2="data/targets/gwas_ms-full.tsv",
    output:
        tsv="data/finngen/clump/finngen_R8_{pheno}.significant.ms-full.tsv",
    shell:
        "Rscript scripts/intersect_finngen_gwas.R"
        " --gwas1 {input.gwas1}"
        " --gwas2 {input.gwas2}"
        " --output {output.tsv}"


rule plink_clump_finngen_chr:
    """
    Perform LD-based clumping using PLINK, for a single chromosome
    """
    input:
        bed="data/1000g/plink/1000G_phase3-chrALL-FIN_GBR_TSI.bed",
        bim="data/1000g/plink/1000G_phase3-chrALL-FIN_GBR_TSI.bim",
        fam="data/1000g/plink/1000G_phase3-chrALL-FIN_GBR_TSI.fam",
        tsv="data/finngen/clump/finngen_R8_{pheno}.significant.ms-full.tsv",
    output:
        clump="data/finngen/clump/finngen_R8_{pheno}.significant.ms-full.clumped",
        nosex=temp("data/finngen/clump/finngen_R8_{pheno}.significant.ms-full.nosex"),
    log:
        log="data/finngen/clump/finngen_R8_{pheno}.significant.ms-full.log",
    params:
        out="data/finngen/clump/finngen_R8_{pheno}.significant.ms-full",
    shell:
        "plink"
        " --bed {input.bed}"
        " --bim {input.bim}"
        " --fam {input.fam}"
        " --clump {input.tsv}"
        " --clump-snp-field rsid"
        " --clump-field pval"
        " --clump-p1 {PLINK_CLUMP_PVAL}"
        " --clump-r2 {PLINK_CLUMP_R2}"
        " --clump-kb {PLINK_CLUMP_KB}"
        " --out {params.out} &> {log}"


rule apply_finngen_clumping:
    """
    Apply the clumped list of statistically independent markers for each FinnGen trait
    """
    input:
        gwas="data/finngen/clump/finngen_R8_{pheno}.significant.ms-full.tsv",
        clump="data/finngen/clump/finngen_R8_{pheno}.significant.ms-full.clumped",
    output:
        gwas="data/targets/gwas_{pheno}-r0.05-kb250.tsv",
        full="data/targets/gwas_{pheno}-full.tsv",
    shell:
        "Rscript scripts/apply_finngen_clumping.R"
        " --gwas {input.gwas}"
        " --clump {input.clump}"
        " --output1 {output.gwas}"
        " --output2 {output.full}"
