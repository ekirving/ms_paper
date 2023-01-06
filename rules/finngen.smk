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


wildcard_constraints:
    pheno="[^.]+",


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
