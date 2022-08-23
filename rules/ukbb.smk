#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2021, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import gzip

"""
NealeLab UKBB GWAS

Rules for downloading and processing the Neale Lab (Round 2) UK Biobank GWAS summary statistics.

See https://www.nealelab.is/uk-biobank
"""

import pandas as pd

# genome-wide significance threshold (see https://www.nature.com/articles/ejhg2015269)
UKBB_PVALUE = "5e-8"


wildcard_constraints:
    sex="both_sexes|female|male",


rule ukbb_nealelab_variants:
    """
    Download the NealeLab UKBB Variant manifest file

    see https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=227859291
    """
    input:
        tsv=ancient("data/ukbb/nealelab/UKBB_GWAS_Imputed_v3_Manifest_201807.tsv"),
    output:
        bgz="data/ukbb/nealelab/variants.tsv.bgz",
    shell:
        """awk -F'\\t' '$5=="variants.tsv.bgz" {{print $7}}' {input.tsv} | xargs wget --quiet -O {output.bgz}"""


checkpoint ukbb_nealelab_phenotypes:
    """
    Download the NealeLab UKBB Phenotype manifest file

    see https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=227859291
    """
    input:
        tsv=ancient("data/ukbb/nealelab/UKBB_GWAS_Imputed_v3_Manifest_201807.tsv"),
    output:
        bgz="data/ukbb/nealelab/phenotypes.{sex}.tsv.bgz",
    shell:
        """awk -F'\\t' '$5=="phenotypes.{wildcards.sex}.tsv.bgz" {{print $7}}' {input.tsv} | xargs wget --quiet -O {output.bgz}"""


rule ukbb_nealelab_gwas_md5:
    """
    Fetch the md5sum for the Per-phenotype file
    """
    input:
        tsv=ancient("data/ukbb/nealelab/UKBB_GWAS_Imputed_v3_Manifest_201807.tsv"),
    output:
        md5="data/ukbb/nealelab/gwas/{pheno}.gwas.imputed_v3.{sex}.tsv.bgz.md5",
    params:
        bgz="data/ukbb/nealelab/gwas/{pheno}.gwas.imputed_v3.{sex}.tsv.bgz",
    shell:
        """awk -F'\\t' '$1=="{wildcards.pheno}" && $4=="{wildcards.sex}" {{print $9" {params.bgz}"}}' {input.tsv} | head -n1 > {output.md5}"""


rule ukbb_nealelab_gwas:
    """
    Download the NealeLab UKBB GWAS summary statistics for a specific phenotype

    see https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?ts=5b5f17db#gid=227859291
    """
    input:
        tsv=ancient("data/ukbb/nealelab/UKBB_GWAS_Imputed_v3_Manifest_201807.tsv"),
        md5="data/ukbb/nealelab/gwas/{pheno}.gwas.imputed_v3.{sex}.tsv.bgz.md5",
    output:
        bgz="data/ukbb/nealelab/gwas/{pheno}.gwas.imputed_v3.{sex}.tsv.bgz",
    shell:
        """awk -F'\\t' '$1=="{wildcards.pheno}" && $4=="{wildcards.sex}" {{print $7}}' {input.tsv} | head -n1 | """
        """xargs wget --quiet -O {output.bgz} && md5sum --status --check {input.md5}"""


rule ukbb_nealelab_gwas_significant:
    """
    Extract the genome-wide significant GWAS variants from the meta analysis group (`pval` < 5e-8)
    """
    input:
        bgz="data/ukbb/nealelab/gwas/{pheno}.gwas.imputed_v3.{sex}.tsv.bgz",
    output:
        bgz="data/ukbb/nealelab/gwas/{pheno}.gwas.imputed_v3.{sex}.significant.tsv.bgz",
    params:
        col=lambda wildcards, input: gzip.open(input[0], "r").readline().decode().strip().split("\t").index("pval") + 1,
    threads: 4
    shell:
        "bgzip --decompress --stdout --threads {threads} {input.bgz} | "
        "awk -F'\\t' 'NR==1 || ${params.col} < {UKBB_PVALUE} {{print $0}}' | "
        "bgzip --stdout --threads {threads} > {output.bgz}"


def ukbb_nealelab_all_phenotypes(wildcards):
    """
    Return a list of all the phenotype association files for the NealeLab UKBB GWAS
    """
    # noinspection PyUnresolvedReferences
    pheno_tsv = checkpoints.ukbb_nealelab_phenotypes.get(sex=wildcards.sex).output.bgz

    phenotypes = pd.read_table(pheno_tsv, compression="gzip")["phenotype"].unique()

    return expand(
        "data/ukbb/nealelab/gwas/{pheno}.gwas.imputed_v3.{sex}.significant.tsv.bgz", pheno=phenotypes, sex=wildcards.sex
    )


rule ukbb_nealelab_download_all:
    """
    Download all the phenotype association files for the NealeLab UKBB GWAS
    """
    input:
        ukbb_nealelab_all_phenotypes,
    output:
        touch("data/ukbb/nealelab/download.{sex}.done"),
