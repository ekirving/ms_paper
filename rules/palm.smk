#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd


"""
Rules for detecting polygenic selection with PALM  
"""


wildcard_constraints:
    trait="[^_]+",


rule reference_metadata:
    """
    Extract the REF, ALT and ancestral alleles from the dataset VCF
    """
    input:
        vcf=lambda wildcards: config["samples"][wildcards.dataset]["genotypes"],
        tsv="data/targets/gwas_{trait}.tsv",
    output:
        pos=temp("data/targets/gwas_{trait}_{dataset}_positions.tsv"),
        var=temp("data/targets/gwas_{trait}_{dataset}_variants.tsv"),
    shell:
        "awk 'NR>1 {{ print $1\"\\t\"$2 }}' {input.tsv} > {output.pos} && "
        "bcftools view --regions-file {output.pos} {input.vcf} | "
        r" bcftools query --format '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AA\n' | "
        r" sed 's/|||//' | tr [a-z] [A-Z] > {output.var}"


checkpoint palm_metadata:
    """
    Convert the GWAS metadata into PALM input format
    """
    input:
        # TODO make this a config call
        tsv="data/targets/gwas_{trait}.tsv",
        var="data/targets/gwas_{trait}_{dataset}_variants.tsv",
    output:
        tsv="data/targets/gwas_{trait}_{dataset}_palm.tsv",
    shell:
        "Rscript scripts/palm_metadata.R"
        " --gwas {input.tsv}"
        " --variants {input.var}"
        " --output {output.tsv}"


def clues_quad_fit(wildcards):
    """
    Resolve the path to the CLUES likelihood quadratic fit
    """
    # noinspection PyUnresolvedReferences
    meta_tsv = checkpoints.palm_metadata.get(**wildcards).output.tsv

    snp = pd.read_table(meta_tsv).set_index("ld_block", drop=False).loc[int(wildcards.block)]

    rsid = snp["rsid"]
    dataset = wildcards.dataset
    ancestry = wildcards.ancestry

    chr, pos, ref, alt = snp["variant"].split(":")
    der = snp["derived_allele"]
    anc = snp["ancestral_allele"]

    # CLUES expects the variant name to be polarized by the ancestral allele
    variant = f"{chr}:{pos}:{anc}:{der}"

    return [
        f"results/clues/{rsid}/{dataset}-{variant}-{ancestry}.quad_fit.npy",
        f"results/clues/{rsid}/{dataset}-{variant}-{ancestry}.json",
    ]


rule palm_organise_clues:
    """
    Copy the CLUES likelihood quadratic fit into the PALM directory structure
    """
    input:
        clues_quad_fit,
    output:
        quad="results/palm/{dataset}/{ancestry}/{trait}/ld_{block}/bp{pos}.quad_fit.npy",
    shell:
        "cp {input[0]} {output}"


def palm_quad_fit(wildcards):
    """
    Return a list of the `.quad_fit` files for each SNP associated with the current trait
    """
    # noinspection PyUnresolvedReferences
    meta_tsv = checkpoints.palm_metadata.get(**wildcards).output.tsv
    meta = pd.read_table(meta_tsv)

    dataset = wildcards.dataset
    ancestry = wildcards.ancestry
    trait = wildcards.trait

    files = []
    for _, snp in meta.iterrows():
        block, pos = snp["ld_block"], snp["variant"].split(":")[1]
        files.append(f"results/palm/{dataset}/{ancestry}/{trait}/ld_{block}/bp{pos}.quad_fit.npy")

    return {"tsv": "data/targets/gwas_{trait}_{dataset}_palm.tsv", "quad_fit": files}


# noinspection PyUnresolvedReferences
rule palm_single_trait:
    """
    Run PALM in single-trait mode
    """
    input:
        unpack(palm_quad_fit),
    output:
        txt="results/palm/{dataset}-{ancestry}-{trait}-palm.txt",
    log:
        log="results/palm/{dataset}-{ancestry}-{trait}-palm.log",
    params:
        dir="results/palm/{dataset}/{ancestry}/{trait}/",
    shell:
        "python bin/palm/palm.py"
        " --traitDir {params.dir}"
        " --metadata {input.tsv}"
        " --maxp 5e-8"
        " --B 1000"
        " 1> {output.txt}"
        " 2> {log}"


rule palm_parse_txt:
    """
    Convert the text output from PALM into a JSON file
    """
    input:
        palm="results/palm/{dataset}-{ancestry}-{trait}-palm.txt",
    output:
        json="results/palm/{dataset}-{ancestry}-{trait}-palm.json",
    shell:
        "python scripts/palm_parse_txt.py"
        " --dataset {wildcards.dataset}"
        " --ancestry {wildcards.ancestry}"
        " --trait {wildcards.trait}"
        " --palm {input.palm}"
        " --out {output.json}"


rule palm_report:
    input:
        tsv="data/targets/gwas_{trait}_{dataset}_palm.tsv",
        json="results/palm/{dataset}-{ancestry}-{trait}-palm.json",
    output:
        tsv="results/palm/{dataset}-{ancestry}-{trait}-palm_report.tsv",
    shell:
        "python scripts/palm_report.py"
        " --data {input.tsv}"
        " --dataset {wildcards.dataset}"
        " --ancestry {wildcards.ancestry}"
        " --output {output.tsv}"


rule palm_plot_trajectory:
    """
    Plot the PALM trajectory as a rater or as stacked lines
    """
    input:
        tsv="results/palm/{dataset}-{ancestry}-{trait}-palm_report.tsv",
        json="results/palm/{dataset}-{ancestry}-{trait}-palm.json",
    output:
        png="results/palm/{dataset}-{ancestry}-{trait}-palm_{type}.png",
    wildcard_constraints:
        type="trajectory|lines",
    resources:
        mem_mb=40 * 1024,
    shell:
        "Rscript scripts/palm_plot_{wildcards.type}.R"
        " --tsv {input.tsv}"
        " --json {input.json}"
        " --trait {wildcards.trait}"
        " --dataset {wildcards.dataset}"
        " --ancestry {wildcards.ancestry}"
        " --output {output.png}"
