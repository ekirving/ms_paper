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


checkpoint palm_metadata:
    """
    Convert the GWAS metadata into PALM input format, and replace any missing SNPs with proxy-SNPs in high LD
    """
    input:
        gwas="data/targets/gwas_{trait}.tsv",
        ld="data/targets/gwas_{trait}_ld.tsv.gz",
        sites="data/sites/{dataset}_sites.tsv.gz",
    output:
        tsv="data/targets/gwas_{trait}_{dataset}_palm.tsv",
    shell:
        "Rscript scripts/palm_metadata.R"
        " --gwas {input.gwas}"
        " --ld {input.ld}"
        " --sites {input.sites}"
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

    # make sure that all the SNP metadata has been built, as we rely on this later
    json = expand("data/metadata/GRCh37/{rsid}.json", rsid=meta["rsid"])

    return {"tsv": "data/targets/gwas_{trait}_{dataset}_palm.tsv", "quad_fit": files, "json": json}


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
        " --maxp 1"
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


rule palm_report_prs:
    input:
        tsv="results/palm/{dataset}-{ancestry}-{trait}-palm_report.tsv",
    output:
        tsv="results/palm/{dataset}-{ancestry}-{trait}-palm_report_prs.tsv",
    shell:
        "Rscript scripts/palm_report_prs.R"
        " --data {input.tsv}"
        " --dataset {wildcards.dataset}"
        " --ancestry {wildcards.ancestry}"
        " --output {output.tsv}"


rule palm_plot_trajectory:
    """
    Plot the PALM trajectory as a raster of the joint posterior densities of all SNPs
    """
    input:
        tsv="results/palm/{dataset}-{ancestry}-{trait}-palm_report.tsv",
        json="results/palm/{dataset}-{ancestry}-{trait}-palm.json",
    output:
        png="results/palm/{dataset}-{ancestry}-{trait}-palm_trajectory.png",
    resources:
        mem_mb=64 * 1024,
    shell:
        "Rscript scripts/palm_plot_trajectory.R"
        " --tsv {input.tsv}"
        " --json {input.json}"
        " --trait {wildcards.trait}"
        " --dataset {wildcards.dataset}"
        " --ancestry {wildcards.ancestry}"
        " --output {output.png}"


rule palm_plot_lines:
    """
    Plot the PALM trajectory as a stacked line plot of each individual SNP trajectory
    """
    input:
        tsv="results/palm/{dataset}-{ancestry}-{trait}-palm_report.tsv",
        json="results/palm/{dataset}-{ancestry}-{trait}-palm.json",
    output:
        png="results/palm/{dataset}-{ancestry}-{trait}-palm_lines-{sort}.png",
    wildcard_constraints:
        sort="pval|prs",
    shell:
        "Rscript scripts/palm_plot_lines_{wildcards.sort}.R"
        " --tsv {input.tsv}"
        " --json {input.json}"
        " --trait {wildcards.trait}"
        " --dataset {wildcards.dataset}"
        " --ancestry {wildcards.ancestry}"
        " --output {output.png}"


rule palm_plot_delta_prs:
    """
    Plot the distribution of the delta PRS for all SNPs in all ancestries  
    """
    input:
        all_json="results/palm/{dataset}-ALL-{trait}-palm.json",
        ana_json="results/palm/{dataset}-ANA-{trait}-palm.json",
        chg_json="results/palm/{dataset}-CHG-{trait}-palm.json",
        ehg_json="results/palm/{dataset}-EHG-{trait}-palm.json",
        whg_json="results/palm/{dataset}-WHG-{trait}-palm.json",
        all_tsv="results/palm/{dataset}-ALL-{trait}-palm_report_prs.tsv",
        ana_tsv="results/palm/{dataset}-ANA-{trait}-palm_report_prs.tsv",
        chg_tsv="results/palm/{dataset}-CHG-{trait}-palm_report_prs.tsv",
        ehg_tsv="results/palm/{dataset}-EHG-{trait}-palm_report_prs.tsv",
        whg_tsv="results/palm/{dataset}-WHG-{trait}-palm_report_prs.tsv",
    output:
        png="results/palm/{dataset}-{trait}-delta_prs.png",
    wildcard_constraints:
        sort="pval|prs",
    shell:
        "Rscript scripts/palm_plot_delta_prs.R"
        " --trait {wildcards.trait}"
        " --dataset {wildcards.dataset}"
        " --all-json {input.all_json}"
        " --ana-json {input.ana_json}"
        " --chg-json {input.chg_json}"
        " --ehg-json {input.ehg_json}"
        " --whg-json {input.whg_json}"
        " --all-tsv {input.all_tsv}"
        " --ana-tsv {input.ana_tsv}"
        " --chg-tsv {input.chg_tsv}"
        " --ehg-tsv {input.ehg_tsv}"
        " --whg-tsv {input.whg_tsv}"
        " --output {output.png}"
