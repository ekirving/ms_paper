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

MHC_CHROM = 6
MHC_START = 28477797
MHC_FINISH = 33448354

ANCESTRIES = ["ALL", "ANA", "CHG", "WHG", "EHG"]

# minimum LD to clump SNPs into the same LD block for the joint PALM analysis
MIN_LD = 0.3

# minimum fraction of genetic overlap between selected variants in the focal trait and a correlated trait in UKBB or FinnGen
MIN_GWAS_OVERLAP = 0.2


checkpoint palm_metadata_single_trait:
    """
    Convert the GWAS metadata into PALM input format, and optionally replace any missing SNPs with proxy-SNPs in high LD
    """
    input:
        gwas="data/targets/gwas_{trait}.tsv",
        # ld="data/targets/gwas_{trait}_ld.tsv.gz",
        sites="data/sites/{dataset}_sites.tsv.gz",
    output:
        tsv="data/targets/gwas_{trait}_{dataset}_palm.tsv",
    log:
        log="data/targets/gwas_{trait}_{dataset}_palm.log",
    shell:
        "Rscript scripts/palm_metadata_single_trait.R"
        " --gwas {input.gwas}"
        " --sites {input.sites}"
        " --output {output.tsv} &> {log}"


checkpoint palm_metadata_multi_trait:
    """
    Merge the two sets of GWAS metadata using pairwise LD in the modern reference panel

    To do this, we need both the independent markers for each marginal trait, as well as the full set of summary 
    statistics for both traits, so we can output beta scores and p-values for all independent SNPs used in both traits
    """
    input:
        gwas1_ind="data/targets/gwas_{trait1}-r0.05-kb250_{dataset}_palm.tsv",
        gwas2_ind="data/targets/gwas_{trait2}-r0.05-kb250_{dataset}_palm.tsv",
        gwas1_all="data/targets/gwas_{trait1}-full_{dataset}_palm.tsv",
        gwas2_all="data/targets/gwas_{trait2}-full_{dataset}_palm.tsv",
        ld="data/targets/gwas_{trait1}-r0.05-kb250_ld.tsv.gz",
    output:
        tsv="data/targets/gwas_{trait1}~{trait2}_{dataset}_palm.tsv",
    log:
        log="data/targets/gwas_{trait1}~{trait2}_{dataset}_palm.log",
    shell:
        "Rscript scripts/palm_metadata_multi_trait.R"
        " --trait1 {wildcards.trait1}"
        " --trait2 {wildcards.trait2}"
        " --gwas1-ind {input.gwas1_ind}"
        " --gwas2-ind {input.gwas2_ind}"
        " --gwas1-all {input.gwas1_all}"
        " --gwas2-all {input.gwas2_all}"
        " --ld {input.ld}"
        " --min-ld {MIN_LD}"
        " --output {output.tsv} &> {log}"


def clues_quad_fit(wildcards):
    """
    Resolve the path to the CLUES likelihood quadratic fit
    """
    if hasattr(wildcards, "trait1"):
        # noinspection PyUnresolvedReferences
        meta_tsv = checkpoints.palm_metadata_multi_trait.get(**wildcards).output.tsv
    else:
        # noinspection PyUnresolvedReferences
        meta_tsv = checkpoints.palm_metadata_single_trait.get(**wildcards).output.tsv

    meta = pd.read_table(meta_tsv)
    snp = (
        meta.loc[(meta["ld_block"] == int(wildcards.block)) & (meta["pos"] == int(wildcards.pos))]
        .to_dict(orient="records")
        .pop()
    )

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
    wildcard_constraints:
        # relax the constraint so we can match paired traits
        trait="[^/]+",
    shell:
        "cp {input[0]} {output}"


def palm_quad_fit(wildcards):
    """
    Return a list of the `.quad_fit` files for each SNP associated with the current trait
    """
    if hasattr(wildcards, "trait1"):
        # noinspection PyUnresolvedReferences
        meta_tsv = checkpoints.palm_metadata_multi_trait.get(**wildcards).output.tsv
    else:
        # noinspection PyUnresolvedReferences
        meta_tsv = checkpoints.palm_metadata_single_trait.get(**wildcards).output.tsv

    meta = pd.read_table(meta_tsv)

    dataset = wildcards.dataset
    ancestry = wildcards.ancestry
    trait = wildcards.trait if hasattr(wildcards, "trait") else "~".join([wildcards.trait1, wildcards.trait2])

    files = []
    for _, snp in meta.iterrows():
        block, pos = snp["ld_block"], snp["variant"].split(":")[1]
        files.append(f"results/palm/{dataset}/{ancestry}/{trait}/ld_{block}/bp{pos}.quad_fit.npy")

    # make sure that all the SNP metadata has been built, as we rely on this later
    json = expand("data/metadata/GRCh37/{rsid}.json", rsid=meta["rsid"])

    return {"quad_fit": files, "json": json}


# noinspection PyUnresolvedReferences
rule palm_single_trait:
    """
    Run PALM in single-trait mode
    """
    input:
        unpack(palm_quad_fit),
        tsv="data/targets/gwas_{trait}_{dataset}_palm.tsv",
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
        " --B 1000"
        " 1> {output.txt}"
        " 2> {log}"


# noinspection PyUnresolvedReferences
rule palm_multi_trait:
    """
    Run PALM in multi-trait mode
    """
    input:
        unpack(palm_quad_fit),
        tsv="data/targets/gwas_{trait1}~{trait2}_{dataset}_palm.tsv",
    output:
        txt="results/palm/{dataset}-{ancestry}-{trait1}~{trait2}-palm.txt",
    log:
        log="results/palm/{dataset}-{ancestry}-{trait1}~{trait2}-palm.log",
    params:
        dir="results/palm/{dataset}/{ancestry}/{trait1}~{trait2}/",
    shell:
        "python bin/palm/palm.py"
        " --traitDir {params.dir}"
        " --metadata {input.tsv}"
        " --B 1000"
        " --maxp 5e-8"
        " --traits {wildcards.trait1},{wildcards.trait2}"
        " 1> {output.txt}"
        " 2> {log}"


ruleorder: palm_multi_trait > palm_single_trait


rule palm_parse_single_trait:
    """
    Convert the text output from PALM into a JSON file
    """
    input:
        palm="results/palm/{dataset}-{ancestry}-{trait}-palm.txt",
    output:
        json="results/palm/{dataset}-{ancestry}-{trait}-palm.json",
    shell:
        "python scripts/palm_parse_single_trait.py"
        " --dataset {wildcards.dataset}"
        " --ancestry {wildcards.ancestry}"
        " --trait {wildcards.trait}"
        " --palm {input.palm}"
        " --out {output.json}"


rule palm_parse_multi_trait:
    """
    Convert the text output from multi trait PALM into a JSON file
    """
    input:
        palm="results/palm/{dataset}-{ancestry}-{trait1}~{trait2}-palm.txt",
    output:
        json="results/palm/{dataset}-{ancestry}-{trait1}~{trait2}-palm.json",
    shell:
        "python scripts/palm_parse_multi_trait.py"
        " --dataset {wildcards.dataset}"
        " --ancestry {wildcards.ancestry}"
        " --trait1 {wildcards.trait1}"
        " --trait2 {wildcards.trait2}"
        " --palm {input.palm}"
        " --out {output.json}"


rule palm_ancestry_report:
    """
    Make a PALM report, summarising all the models for this ancestry
    """
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


rule palm_ancestry_report_prs:
    """
    Calculate the change in PRS for each SNP and add it to the report
    """
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


rule palm_report_prs:
    """
    Concatenate all the ancestry specific PALM reports together
    """
    input:
        expand("results/palm/{dataset}-{ancestry}-{trait}-palm_report_prs.tsv", ancestry=ANCESTRIES, allow_missing=True),
    output:
        tsv="results/palm/{dataset}-{trait}-palm_report_prs.tsv",
    shell:
        "head -n1 {input[0]} > {output.tsv} && "
        "tail -n +2 -q {input} >> {output.tsv}"


ruleorder: palm_ancestry_report_prs > palm_report_prs


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


rule palm_plot_scatter:
    """
    Plot a scatter plot of -log10(pval) vs delra_prs
    """
    input:
        tsv="results/palm/{dataset}-{trait}-palm_report_prs.tsv",
    output:
        png="results/palm/{dataset}-{trait}-scatter.png",
    shell:
        "Rscript scripts/palm_plot_scatter.R "
        " --palm {input.tsv}"
        " --output {output.png}"


def all_overlapping_traits(wildcards):
    """
    Run all the overlapping traits in UKBB and FinnGen
    """
    files = []

    # noinspection PyUnresolvedReferences
    ukbb = pd.read_table(checkpoints.ukbb_compare.get(**wildcards, polarize="marginal").output.tsv)
    ukbb = ukbb[ukbb["frac_snps"] > MIN_GWAS_OVERLAP]

    # truncate the trait name
    trait1 = wildcards.trait.split("-")[0]

    # run all the overlapping UKBB traits
    files += expand(
        "results/palm/{dataset}-{ancestry}-{trait1}~{pheno}-ukbb-palm.json",
        ancestry=ANCESTRIES,
        trait1=trait1,
        pheno=ukbb["phenotype"],
        allow_missing=True,
    )

    # noinspection PyUnresolvedReferences
    finngen = pd.read_table(checkpoints.finngen_compare.get(**wildcards, polarize="marginal").output.tsv)
    finngen = finngen[finngen["frac_snps"] > MIN_GWAS_OVERLAP]

    # run all the overlapping FinnGen traits
    files += expand(
        "results/palm/{dataset}-{ancestry}-{trait1}~{pheno}-finngen-palm.json",
        ancestry=ANCESTRIES,
        trait1=trait1,
        pheno=finngen["phenotype"],
        allow_missing=True,
    )

    return {"models": files}


rule palm_all_overlapping_traits_report:
    """
    Make a PALM multi-trait report for all overlapping traits in UKBB and FinnGEN  
    """
    input:
        unpack(all_overlapping_traits),
        ukbb="data/ukbb/nealelab/phenotypes.both_sexes.tsv.bgz",
        finngen="data/finngen/finngen_R8_manifest.tsv",
    params:
        models=lambda wildcards, input: [f"--model {model}" for model in input.models],
    output:
        tsv="results/palm/{dataset}-{trait}-palm_report_multi_trait.tsv",
    shell:
        "python scripts/palm_report_multi_trait.py"
        " {params.models}"
        " --ukbb {input.ukbb}"
        " --finngen {input.finngen}"
        " --output {output.tsv}"
