#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

"""
CLUES analyses

Rules for inference of allele frequency trajectories and selection coefficients with CLUES.

See https://github.com/35ajstern/clues/
"""

ANCESTRIES = ["ALL", "ANA", "CHG", "WHG", "EHG"]

CHROM = r"(\d+|X|Y|MT)"


wildcard_constraints:
    # `chr:pos:ref:alt`
    variant="\d+:\d+:[ACGT]:[ACGT]",
    rsid="rs\d+",


rule clues_ancient_samples:
    input:
        vcf=lambda wildcards: config["samples"][wildcards.dataset]["genotypes"],
    output:
        anc="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.ancient",
        frq="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.modfreq",
    log:
        log="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.ancient.log",
    shell:
        "python scripts/clues_ancient_samples.py"
        " --vcf {input.vcf}"
        " --dataset {wildcards.dataset}"
        " --variant {wildcards.variant}"
        " --ancestry {wildcards.ancestry}"
        " --gen-time {config[gen_time]}"
        " --mod-freq {output.frq}"
        " --output {output.anc} 2> {log}"


rule clues_time_bins:
    output:
        bins="data/clues/{dataset}-time.bins",
    shell:
        "python scripts/clues_time_bins.py"
        " --dataset {wildcards.dataset}"
        " --gen-time {config[gen_time]}"
        " --output {output}"


rule clues_inference_ancient:
    input:
        coal="data/relate/1000G_phase3-FIN_GBR_TSI-popsize.coal",
        bins="data/clues/{dataset}-time.bins",
        anct="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.ancient",
        freq="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.modfreq",
    output:
        epochs="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.epochs.npy",
        freqs="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.freqs.npy",
        post="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.post.npy",
        quad="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.quad_fit.npy",
    log:
        "results/clues/{rsid}/{dataset}-{variant}-{ancestry}.epochs.log",
    params:
        out="results/clues/{rsid}/{dataset}-{variant}-{ancestry}",
        daf=lambda wildcards, input: open(str(input.freq)).read().strip(),
        flg=lambda wildcards: "ancientSamps" if wildcards.ancestry == "ALL" else "ancientHaps",
    shell:
        "python bin/clues/inference.py"
        " --lik"
        " --popFreq {params.daf}"
        " --coal {input.coal}"
        " --{params.flg} {input.anct}"
        " --timeBins {input.bins}"
        " --out {params.out} &> {log}"


rule clues_parse_log:
    input:
        log="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.epochs.log",
    output:
        json="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.json",
    shell:
        "python scripts/clues_parse_log.py"
        " --rsid {wildcards.rsid}"
        " --ancestry {wildcards.ancestry}"
        " --log {input.log}"
        " --out {output.json}"


rule clues_plot_trajectory:
    input:
        epch="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.epochs.npy",
        freq="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.freqs.npy",
        post="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.post.npy",
        json="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.json",
        label="data/labels/{dataset}/{dataset}-{rsid}-label.json",
    output:
        png="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.png",
    log:
        log="results/clues/{rsid}/{dataset}-{variant}-{ancestry}.png.log",
    params:
        input="results/clues/{rsid}/{dataset}-{variant}-{ancestry}",
        output="results/clues/{rsid}/{dataset}-{variant}-{ancestry}",
    shell:
        "python scripts/clues_plot_trajectory.py"
        " --gen-time {config[gen_time]}"
        " --params {input.json}"
        " --label {input.label}"
        " --ancestry {wildcards.ancestry}"
        " --ext png"
        " {params.input}"
        " {params.output} 2> {log}"
