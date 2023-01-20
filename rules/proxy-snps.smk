#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd


"""
Rules for finding SNPs in high LD with missing GWAS tag SNPs
"""

TGP_FTP = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"


rule download_1000g_vcf:
    """
    Download the GRCh37 build of the 1000G phase 3 dataset
    """
    output:
        vcf=temp("data/1000g/vcf/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"),
        tbi=temp("data/1000g/vcf/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"),
    shell:
        "wget --quiet -O {output.vcf} -o /dev/null {TGP_FTP}/ALL.chr{wildcards.chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz && "
        "wget --quiet -O {output.tbi} -o /dev/null {TGP_FTP}/ALL.chr{wildcards.chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi && "
        "gunzip --test {output.vcf}"


rule download_1000g_sample_list:
    """
    Download the list of samples from the 1000G phase 3 dataset
    """
    output:
        panel="data/1000g/integrated_call_samples_v3.20130502.ALL.panel",
    shell:
        "wget --quiet -O {output.panel} -o /dev/null {TGP_FTP}/integrated_call_samples_v3.20130502.ALL.panel"


rule filter_1000g_sample_list:
    """
    Get the list of samples belonging to the FIN, GBR, and TSI populations 
    """
    input:
        panel="data/1000g/integrated_call_samples_v3.20130502.ALL.panel",
    output:
        list="data/1000g/1000G_phase3-samples-FIN_GBR_TSI.list",
    shell:
        r"grep -P '\b(FIN|GBR|TSI)\b' {input.panel} | awk '{{ print $1 }}'> {output.list}"


rule filter_1000g_populations:
    """
    Extract all the samples from the FIN, GBR, and TSI populations, and annotate variant IDs with dbSNP
    """
    input:
        list="data/1000g/1000G_phase3-samples-FIN_GBR_TSI.list",
        vcf="data/1000g/vcf/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        tbi="data/1000g/vcf/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi",
        dbsnp="data/dbsnp/GRCh37.dbSNP155.vcf.gz",
    output:
        vcf=temp("data/1000g/vcf/1000G_phase3-chr{chr}-FIN_GBR_TSI.vcf.gz"),
        tbi=temp("data/1000g/vcf/1000G_phase3-chr{chr}-FIN_GBR_TSI.vcf.gz.tbi"),
    threads: 4
    shell:
        "bcftools annotate --threads {threads} --annotations {input.dbsnp} --columns ID --output-type u {input.vcf} | "
        "bcftools view --threads {threads} --samples-file {input.list} --output-type z --output {output.vcf} && "
        "bcftools index --tbi --threads {threads} {output.vcf}"


rule plink_convert:
    """
    Convert the subset VCF into plink format
    """
    input:
        vcf="data/1000g/vcf/1000G_phase3-chr{chr}-FIN_GBR_TSI.vcf.gz",
        tbi="data/1000g/vcf/1000G_phase3-chr{chr}-FIN_GBR_TSI.vcf.gz.tbi",
    output:
        bed=temp("data/1000g/plink/1000G_phase3-chr{chr}-FIN_GBR_TSI.bed"),
        bim=temp("data/1000g/plink/1000G_phase3-chr{chr}-FIN_GBR_TSI.bim"),
        fam=temp("data/1000g/plink/1000G_phase3-chr{chr}-FIN_GBR_TSI.fam"),
    log:
        log="data/1000g/plink/1000G_phase3-chr{chr}-FIN_GBR_TSI.log",
    params:
        prefix="data/1000g/plink/1000G_phase3-chr{chr}-FIN_GBR_TSI",
    shell:
        "plink"
        " --make-bed"
        " --snps-only"
        " --set-missing-var-ids 'chr@:#'"
        " --vcf {input.vcf}"
        " --out {params.prefix} &> {log}"


rule plink_merge:
    """
    Merge all the chromosomes back together
    """
    input:
        bed=expand("data/1000g/plink/1000G_phase3-chr{chr}-FIN_GBR_TSI.bed", chr=config["chroms"]),
        bim=expand("data/1000g/plink/1000G_phase3-chr{chr}-FIN_GBR_TSI.bim", chr=config["chroms"]),
        fam=expand("data/1000g/plink/1000G_phase3-chr{chr}-FIN_GBR_TSI.fam", chr=config["chroms"]),
    output:
        bed="data/1000g/plink/1000G_phase3-chrALL-FIN_GBR_TSI.bed",
        bim="data/1000g/plink/1000G_phase3-chrALL-FIN_GBR_TSI.bim",
        fam="data/1000g/plink/1000G_phase3-chrALL-FIN_GBR_TSI.fam",
        list=temp("data/1000g/plink/1000G_phase3-chrALL-FIN_GBR_TSI.list"),
    log:
        log="data/1000g/plink/1000G_phase3-chrALL-FIN_GBR_TSI.log",
    params:
        bed="\\n".join([f"data/1000g/plink/1000G_phase3-chr{chr}-FIN_GBR_TSI" for chr in config["chroms"]]),
        out="data/1000g/plink/1000G_phase3-chrALL-FIN_GBR_TSI",
    shell:
        "printf '{params.bed}' > {output.list} && "
        "plink"
        " --make-bed"
        " --merge-list {output.list}"
        " --out {params.out} &> {log}"


rule plink_pairwise_ld:
    """
    Calculate pairwise LD between each SNP and all other SNPs within 1Mb
    """
    input:
        bed="data/1000g/plink/1000G_phase3-chrALL-FIN_GBR_TSI.bed",
        bim="data/1000g/plink/1000G_phase3-chrALL-FIN_GBR_TSI.bim",
        fam="data/1000g/plink/1000G_phase3-chrALL-FIN_GBR_TSI.fam",
    output:
        ld=temp("data/1000g/ld/1000G_phase3-chrALL-{rsid}.ld"),
    log:
        log="data/1000g/ld/1000G_phase3-chrALL-{rsid}.log",
    params:
        out="data/1000g/ld/1000G_phase3-chrALL-{rsid}",
    shell:
        "plink"
        " --bed {input.bed}"
        " --bim {input.bim}"
        " --fam {input.fam}"
        " --r2 in-phase"
        " --ld-snp {wildcards.rsid}"
        " --ld-window-kb 1000"
        " --ld-window 99999"
        " --ld-window-r2 0"
        " --out {params.out}"
        " &> {log}"


def concatenate_pairwise_ld_input(wildcards):
    """
    List all the SNPs for which we need to calculate the pair-wise LD
    """
    # noinspection PyUnresolvedReferences
    gwas_tsv = checkpoints.gwas_metadata.get(trait=wildcards.trait).output.tsv

    snp = pd.read_table(gwas_tsv)

    # drop SNPs that are not present in the 1000G callset
    missing = [
        "rs736160",
        "rs2853952",
        "rs372968977",
        "rs28637341",
        "rs1610628",
        "rs376236557",
        "rs28480108",
        "rs9689804",
        "rs9258716",
        "rs59396305",
        "rs28693951",
        "rs9277575",
        "rs67375766",
        "rs1631950",
        "rs55873403",
        "rs34341880",
    ]

    files = [f"data/1000g/ld/1000G_phase3-chrALL-{rsid}.ld" for rsid in snp["SNP"] if rsid not in missing]

    return files


rule concatenate_pairwise_ld:
    """
    Concatenate the LD tables into a single sheet
    """
    input:
        concatenate_pairwise_ld_input,
    output:
        tsv="data/targets/gwas_{trait}_ld.tsv.gz",
    shell:
        "head -n1 {input[0]} | gzip -c > {output.tsv} && "
        "tail -n +2 -q {input} | gzip -c >> {output.tsv}"


rule list_callable_sites:
    """
    Make a list of all the callable sites in the current dataset
    """
    input:
        vcf=lambda wildcards: config["samples"][wildcards.dataset]["genotypes"],
    output:
        tsv="data/sites/{dataset}_sites.tsv.gz",
    shell:
        r"bcftools query --print-header --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AA\n' {input.vcf} | "
        r"sed 's/|||//' | tr [acgt] [ACGT] | bgzip -c > {output.tsv}"
