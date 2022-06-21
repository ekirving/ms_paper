#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"


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
        "bcftools view --threads {threads} --samples-file {input.list} --output-type u {input.vcf} | "
        "bcftools annotate --threads {threads} --annotations {input.dbsnp} -c ID --output-type z --output {output.vcf} && "
        "bcftools index --tbi --threads {threads} {output.vcf}"


rule plink_convert:
    """
    Convert the subset VCF into plink format
    """
    input:
        vcf="data/1000g/vcf/1000G_phase3-chr{chr}-FIN_GBR_TSI.vcf.gz",
        tbi="data/1000g/vcf/1000G_phase3-chr{chr}-FIN_GBR_TSI.vcf.gz.tbi",
    output:
        bed="data/1000g/plink/1000G_phase3-chr{chr}-FIN_GBR_TSI.bed",
        bim="data/1000g/plink/1000G_phase3-chr{chr}-FIN_GBR_TSI.bim",
        fam="data/1000g/plink/1000G_phase3-chr{chr}-FIN_GBR_TSI.fam",
    log:
        log="data/1000g/plink/1000G_phase3-chr{chr}-FIN_GBR_TSI.log",
    params:
        prefix="data/1000g/plink/1000G_phase3-chr{chr}-FIN_GBR_TSI",
    shell:
        "plink --make-bed --vcf {input.vcf} --out {params.prefix} &> {log}"
