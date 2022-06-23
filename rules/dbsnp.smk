#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import temp, expand

"""
Rules to fetch the dbSNP database
"""


rule reference_grch37_assembly_report:
    """
    Fetch the GRCh37 assembly report
    """
    output:
        txt=temp("data/GCF_000001405.25_GRCh37.p13_assembly_report.txt"),
    shell:
        "wget --quiet -O {output.txt} -o /dev/null ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt"


rule reference_dbsnp_b155:
    """
    Fetch the dbSNP database, build 155 (2021-05-13)
    """
    output:
        vcf=temp("data/dbsnp/GCF_000001405.25.gz"),
        tbi=temp("data/dbsnp/GCF_000001405.25.gz.tbi"),
    shell:
        "wget --quiet -O {output.vcf} -o /dev/null ftp://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.25.gz && "
        "wget --quiet -O {output.tbi} -o /dev/null ftp://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.25.gz.tbi"


rule reference_dbsnp_b155_reheader:
    """
    Change the chromosome names in the dbsnp VCF from `RefSeq ID` style to `UCSC-style-name`
    """
    input:
        vcf="data/dbsnp/GCF_000001405.25.gz",
        tbi="data/dbsnp/GCF_000001405.25.gz.tbi",
        txt="data/GCF_000001405.25_GRCh37.p13_assembly_report.txt",
    output:
        vcf="data/dbsnp/GRCh37.dbSNP155.vcf.gz",
        tbi="data/dbsnp/GRCh37.dbSNP155.vcf.gz.tbi",
        txt=temp("data/GCF_000001405.25_GRCh37.p13_assembly_report.chroms"),
    shell:
        "grep -v '^#' {input.txt} | awk -v FS='\\t' '{{ if ($2==\"assembled-molecule\") {{ print $7, $1 }} else {{ print $7, $5 }} }}' | grep -vw 'na' > {output.txt} && "
        "bcftools annotate --rename-chrs {output.txt} -Oz -o {output.vcf} {input.vcf} && "
        "bcftools index --tbi {output.vcf}"
