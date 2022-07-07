#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

# FTP server for the GWAS Catalog
GWASCAT_FTP = "ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases"

# https://www.nature.com/articles/ejhg2015269
GWASCAT_PVALUE = "5e-8"

import gzip


rule gwascat_fetch:
    """
    Download a recent build of the GWAS Catalog
    """
    output:
        "data/gwascat/gwas_catalog_v1.0.2-associations_e105_r2022-02-21.tsv.gz",
    shell:
        "wget --quiet -O - {GWASCAT_FTP}/2022/02/21/gwas-catalog-associations_ontology-annotated.tsv | gzip > {output}"


rule gwascat_genome_wide_significant:
    """
    Filter the GWAS Catalog to only retain genome-wide significant associations (i.e., p<5e-8)
    """
    input:
        "data/gwascat/gwas_catalog_v1.0.2-associations_e105_r2022-02-21.tsv.gz",
    output:
        "data/gwascat/gwas_catalog_significant.tsv.gz",
    params:
        col=lambda wildcards, input: gzip.open(input[0], "rb").readline().split(b"\t").index(b"P-VALUE") + 1,
    shell:
        r"gunzip -c {input} | awk -F'\t' 'NR==1 || ${params.col} < {GWASCAT_PVALUE} {{print $0}}' | gzip > {output}"


rule convert_ms_metadata:
    """
    Convert the GWAS metadata from (International Multiple Sclerosis Genetics Consortium, 2019), for Multiple sclerosis 

    https://doi.org/10.1126/science.aav7188
    """
    input:
        tsv="data/targets/ms_snps_final_discovery_0.7_combined.csv",
    output:
        tsv="data/targets/gwas_ms.tsv",
    shell:
        "Rscript scripts/convert_metadata.R"
        " --gwas {input.tsv}"
        " --output {output.tsv}"


rule convert_cd_metadata:
    """
    Convert the GWAS metadata for Celiac disease 
    """
    input:
        tsv="data/targets/cd_snps_final_0.7_combined.csv",
    output:
        tsv="data/targets/gwas_cd.tsv",
    shell:
        "Rscript scripts/convert_metadata.R"
        " --gwas {input.tsv}"
        " --output {output.tsv}"


rule convert_ra_metadata:
    """
    Convert the GWAS metadata for  
    """
    input:
        tsv="data/targets/ra_snps_final_0.7_combined.csv",
    output:
        tsv="data/targets/gwas_ra.tsv",
    shell:
        "Rscript scripts/convert_metadata.R"
        " --gwas {input.tsv}"
        " --output {output.tsv}"


# rule convert_ms_patsopoulos_metadata:
#     """
#     Convert the four different SNP ascertainments from Patsopoulos et. al. 2019
#
#     1. Genome-wide significant autosomal SNPs (n=200)
#     2. MHC SNPs, classical model (n=32)
#     3. Strongly suggestive SNPs [p < 1e-5] (n=117)
#     4. Weakly suggestive SNPs (n=299)
#
#     https://doi.org/10.1126/science.aav7188
#     """
#     input:
#         auto="data/targets/Patsopoulos_et_al_2019_ST7.tsv",
#         mhc="data/targets/Patsopoulos_et_al_2019_ST11.tsv",
#         sS="data/targets/Patsopoulos_et_al_2019_ST14_sS.tsv",
#         wS="data/targets/Patsopoulos_et_al_2019_ST14_wS.tsv",
#     output:
#         auto="data/targets/gwas_ms-auto.tsv",
#         mhc="data/targets/gwas_ms-mhc.tsv",
#         sig="data/targets/gwas_ms-auto-mhc.tsv",
#         sug="data/targets/gwas_ms-auto-mhc-sS.tsv",
#         all="data/targets/gwas_ms-auto-mhc-sS-wS.tsv",
#     shell:
#         "Rscript scripts/convert_Patsopoulos_metadata.R"
#         " --auto {input.auto}"
#         " --mhc {input.mhc}"
#         " --sS {input.sS}"
#         " --wS {input.wS}"
#         " --out-auto {output.auto}"
#         " --out-mhc {output.mhc}"
#         " --out-sig {output.sig}"
#         " --out-sug {output.sug}"
#         " --out-all {output.all}"


rule convert_ms_shams_metadata:
    """
    Convert the metadata from Shams et al., 2022

    https://doi.org/10.1093/brain/awac092
    """
    input:
        tsv="data/targets/Shams_et_al_2022_S15.tsv",
    output:
        tsv="data/targets/gwas_ms-mr.tsv",
    shell:
        "Rscript scripts/convert_Shams_metadata.R"
        " --gwas {input.tsv}"
        " --output {output.tsv}"


rule gwas_metadata_mhc:
    """
    Split the MHC and non-MHC SNPs into two files, so we can compare results  
    """
    input:
        gwas="data/targets/gwas_{trait}.tsv",
    output:
        mhc="data/targets/gwas_{trait}-mhc.tsv",
        auto="data/targets/gwas_{trait}-auto.tsv",
    shell:
        "awk 'NR==1 || ($1=={MHC_CHROM} && $2>={MHC_START} && $2<={MHC_FINISH})' {input.gwas} > {output.mhc} && "
        "awk 'NR==1 || ($1!={MHC_CHROM} || $2< {MHC_START} || $2> {MHC_FINISH})' {input.gwas} > {output.auto}"


checkpoint gwas_metadata:
    """
    Wrap the different sources of GWAS metadata, so we can access the list of SNPs as a checkpoint 
    """
    input:
        tsv="data/targets/gwas_{trait}.tsv",
    output:
        tsv=temp("data/targets/gwas_{trait}_static.tsv"),
    shell:
        "cp {input} {output}"
