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


rule convert_gwas_catalog_metadata:
    input:
        tsv="data/targets/gwas-association-accessionId_{accession}.tsv",
    output:
        tsv="data/targets/gwas_{accession}.tsv",
    shell:
        "Rscript scripts/convert_gwas_catalog.R"
        " --gwas {input.tsv}"
        " --output {output.tsv}"


rule convert_t1d_metadata:
    """
    Convert the GWAS metadata from (Robertson et al., 2021), for Type 1 Diabetes 

    https://doi.org/10.1038/s41588-021-00880-5
    """
    input:
        tsv="data/targets/Robertson_et_al_2021_ST7.tsv",
    output:
        tsv="data/targets/gwas_t1d.tsv",
    shell:
        "Rscript scripts/convert_t1d_metadata.R"
        " --gwas {input.tsv}"
        " --output {output.tsv}"


rule convert_t2d_metadata:
    """
    Convert the GWAS metadata from (Vujkovic et al., 2020), for Type 2 Diabetes 

    https://doi.org/10.1038/s41588-020-0637-y
    """
    input:
        tsv="data/targets/Vujkovic_et_al_2020_ST6.tsv",
    output:
        tsv="data/targets/gwas_t2d.tsv",
    shell:
        "Rscript scripts/convert_t2d_metadata.R"
        " --gwas {input.tsv}"
        " --output {output.tsv}"
