#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

# FTP server for the GWAS Catalog
GWASCAT_FTP = "ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases"

# https://www.nature.com/articles/ejhg2015269
GWASCAT_PVALUE = "5e-8"

import gzip


rule gwascat_fetch:
    output:
        "data/gwascat/gwas_catalog_v1.0.2-associations_e105_r2022-02-21.tsv.gz",
    shell:
        "wget --quiet -O - {GWASCAT_FTP}/2022/02/21/gwas-catalog-associations_ontology-annotated.tsv | gzip > {output}"


rule gwascat_genome_wide_significant:
    input:
        "data/gwascat/gwas_catalog_v1.0.2-associations_e105_r2022-02-21.tsv.gz",
    output:
        "data/gwascat/gwas_catalog_significant.tsv.gz",
    params:
        col=lambda wildcards, input: gzip.open(input[0], "rb").readline().split(b"\t").index(b"P-VALUE") + 1,
    shell:
        r"gunzip -c {input} | awk -F'\t' 'NR==1 || ${params.col} < {GWASCAT_PVALUE} {{print $0}}' | gzip > {output}"
