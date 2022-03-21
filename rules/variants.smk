#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

# maximum concurrent requests
ENSEMBL_MAX_CONCURRENT = 15


rule variant_ensembl:
    output:
        var="data/ensembl/{reference}/var/{rsid}.json",
        vep="data/ensembl/{reference}/vep/{rsid}.json",
    resources:
        ensembl_api=1,
    shell:
        "python scripts/variant_ensembl.py"
        " --ref {wildcards.reference}"
        " --rsid {wildcards.rsid}"
        " --var {output.var}"
        " --vep {output.vep}"


rule variant_metadata:
    input:
        var="data/ensembl/{reference}/var/{rsid}.json",
        vep="data/ensembl/{reference}/vep/{rsid}.json",
        gwas="data/gwascat/gwas_catalog_significant.tsv.gz",
    output:
        "data/metadata/{reference}/{rsid}.json",
    resources:
        mem_mb=8000,
    wildcard_constraints:
        rsid="rs\d+",
    shell:
        "python scripts/variant_metadata.py"
        " --rsid {wildcards.rsid}"
        " --var {input.var}"
        " --vep {input.vep}"
        " --gwas {input.gwas}"
        " --out {output}"


rule variant_label:
    input:
        meta=lambda wildcards: "data/metadata/{reference}/{{rsid}}.json".format(
            reference=config["samples"][wildcards.dataset]["reference"]
        ),
    output:
        "data/labels/{dataset}/{dataset}-{rsid}-label.json",
    params:
        vcf=lambda wildcards: config["samples"][wildcards.dataset]["genotypes"],
    shell:
        "python scripts/variant_label.py"
        " --vcf {params.vcf}"
        " --meta {input.meta}"
        " --out {output}"
