#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import json
import os
import re
import sys
from collections import OrderedDict

import click
import pandas as pd

sys.path.append(os.getcwd())

from scripts.utils import gen_dict_extract


@click.command()
@click.option("--rsid", metavar="<string>", help="RefSeq ID", required=True)
@click.option("--var", "var_file", metavar="<file>", type=click.File("r"), help="Ensembl VAR record", required=True)
@click.option("--vep", "vep_file", metavar="<file>", type=click.File("r"), help="Ensembl VEP record", required=True)
@click.option(
    "--gwas", "gwas_file", metavar="<file>", type=click.Path(exists=True), help="GWAS Catalog file", required=True
)
@click.option("--out", "output_file", metavar="<file>", type=click.File("w"), help="Output file", required=True)
def variant_metadata(rsid, var_file, vep_file, gwas_file, output_file):
    """
    Parse the Ensembl record to retrieve the rsID metadata
    """
    try:
        var = json.load(var_file)
        vep = json.load(vep_file)
    except json.decoder.JSONDecodeError as e:
        print("ERROR loading {} and {}".format(var_file, vep_file))
        raise e

    # get the clinical significance data from the VEP record
    clin_sig = set()
    clin_sig_allele = set()
    severity = set()
    gwas_records = []

    for val in gen_dict_extract("clin_sig_allele", vep):
        for assoc in val.split(";"):
            allele, significance = assoc.split(":")
            clin_sig.add(significance)
            clin_sig_allele.add(allele)

    for val in gen_dict_extract("most_severe_consequence", vep):
        severity.add(val)

    try:
        # get all the GWAS Catalog records for the target rsID
        gwas_target = pd.read_table(gwas_file, dtype=str).set_index("SNPS", drop=False).fillna("").loc[[rsid]]
    except KeyError:
        # no matching rows (i.e. this is a neutral SNP)
        gwas_ra = None
        gwas_genes = ""

    else:
        # handle multiple GWAS associations for this rsID
        for idx, gwascat in gwas_target.iterrows():
            # extract the GWAS Catalog allele
            gwascat_allele = gwascat["STRONGEST SNP-RISK ALLELE"].split("-").pop()

            # drop any GWAS studies where the RA is not listed
            if gwascat_allele not in ["A", "C", "G", "T", "?"]:
                continue

            # deduplicate within studies (due to duplicates in the GWAS Catalog), but keep the original order
            genes = re.split("[,;] *", gwascat["MAPPED_GENE"])
            genes = list(OrderedDict.fromkeys(genes))

            gwas_records.append(
                {
                    "pubmedid": gwascat["PUBMEDID"],
                    "genes": genes,
                    "allele": gwascat_allele,
                    "raf": gwascat["RISK ALLELE FREQUENCY"],
                    "p": gwascat["P-VALUE"],
                    "ororbeta": gwascat["OR or BETA"],
                    "phenotype": gwascat["DISEASE/TRAIT"],
                    "samples": gwascat["INITIAL SAMPLE SIZE"],
                    "trait": gwascat["MAPPED_TRAIT"],
                    "trait_uri": gwascat["MAPPED_TRAIT_URI"],
                }
            )

        # get the most common risk allele (some rsIDs have multiple)
        gwas_ra = gwas_target["STRONGEST SNP-RISK ALLELE"].value_counts().index[0].split("-").pop().replace("?", "")

        # deduplicate across studies, but keep the original order
        gwas_genes = list(OrderedDict.fromkeys([gene for record in gwas_records for gene in record["genes"]]))
        gwas_genes = ", ".join(sorted(gwas_genes))

    # always use the first mapping
    mappings = var.get("mappings", [])
    mapping = mappings[0] if len(mappings) > 0 else {}

    # get the ancestral and derived alleles
    ancestral = mapping.get("ancestral_allele")

    # use the Ensembl annotation
    alleles = mapping.get("allele_string", "").split("/")

    if ancestral is None or len(alleles) > 2:
        derived = None
    else:
        # the derived is whichever allele is left
        derived = (set(alleles) - {ancestral}).pop()

    # save the metadata for the variant
    metadata = {
        "rsid": var.get("name"),
        "type": var.get("var_class", var.get("error")),
        "chrom": mapping.get("seq_region_name"),
        "start": mapping.get("start"),
        "end": mapping.get("end"),
        "allele": gwas_ra,
        "genes": gwas_genes,
        "alleles": mapping.get("allele_string"),
        "ancestral": ancestral,
        "derived": derived,
        "minor": var.get("minor_allele"),
        "clin_allele": "/".join(set(clin_sig_allele)),
        "clin_sig": "/".join(set(clin_sig)),
        "severity": "/".join(set(severity)),
        "gwascat": gwas_records,
    }

    # save the metadata
    json.dump(metadata, output_file, indent=2)


if __name__ == "__main__":
    variant_metadata()
