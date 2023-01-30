#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import json

import click
import pandas as pd


@click.command()
@click.option("--model", "models", metavar="<file>", type=click.Path(exists=True), multiple=True, required=True)
@click.option("--ukbb", "ukbb_tsv", metavar="<file>", type=click.Path(exists=True), required=True)
@click.option("--finngen", "finngen_tsv", metavar="<file>", type=click.Path(exists=True), required=True)
@click.option("--output", metavar="<file>", type=click.Path(writable=True), help="Output filename", required=True)
def palm_report(models, ukbb_tsv, finngen_tsv, output):
    """
    Make a PALM multi-trait report for all overlapping traits in UKBB and FinnGEN
    """

    ukbb = pd.read_table(ukbb_tsv, compression="gzip").set_index("phenotype", drop=False)
    finngen = pd.read_table(finngen_tsv).set_index("phenocode", drop=False)

    rows = []

    for model_json in models:

        # load the SNP metadata
        with open(model_json) as fin:
            model = json.load(fin)

            # get the results for the two traits
            results = model.pop("results")

            trait1 = results[model["trait1"]]
            trait1 = dict(zip([f"trait1_{key}" for key in trait1.keys()], trait1.values()))

            trait2 = results[model["trait2"]]
            trait2 = dict(zip([f"trait2_{key}" for key in trait2.keys()], trait2.values()))

            model["trait2"], model["trait2_gwas"] = model["trait2"].split("-")

            if model["trait2_gwas"] == "ukbb":
                model["trait2_pheno"] = ukbb.loc[model["trait2"]]["description"]
            else:
                model["trait2_pheno"] = finngen.loc[model["trait2"]]["name"]

            rows.append({**model, **trait1, **trait2})

    # convert to a df
    df = pd.DataFrame(rows)
    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    palm_report()
