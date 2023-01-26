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
@click.option("--output", metavar="<file>", type=click.Path(writable=True), help="Output filename", required=True)
def palm_report(models, output):
    """
    Make a PALM multi-trait report for all overlapping traits in UKBB and FinnGEN
    """

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

            rows.append({**model, **trait1, **trait2})

    # convert to a df
    df = pd.DataFrame(rows)
    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    palm_report()
