#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import json

import click
from scipy.stats import norm


@click.command()
@click.option("--dataset", metavar="<string>", help="Dataset", required=True)
@click.option("--ancestry", metavar="<string>", help="Ancestral path", required=True)
@click.option("--trait1", metavar="<string>", help="GWAS trait", required=True)
@click.option("--trait2", metavar="<string>", help="GWAS trait", required=True)
@click.option(
    "--palm", "palm_file", metavar="<file>", type=click.Path(readable=True), help="PALM output file", required=True
)
@click.option("--out", "output", metavar="<file>", type=click.File("w"), help="Output filename", required=True)
def palm_parse_txt(dataset, ancestry, trait1, trait2, palm_file, output):
    """
    Parse the PALM log file to extract the information we want.
    """
    with open(palm_file) as fin:
        result = False

        results = {}

        for line in fin:
            words = line.split()
            num_words = len(words)
            if num_words == 3 and words[0] == "Analyzing" and words[2] == "loci...":
                num_loci = words[1]
            elif num_words == 1 and words[0][0] == "=":
                result = True
            elif result and len(words) == 6:
                trait, sel, se, z, zmarg, r = words

                # the PALM output truncates long trait names
                trait = trait1 if trait1.startswith(trait) else trait2

                results[trait] = {
                    "sel": sel,
                    "se": se[1:-1],
                    "z": z,
                    "zmarg": zmarg,
                    "r": r,
                    "pjoint": "{:.2E}".format(norm.sf(abs(float(z))) * 2),
                    "pmarg": "{:.2E}".format(norm.sf(abs(float(zmarg))) * 2),
                    "pR": "{:.2E}".format(norm.sf(abs(float(r))) * 2),
                }

        data = {
            "dataset": dataset,
            "ancestry": ancestry,
            "trait1": trait1,
            "trait2": trait2,
            "num_loci": num_loci,
            "results": results,
        }

    json.dump(data, output, indent=2)


if __name__ == "__main__":
    palm_parse_txt()
