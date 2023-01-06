#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import json
from collections import defaultdict

import click
import pandas as pd
from scipy.stats.distributions import chi2

ANCESTRIES = ["ALL", "ANA", "CHG", "WHG", "EHG"]


@click.command()
@click.option("--data", "data_tsv", metavar="<file>", help="SNP data", type=click.Path(exists=True), required=True)
@click.option("--dataset", metavar="<string>", help="Name of the dataset", required=True)
@click.option("--ancestry", metavar="<string>", help="Ancestral path", type=click.Choice(ANCESTRIES), required=True)
@click.option("--output", metavar="<file>", type=click.Path(writable=True), help="Output filename", required=True)
def palm_report(data_tsv, dataset, ancestry, output):
    """
    Generate a PALM report
    """
    # get the list of SNPs to load
    data = pd.read_table(data_tsv)

    rows = []

    for _, snp in data.iterrows():
        rsid = snp["rsid"]

        # load the SNP metadata
        with open(f"data/metadata/GRCh37/{rsid}.json") as fin:
            meta = json.load(fin)

            gwascat = defaultdict(list)

            # extract the phenotypes
            for gwas in meta.get("gwascat", []):
                if gwas["trait"] != "":
                    # use the standardised trait ontology name
                    for trait in gwas["trait"].split(", "):
                        gwascat[trait].append(gwas["pubmedid"])
                else:
                    # use the free-text trait name reported in the submission
                    gwascat[gwas["phenotype"]].append(gwas["pubmedid"])

            snp["genes"] = meta["genes"]
            snp["gwascat"] = "; ".join(
                sorted([f"{trait} (PMID: {', '.join(studies)})" for trait, studies in gwascat.items()])
            )

        # compose the variant identifier
        variant = f"{snp['chrom']}:{snp['pos']}:{snp['ancestral_allele']}:{snp['derived_allele']}"

        # load the model parameters
        mod_file = f"results/clues/{rsid}/{dataset}-{variant}-{ancestry}.json"

        with open(mod_file) as fin:
            model = json.load(fin)

        # extract the first (and only) epoch
        epoch, s = list(model.pop("epochs").items())[0]

        model["epoch"] = epoch
        model["s"] = s

        rows.append({**snp, **model})

    # get the Bonferroni corrected significance threshold
    num_snps = data.shape[0]
    bonferroni = 0.05 / num_snps

    # convert to a df
    df = pd.DataFrame(rows)
    df["chrom"] = df.chrom.astype("int64")
    # convert the log-likelihood ratio into a p-value
    # https://en.wikipedia.org/wiki/Wilks%27_theorem
    df["p.value"] = df["logLR"].apply(lambda logLR: chi2.sf(2 * logLR, 1))
    df["significant"] = df["p.value"].apply(lambda p: p <= bonferroni)
    df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    palm_report()
