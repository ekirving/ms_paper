#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import pandas as pd

# parameters for PLINK clumping
PLINK_CLUMP_PVAL = 5e-8
PLINK_CLUMP_R2 = 0.05
PLINK_CLUMP_KB = 250


def get_samples(config, dataset):
    """
    Get the samples in this dataset
    """
    return pd.read_table(config["samples"][dataset]["metadata"]).set_index("sampleId", drop=False)


def get_ancient_samples(config, wildcards):
    """
    Get list of sample IDs for all the ancient samples in the current analysis group
    """
    samples = get_samples(config, wildcards.dataset)
    return [str(s) for s in samples[samples["age"] != 0]["sampleId"].tolist()]


def get_modern_samples(config, wildcards):
    """
    Get list of sample IDs for all the modern samples in the current analysis group
    """
    samples = get_samples(config, wildcards.dataset)
    return [str(s) for s in samples[samples["age"] == 0]["sampleId"].tolist()]


def get_modern_pops(config, wildcards):
    """
    Get list of population IDs for all the modern samples in the current analysis group
    """
    samples = get_samples(config, wildcards.dataset)
    return sorted(set(samples[samples["age"] == 0]["popId"].tolist()))


def gen_dict_extract(key, var: dict):
    """Find keys in nested dictionaries"""
    if hasattr(var, "items"):
        for k, v in var.items():
            if k == key:
                yield v
            if isinstance(v, dict):
                for result in gen_dict_extract(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract(key, d):
                        yield result
    elif isinstance(var, list):
        for d in var:
            for result in gen_dict_extract(key, d):
                yield result
