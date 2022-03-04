#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2022, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"


configfile: "config.yaml"


include: "rules/gwascat.smk"
include: "rules/variants.smk"
include: "rules/clues.smk"


rule all:
    input:
        "",
