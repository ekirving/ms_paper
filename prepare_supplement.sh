#!/usr/bin/env bash

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

mkdir -p supplement/

# ----------------------------------------------------------------------------------------------------------------------
# prepare all the supplementary figures
# ----------------------------------------------------------------------------------------------------------------------

# MS / PALM analyses
cp results/palm/ancestral_paths_v3-ALL-ms-r0.05-kb250-palm_lines-pval.png supplement/fig_S6.1_ancestral_paths_v3-ALL-ms-r0.05-kb250-palm_lines-pval.png
cp results/palm/ancestral_paths_v3-WHG-ms-r0.05-kb250-palm_lines-pval.png supplement/fig_S6.2_ancestral_paths_v3-WHG-ms-r0.05-kb250-palm_lines-pval.png
cp results/palm/ancestral_paths_v3-EHG-ms-r0.05-kb250-palm_lines-pval.png supplement/fig_S6.3_ancestral_paths_v3-EHG-ms-r0.05-kb250-palm_lines-pval.png
cp results/palm/ancestral_paths_v3-CHG-ms-r0.05-kb250-palm_lines-pval.png supplement/fig_S6.4_ancestral_paths_v3-CHG-ms-r0.05-kb250-palm_lines-pval.png
cp results/palm/ancestral_paths_v3-ANA-ms-r0.05-kb250-palm_lines-pval.png supplement/fig_S6.5_ancestral_paths_v3-ANA-ms-r0.05-kb250-palm_lines-pval.png

# MS / Cross ancestry comparisons
cp results/palm/ancestral_paths_v3-ms-r0.05-kb250-delta_prs.png supplement/fig_S6.6_ancestral_paths_v3-ms-r0.05-kb250-delta_prs.png
cp results/palm/ancestral_paths_v3-ms-r0.05-kb250-scatter.png   supplement/fig_S6.7_ancestral_paths_v3-ms-r0.05-kb250-scatter.png

# MS / Pleiotropic UK Biobank traits
cp results/compare/ancestral_paths_v3-ms-r0.05-kb250-ukbb-marginal-001.png supplement/fig_S6.8_ancestral_paths_v3-ms-r0.05-kb250-ukbb-marginal-001.png
cp results/compare/ancestral_paths_v3-ms-r0.05-kb250-ukbb-marginal-002.png supplement/fig_S6.9_ancestral_paths_v3-ms-r0.05-kb250-ukbb-marginal-002.png
cp results/compare/ancestral_paths_v3-ms-r0.05-kb250-ukbb-marginal-003.png supplement/fig_S6.10_ancestral_paths_v3-ms-r0.05-kb250-ukbb-marginal-003.png
cp results/compare/ancestral_paths_v3-ms-r0.05-kb250-ukbb-upset-top.png    supplement/fig_S6.11_ancestral_paths_v3-ms-r0.05-kb250-ukbb-upset-top.png

# MS / Pleiotropic FinnGen traits
cp results/compare/ancestral_paths_v3-ms-r0.05-kb250-finngen-marginal-001.png     supplement/fig_S6.12_ancestral_paths_v3-ms-r0.05-kb250-finngen-marginal-001.png
cp results/compare/ancestral_paths_v3-ms-r0.05-kb250-finngen-marginal-002.png     supplement/fig_S6.13_ancestral_paths_v3-ms-r0.05-kb250-finngen-marginal-002.png
cp results/compare/ancestral_paths_v3-ms-r0.05-kb250-finngen-marginal-003.png     supplement/fig_S6.14_ancestral_paths_v3-ms-r0.05-kb250-finngen-marginal-003.png
cp results/compare/ancestral_paths_v3-ms-r0.05-kb250-finngen-upset-top.png        supplement/fig_S6.15_ancestral_paths_v3-ms-r0.05-kb250-finngen-upset-top.png
cp results/compare/ancestral_paths_v3-ms-r0.05-kb250-finngen-upset-infectious.png supplement/fig_S6.16_ancestral_paths_v3-ms-r0.05-kb250-finngen-upset-infectious.png

# MS / Joint polygenic selection analysis
cp results/palm/ancestral_paths_v3-ms-r0.05-kb250-palm_report_multi_trait-ukbb.png    supplement/fig_S6.17_ancestral_paths_v3-ms-r0.05-kb250-palm_report_multi_trait-ukbb.png
cp results/palm/ancestral_paths_v3-ms-r0.05-kb250-palm_report_multi_trait-finngen.png supplement/fig_S6.18_ancestral_paths_v3-ms-r0.05-kb250-palm_report_multi_trait-finngen.png

#  ------------

# RA / PALM analyses
cp results/palm/ancestral_paths_v3-ALL-ra-r0.05-kb250-palm_lines-pval.png supplement/fig_S6.19_ancestral_paths_v3-ALL-ra-r0.05-kb250-palm_lines-pval.png
cp results/palm/ancestral_paths_v3-WHG-ra-r0.05-kb250-palm_lines-pval.png supplement/fig_S6.20_ancestral_paths_v3-WHG-ra-r0.05-kb250-palm_lines-pval.png
cp results/palm/ancestral_paths_v3-EHG-ra-r0.05-kb250-palm_lines-pval.png supplement/fig_S6.21_ancestral_paths_v3-EHG-ra-r0.05-kb250-palm_lines-pval.png
cp results/palm/ancestral_paths_v3-CHG-ra-r0.05-kb250-palm_lines-pval.png supplement/fig_S6.22_ancestral_paths_v3-CHG-ra-r0.05-kb250-palm_lines-pval.png
cp results/palm/ancestral_paths_v3-ANA-ra-r0.05-kb250-palm_lines-pval.png supplement/fig_S6.23_ancestral_paths_v3-ANA-ra-r0.05-kb250-palm_lines-pval.png

# RA / Cross ancestry comparisons
cp results/palm/ancestral_paths_v3-ra-r0.05-kb250-delta_prs.png supplement/fig_S6.24_ancestral_paths_v3-ra-r0.05-kb250-delta_prs.png
cp results/palm/ancestral_paths_v3-ra-r0.05-kb250-scatter.png   supplement/fig_S6.25_ancestral_paths_v3-ra-r0.05-kb250-scatter.png

# RA / Pleiotropic UK Biobank traits
cp results/compare/ancestral_paths_v3-ra-r0.05-kb250-ukbb-marginal-001.png supplement/fig_S6.26_ancestral_paths_v3-ra-r0.05-kb250-ukbb-marginal-001.png
cp results/compare/ancestral_paths_v3-ra-r0.05-kb250-ukbb-marginal-002.png supplement/fig_S6.27_ancestral_paths_v3-ra-r0.05-kb250-ukbb-marginal-002.png
cp results/compare/ancestral_paths_v3-ra-r0.05-kb250-ukbb-marginal-003.png supplement/fig_S6.28_ancestral_paths_v3-ra-r0.05-kb250-ukbb-marginal-003.png
cp results/compare/ancestral_paths_v3-ra-r0.05-kb250-ukbb-upset-top.png    supplement/fig_S6.29_ancestral_paths_v3-ra-r0.05-kb250-ukbb-upset-top.png

# RA / Pleiotropic FinnGen traits
cp results/compare/ancestral_paths_v3-ra-r0.05-kb250-finngen-marginal-001.png     supplement/fig_S6.30_ancestral_paths_v3-ra-r0.05-kb250-finngen-marginal-001.png
cp results/compare/ancestral_paths_v3-ra-r0.05-kb250-finngen-marginal-002.png     supplement/fig_S6.31_ancestral_paths_v3-ra-r0.05-kb250-finngen-marginal-002.png
cp results/compare/ancestral_paths_v3-ra-r0.05-kb250-finngen-marginal-003.png     supplement/fig_S6.32_ancestral_paths_v3-ra-r0.05-kb250-finngen-marginal-003.png
cp results/compare/ancestral_paths_v3-ra-r0.05-kb250-finngen-upset-top.png        supplement/fig_S6.33_ancestral_paths_v3-ra-r0.05-kb250-finngen-upset-top.png
cp results/compare/ancestral_paths_v3-ra-r0.05-kb250-finngen-upset-infectious.png supplement/fig_S6.34_ancestral_paths_v3-ra-r0.05-kb250-finngen-upset-infectious.png

# ----------------------------------------------------------------------------------------------------------------------
# prepare all the supplementary tables
# ----------------------------------------------------------------------------------------------------------------------

# MS `all` and `pruned` tables
cp results/palm/ancestral_paths_v3-ms-all-palm_report_prs.tsv         supplement/table_6.i_ancestral_paths_v3-ms-all-palm_report_prs.tsv
cp results/palm/ancestral_paths_v3-ms-r0.05-kb250-palm_report_prs.tsv supplement/table_6.ii_ancestral_paths_v3-ms-r0.05-kb250-palm_report_prs.tsv

# RA `all` and `pruned` tables
cp results/palm/ancestral_paths_v3-ra-all-palm_report_prs.tsv         supplement/table_6.iii_ancestral_paths_v3-ra-all-palm_report_prs.tsv
cp results/palm/ancestral_paths_v3-ra-r0.05-kb250-palm_report_prs.tsv supplement/table_6.iv_ancestral_paths_v3-ra-r0.05-kb250-palm_report_prs.tsv

# MS joint-PALM
cp results/palm/ancestral_paths_v3-ms-r0.05-kb250-palm_report_multi_trait.tsv supplement/table_6.v_ancestral_paths_v3-ms-r0.05-kb250-palm_report_multi_trait.tsv

