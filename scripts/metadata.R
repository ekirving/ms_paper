#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

library(tidyverse) # v1.3.1

setwd("/Users/evan/Dropbox/Code/ms")

ms <- read_tsv("data/targets/all_clumped_annotated_ms.tsv", col_types = cols(.default = "c"))
ra <- read_tsv("data/targets/all_clumped_annotated_ra.tsv", col_types = cols(.default = "c"))

table_cols <- colnames(ms)

# convert the `ra` columns to the same format as `ms`
ra <- ra %>%
    rename(effect_allele=A1, other_allele=A2, OR=`OR(A1)`) %>%
    separate(`OR_95%CIup-OR_95%CIlow`, into=c("OR_upper", "OR_lower"), sep="-") %>%
    mutate(beta=log(as.numeric(OR)), se=NA) %>%
    select(all_of(table_cols))
    
write_tsv(ra, "data/targets/all_clumped_annotated_ra.tsv")
