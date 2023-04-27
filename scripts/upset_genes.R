library(tidyverse)

cd_finngen <- read_tsv("results/compare/ancestral_paths_v3-cd-r0.05-kb250-finngen-upset-top.tsv")
ms_finngen <- read_tsv("results/compare/ancestral_paths_v3-ms-r0.05-kb250-finngen-upset-top.tsv")
ra_finngen <- read_tsv("results/compare/ancestral_paths_v3-ra-r0.05-kb250-finngen-upset-top.tsv")

data_finngen <- bind_rows(
    cd_finngen %>% mutate(trait="CD"),
    ms_finngen %>% mutate(trait="MS"),
    ra_finngen %>% mutate(trait="RA")
) %>%
    group_by(nearest_genes) %>%
    summarize(traits=paste0(sort(unique(trait)), collapse = "-")) %>%
    group_by(traits) %>%
    summarize(genes=paste0(sort(unique(nearest_genes)), collapse = "; "))

write_tsv(data_finngen, "astrid/ancestral_paths_v3-finngen-upset-top.tsv")
