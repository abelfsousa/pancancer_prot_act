# get TF phosphosites with respective functional score

library(tidyverse)


# load TF activities and get TF IDs
tfs <- data.table::fread("./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  as_tibble() %>%
  rename(tf=V1) %>%
  pivot_longer(-tf, names_to = "sample", values_to = "tf_activity") %>%
  pull(tf) %>%
  unique()


# load CPTAC phosphosites with respective functional score
cptac_psites_scores <- read_tsv("./output/files/cptac_psites_funscoR.txt.gz")


# select phosphosites in TFs
tf_psites <- cptac_psites_scores %>%
  filter(gene %in% tfs) %>%
  select(-acc)

write_tsv(tf_psites, "./output/files/cptac_tf_psites_funscoR.txt.gz")

