# Kinase-activity inference using Claudia's benchmark phosphorylation data
# Use the same kinase-substrate lists that were used in the CPTAC KA


library(tidyverse)
library(parallel)


# source kinase-activity inference pipeline
source("./src/utils/infer_kinase_activity.R")


# load kinase-substrate list
ks <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz")


# prepare kinase-substrate lists for inference
ks_lists <- ks %>%
  select(-pair, -source) %>%
  distinct() %>%
  mutate(source_type = str_replace(source_type, "database|text-mining", "DB_text-mining")) %>%
  distinct()


# load benchmark phosphorylation data
bk_pho <- read_tsv(file = "./output/files/bk_pho.txt.gz")



# infer kinase activity
# use z-test without weights
inference <- inferKA_W(
  ks_lists = ks_lists,
  phospho = bk_pho,
  with_weights = F,
  multi_core = F)

write_tsv(inference, "./output/files/BK_phospho_KA.txt.gz")
