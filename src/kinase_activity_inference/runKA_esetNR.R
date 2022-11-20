# Kinase-activity inference using David/Danish compilation of phosphorylation data
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


# load phosphorylation data
phospho <- data.table::fread(file = "./output/files/esetNR_phospho.txt.gz") %>%
  as_tibble()


# infer kinase activity
# use z-test without weights
inference <- inferKA_W(
  ks_lists = ks_lists,
  phospho = phospho,
  with_weights = F,
  multi_core = T)

write_tsv(inference, "./output/files/KA_esetNR.txt.gz")
