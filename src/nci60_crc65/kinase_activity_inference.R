# Substrate-based kinase activity inference


# load R packages
library(tidyverse)
library(parallel)

# source kinase-activity inference pipeline
source("./src/utils/infer_kinase_activity.R")


# load paths to phosphorylation data and output file
#args <- commandArgs(trailingOnly = TRUE)
args <- c("./output/files/nci60_crc65_phospho_log2fc_ProtRegOut.txt.gz", "./output/files/nci60_crc65_kinase_activities.txt.gz")


# load phosphorylation data
nci60_phospho <- read_tsv(file = args[1]) %>%
  filter(batch == "NCI60") %>%
  select(-batch) %>%
  pivot_wider(names_from = "sample", values_from = "log2fc") %>%
  select(gene, position, residue, everything())

crc65_phospho <- read_tsv(file = args[1]) %>%
  filter(batch == "CRC65") %>%
  select(-batch) %>%
  pivot_wider(names_from = "sample", values_from = "log2fc") %>%
  select(gene, position, residue, everything())


# load kinase-substrate list
# prepare kinase-substrate lists for inference
ks_lists <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz") %>%
  select(-pair, -source) %>%
  distinct() %>%
  mutate(source_type = str_replace(source_type, "database|text-mining", "DB_text-mining")) %>%
  distinct()


# infer kinase activities
# use z-test without weights
nci60_KA <- inferKA_W(
  ks_lists = ks_lists,
  phospho = nci60_phospho,
  with_weights = F,
  multi_core = T)

nci60_KA <- nci60_KA %>%
  mutate(batch = "NCI60") %>%
  select(batch, sample, everything())

crc65_KA <- inferKA_W(
  ks_lists = ks_lists,
  phospho = crc65_phospho,
  with_weights = F,
  multi_core = T)

crc65_KA <- crc65_KA %>%
  mutate(batch = "CRC65") %>%
  select(batch, sample, everything())

kin_activity <- bind_rows(nci60_KA, crc65_KA)

write_tsv(kin_activity, args[2])
