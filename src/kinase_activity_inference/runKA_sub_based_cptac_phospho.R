# Substrate-based kinase activity inference using CPTAC phosphorylation data


# load R packages
library(tidyverse)
library(parallel)

# source kinase-activity inference pipeline
source("./src/utils/infer_kinase_activity.R")


# load phospho data and respective flag
args <- commandArgs(trailingOnly = TRUE)
#args <- c("./output/files/phosphoproteomicsQ_Protreg_allSamp_withZtransf.txt.gz", "./output/files/", "phosphoQ_Protregout_allSamp_withZtransf", "FALSE")

if(length(args) < 4) stop(
  "Four command line arguments required:
  The first the path to the file containing the phosphorylation data.
  The second the output directory.
  The third an identifier string to append to the output file name.
  The fourth a logical to indicate whether the kinase auto-regulatory phosphosites should be removed from the kinase-substrate list.")


# load cptac phosphorylation data
cptac_pho <- read_tsv(file = args[1]) %>%
  filter(psites == 1) %>%
  select(-gene, -psites) %>%
  separate(col = "psite", into = c("gene", "psite"), sep = "_") %>%
  mutate(residue = str_extract(psite, "^[STY]{1}")) %>%
  mutate(position = as.numeric(str_extract(psite, "[0-9]+"))) %>%
  select(-psite) %>%
  select(gene, position, residue, everything())


# load kinase-substrate list
ks <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz")


# prepare kinase-substrate lists for inference
ks_lists <- ks %>%
  select(-pair, -source) %>%
  distinct() %>%
  mutate(source_type = str_replace(source_type, "database|text-mining", "DB_text-mining")) %>%
  distinct()


# remove auto-regulated phosphosites?
if(as.logical(args[4])){
  ks_lists <- ks_lists %>%
    filter(map2_lgl(.x=kinase, .y=substrate, .f = ~ .x != .y))
}


# infer kinase activity
# use z-test without weights
inference <- inferKA_W(
  ks_lists = ks_lists,
  phospho = cptac_pho,
  with_weights = F,
  multi_core = T)


# export inference data
if(!as.logical(args[4])){
  file_name <- paste0(args[2], "CPTAC_KA_", args[3], ".txt.gz")
} else {
  file_name <- paste0(args[2], "CPTAC_KA_", args[3], "_notKinSites.txt.gz")
}

gz1 <- gzfile(file_name, "w")
write.table(inference, gz1, sep="\t", quote=F, row.names=F)
close(gz1)

