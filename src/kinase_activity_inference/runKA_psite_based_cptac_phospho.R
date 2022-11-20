# Regulatory phosphosite-based kinase activity inference using CPTAC phosphorylation data


# load R packages
library(tidyverse)
library(parallel)

# source kinase-activity inference pipeline
source("./src/utils/infer_kinase_activity2.R")


# load phospho data and respective flag
args <- commandArgs(trailingOnly = TRUE)
#args <- c("./output/files/phosphoproteomics.txt.gz", "./output/files/", "phospho")

if(length(args) < 3) stop(
  "Three command line arguments required:
  The first the path to the file containing the phosphorylation data.
  The second the output directory.
  The third an identifier string to append to the output file name.")


# load cptac phosphorylation data
cptac_pho <- read_tsv(file = args[1]) %>%
  filter(psites == 1) %>%
  select(-gene, -psites) %>%
  separate(col = "psite", into = c("gene", "psite"), sep = "_") %>%
  mutate(residue = str_extract(psite, "^[STY]{1}")) %>%
  mutate(position = as.numeric(str_extract(psite, "[0-9]+"))) %>%
  select(-psite) %>%
  select(gene, position, residue, everything())


# load functional score of kinases phosphosites from cptac
cptac_kin_funscoR <- read_tsv(file = "./output/files/cptac_kin_psites_funscoR.txt.gz") %>%
  filter(map2_lgl(.x = PSP_reg_status, .y = probabilities, .f = ~ if(.x == "known"){TRUE}else{if(.x == "unknown" & .y > 0.4){TRUE}else{FALSE}})) %>%
  select(-probabilities, -PSP_reg_status)


# infer kinase activity
inference <- infer_kinase_activity(phospho = cptac_pho, kin_list = cptac_kin_funscoR, multi_core = T)


# export inference data
file_name <- paste0(args[2], "CPTAC_KA_kin_psites_", args[3], ".txt.gz")
write_tsv(inference, file_name)

