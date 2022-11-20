# build a position weighted matrix (PWM) for each kinase


library(tidyverse)

source("./src/utils/matchTool.R")


# load kinase-substrate lists
ks_lists <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz")

inVT_ks <- ks_lists %>%
  filter(source_type == "in_vitro") %>%
  select(kinase, substrate, position, residue)

all_ks <- ks_lists %>%
  select(kinase, substrate, position, residue) %>%
  distinct()


# load protein sequences for the canonical uniprot transcripts
gene_prot_seq <- read_tsv(file = "./output/files/uniprot_canTrpt_protSeqGRCh37.txt")


# build PWMs for each kinase
# select k-mers with 5 aa upstream and downstream the phosphosite
# remove k-mers with less than 11 aa (l = k*2+1)
# remove kinases with less than 5 kmers to build a PWM
kinase_pwm <- getKinasePWMs(all_ks, gene_prot_seq, k=5, rm_kmer=T, n_kmers=5, rm_kinase=T)

write_rds(kinase_pwm, "./output/r_objects/kinase_PWMs_allLists.rds", compress = "gz")
