# Kinase-activity inference using Claudia's benchmark phosphorylation data


library(tidyverse)
library(parallel)


# source kinase-activity inference pipeline
source("./src/utils/infer_kinase_activity.R")

# source match pipeline
source("./src/utils/matchTool.R")


# load kinase-substrate list
ks <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz")


# just in vitro kinase-substrate pairs
inVT_ks <- ks %>%
  filter(source_type == "in_vitro") %>%
  select(kinase, substrate, position, residue)


# all kinase-substrate pairs from all lists
all_ks <- ks %>%
  select(kinase, substrate, position, residue) %>%
  distinct()


# load protein sequences for the canonical uniprot transcripts
prot_seqs <- read_tsv(file = "./output/files/uniprot_canTrpt_protSeqGRCh37.txt")


# build PWMs for each kinase
# select k-mers with 5 aa upstream and downstream the phosphosite
# remove k-mers with less than 11 aa (l = k*2+1)
# remove kinases with less than 5 kmers to build a PWM
kinase_pwm <- getKinasePWMs(all_ks, prot_seqs, k=7, rm_kmer=T, n_kmers=10, rm_kinase=T)


# load benchmark phosphorylation data
bk_pho <- read_tsv(file = "./output/files/bk_pho.txt.gz")


# multiply the log2fc by the phosphosite functional score

# load the functional score reported on David's paper
#fun_score <- read_tsv(file = "./output/files/fun_score_David_Paper.txt")

#bk_pho <- bk_pho %>%
#  pivot_longer(cols = -c(gene, position, residue), names_to = "sample", values_to = "log2fc") %>%
#  inner_join(fun_score, by = c("gene", "position")) %>%
#  mutate(log2fc = log2fc*functional_score) %>%
#  select(-uniprot, -functional_score) %>%
#  pivot_wider(names_from = "sample", values_from = "log2fc")


# prepare kinase-substrate lists for inference
ks_lists <- ks %>%
  select(-pair, -source) %>%
  distinct()

db <- ks_lists %>% filter(source_type == "database")
tm <- ks_lists %>% filter(source_type == "text-mining")
tmNoDB <- tm %>%
  anti_join(db, by = c("kinase", "substrate", "position", "residue")) %>%
  mutate(source_type = "text-mining-NoDB")

in_vv <- ks_lists %>% filter(source_type == "in_vivo") %>% select(1:2) %>% distinct()
in_vt <- ks_lists %>% filter(source_type == "in_vitro") %>% select(1:2) %>% distinct()
in_vivt <- dplyr::intersect(in_vv, in_vt)

in_vv <- ks_lists %>%
  filter(source_type == "in_vivo") %>%
  semi_join(in_vivt, by = c("kinase", "substrate")) %>%
  mutate(source_type = "in_vivo-KSin_VVVT")

in_vt <- ks_lists %>%
  filter(source_type == "in_vitro") %>%
  semi_join(in_vivt, by = c("kinase", "substrate")) %>%
  mutate(source_type = "in_vitro-KSin_VVVT")

ks_lists <- bind_rows(ks_lists, tmNoDB, in_vv, in_vt)



# infer kinase activity with weights
inference <- inferKA_W(
  ks_lists = ks_lists,
  phospho = bk_pho,
  with_weights = T,
  prot_seqs = prot_seqs,
  kinase_pwm = kinase_pwm,
  multi_core = T)

gz1 <- gzfile("./output/files/BKphospho_kinase_activity_ztestW.txt.gz", "w")
write.table(inference, gz1, sep="\t", quote=F, row.names=F)
close(gz1)



# infer kinase activity without weights
# use the same kinases used on the weighted test
inference <- inferKA_W(
  ks_lists = ks_lists %>% filter(kinase %in% kinase_pwm$kinase),
  phospho = bk_pho,
  with_weights = F,
  multi_core = T)

gz1 <- gzfile("./output/files/BKphospho_kinase_activity_ztestNW.txt.gz", "w")
write.table(inference, gz1, sep="\t", quote=F, row.names=F)
close(gz1)
