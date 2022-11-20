# Kinase-activity inference using Claudia's benchmark phosphorylation data


library(tidyverse)
library(parallel)
library(viridis)
library(ROCR)



# source kinase-activity inference pipeline
source("./src/utils/infer_kinase_activity.R")

# source match pipeline
source("./src/utils/matchTool.R")


# load kinase-substrate list
ks <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz")


# just in vitro kinase-substrate pairs
psp_ks <- ks %>%
  filter(source == "psp") %>%
  select(kinase, substrate, position, residue) %>%
  distinct()


# load protein sequences for the canonical uniprot transcripts
prot_seqs <- read_tsv(file = "./output/files/uniprot_canTrpt_protSeqGRCh37.txt")


# build PWMs for each kinase
# select k-mers with 5 aa upstream and downstream the phosphosite
# remove k-mers with less than 11 aa (l = k*2+1)
# remove kinases with less than 5 kmers to build a PWM
kinase_pwm <- getKinasePWMs(psp_ks, prot_seqs, k=7, rm_kmer=T, n_kmers=10, rm_kinase=T)


# load benchmark phosphorylation data
bk_pho <- read_tsv(file = "./output/files/bk_pho.txt.gz")


# infer kinase activity with weights
inference1 <- infer_kinase_activity(
  ks_list = psp_ks,
  phospho = bk_pho,
  with_weights = T,
  prot_seqs = prot_seqs,
  kinase_pwm = kinase_pwm,
  multi_core = T)

# infer kinase activity without weights
inference2 <- infer_kinase_activity(
  ks_list = psp_ks %>% filter(kinase %in% kinase_pwm$kinase),
  phospho = bk_pho,
  with_weights = F,
  multi_core = T)



# ROC curves

# load known condition-specific kinase regulation
kreg <- read_tsv(file = "./data/kinase_substrate/benchmark_dataset/kinase_condition_pairs.txt") %>%
  select(sample=Condition,kinase=Kinase) %>%
  mutate(regulated = 1)


inference1_r <- compute_roc(inference1, kreg) %>% mutate(test = "z-test-weighted")
inference2_r <- compute_roc(inference2, kreg) %>% mutate(test = "z-test")


# prepare data for plotting
roc_data <- bind_rows(
  inference1_r,
  inference2_r) %>%
  mutate(auc = round(auc, 2)) %>%
  mutate(test = str_c(test, auc, sep = ": ")) %>%
  mutate(test = fct_rev(fct_reorder(test, auc, median)))


# ROC curves
roc_curve1 <- roc_data %>%
  ggplot(mapping=aes(x = fpr, y = tpr, color = test)) +
  geom_line() +
  theme_classic() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    strip.text=element_text(colour="black", size=13),
    strip.background=element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=13)) +
  scale_color_viridis(discrete = T) +
  labs(x = "FPR", y = "TPR", title = "")

roc_curve1

ggsave("bk_phospho_psp_roc.png", plot=roc_curve1, path="./output/plots/kinase_activity_inference/", height=4, width=6)
unlink("bk_phospho_psp_roc.png")
