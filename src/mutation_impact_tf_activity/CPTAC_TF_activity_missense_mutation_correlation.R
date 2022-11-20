#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)

#' load R packages
library(tidyverse)

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load CPTAC samples
source("./src/utils/getSamples.R")
cptac_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")
cptac_samples <- getSamples(cptac_samples, c("protein","mRNA")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = str_c(batch, tissue, sep = "-")) %>%
  select(-tissue)


#' load CPTAC TF activities
tf_activity <- read_csv(file = "./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "tf_activity") %>%
  rename(tf=X1) %>%
  select(sample, tf, tf_activity)


#' load CPTAC mutations
cptac_mutations <- read_tsv(file = "./output/files/mutations_protpos.txt.gz")


#' load SIFT data
sift_map <- read_tsv(file = "./data/mutfunc/human/sift_uniprot_acc_genename.tab") %>%
  rename(acc = From, gene = To)

sift <- data.table::fread(file = "./data/mutfunc/human/sift.tab.gz") %>%
  as_tibble() %>%
  inner_join(sift_map, by = c("acc")) %>%
  select(acc, gene, everything())

exclude <- sift %>%
  select(gene, acc) %>%
  distinct() %>%
  group_by(gene) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  filter(n > 1)

sift <- sift %>%
  filter(!(gene %in% exclude$gene)) %>%
  filter(gene %in% unique(tf_activity$tf)) %>%
  select(gene, pos, ref, alt, score)


#' select only missense mutations\
#' select only the gene, sample, protein position, wt aa and mutated aa\
#' remove duplicates (same gene in the same sample can have multiple missense mutations that map to the same protein position and alteration)
missense <- cptac_mutations %>%
  filter(variant_class == "Missense_Mutation") %>%
  select(sample, gene=gene_symbol, prot_pos_i = prot_pos, aa_wt, aa_mut) %>%
  distinct()


#' # -- assess the impact of missense mutations on TF activities
missense_impact <- missense %>%
  inner_join(tf_activity, by = c("sample", "gene" = "tf")) %>%
  inner_join(sift, by = c("gene", "prot_pos_i"="pos","aa_wt"="ref","aa_mut"="alt")) %>%
  group_by(sample, gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  mutate(sift = -log10(score+1e-10)) %>%
  #filter(sift < 5) %>%
  ggplot(mapping = aes(x = sift, y = tf_activity)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(method = "pearson")

#+ fig.width=6, fig.height=6
missense_impact

ggsave(filename = "tf_act_missense_mutation_sift_correlation.png", plot = missense_impact, path = "./output/plots/mutation_impact_tf_activity/", width = 6, height = 6)
unlink("tf_act_missense_mutation_sift_correlation.png")

