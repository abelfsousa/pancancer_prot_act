#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)

#' load R packages
library(tidyverse)

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load CPTAC/CCLE samples
source("./src/utils/getSamples.R")
cptac_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")
cptac_samples <- getSamples(cptac_samples, c("protein","mRNA"))


#' load NCI60/CRC65 cell lines
nci60_crc65 <- read_tsv(file = "./output/files/nci60_crc65_cell_lines.txt")

shared_cell_lines <- nci60_crc65 %>%
  group_by(batch) %>%
  summarise(cell_line = list(cell_line)) %>%
  pull(cell_line) %>%
  reduce(intersect)

alternative_names <- nci60_crc65 %>%
  filter(!is.na(alternative_names)) %>%
  mutate(alternative_names = str_split(alternative_names, ";")) %>%
  unnest() %>%
  distinct()


#' load TF activity inference data from NCI60/CRC65 dataset
crc65_TF <- read_csv(file = "./data/dorothea/nci60_crc65/TF_activity_log2FC_allsamples_CRC65.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "tf_activity") %>%
  rename(tf=X1) %>%
  filter(!sample %in% shared_cell_lines) %>%
  select(sample, tf, tf_activity)

nci60_TF <- read_csv(file = "./data/dorothea/nci60_crc65/TF_activity_log2FC_allsamples_NCI60.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "tf_activity") %>%
  rename(tf=X1) %>%
  select(sample, tf, tf_activity)

nci60_crc65_TF <- bind_rows(crc65_TF, nci60_TF)


#' load TF activity inference data from CPATC/CCLE dataset
ccle_TF <- read_csv(file = "./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "tf_activity") %>%
  rename(tf=X1) %>%
  select(sample, tf, tf_activity) %>%
  inner_join(cptac_samples[, c("sample", "batch")], by = "sample") %>%
  filter(batch == "ccle") %>%
  select(-batch) %>%
  anti_join(nci60_crc65_TF, by = "sample") #%>%
  #anti_join(alternative_names, by = c("sample" = "alternative_names"))


#' bind the TF activities of both datasets
tf_activity <- bind_rows(ccle_TF, nci60_crc65_TF)


#' load NCI60/CRC65 mutations
nci60_crc65_mutation <- read_tsv(file = "./output/files/nci60_crc65_mutations.txt.gz")


#' load Roumeliotis mutations
roum_mutation <- read_tsv(file = "./output/files/mutations_protpos.txt.gz") %>%
  inner_join(cptac_samples[, c("sample", "batch")], by = "sample") %>%
  filter(batch == "ccle") %>%
  select(-batch) %>%
  anti_join(nci60_crc65_TF, by = "sample") %>%
  anti_join(alternative_names, by = c("sample" = "alternative_names"))


#' load protein sequences
prot_seqs <- read_tsv(file = "./output/files/protein_sequences_GRCh37.txt") %>%
  filter(state == "known" & gene_biotype == "protein_coding" & trpt_biotype == "protein_coding") %>%
  select(trpt_id, protein_id, seq)


#' load canonical uniprot transcripts
canonical_trpt <- read_tsv(file = "./data/uniprot/canonical_trpt/isoform_overrides_uniprot.txt", skip = 1)


#' # -- assess the impact of missense mutations on TF activities

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


#' select missense mutations
missense_nci60_crc65 <- nci60_crc65_mutation %>%
  filter(mutation_type == "Missense") %>%
  select(sample, gene, identifier, HGVSp, prot_pos_i, aa_wt, aa_mut)

missense_roum <- roum_mutation %>%
  filter(variant_class == "Missense_Mutation") %>%
  select(sample, gene=gene_symbol, CHROM, POS, REF, ALT, HGVSp, prot_pos_i = prot_pos, aa_wt, aa_mut) %>%
  unite(col = "identifier","CHROM", "POS", "REF", "ALT")


# function to convert sift scores to a discrete scale
sift_to_factor <- function(x){
  if(x < 0 | (x >= 0 & x < 1)){
    val <- "1"
  } else if(x >= 1 & x < 2){
    val <- "2"
  } else if(x >= 2 & x < 3){
    val <- "3"
  } else if(x >= 3 & x < 4){
    val <- "4"
  } else if(x >= 4 & x < 5){
    val <- "5"
  } else if(x >= 5){
    val <- "6"
  }
  val
}


#' scatterplot
missense_impact <- missense_nci60_crc65 %>%
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

ggsave(filename = "tf_act_missense_mutation_sift_correlation1.png", plot = missense_impact, path = "./output/plots/nci60_crc65/", width = 6, height = 6)
unlink("tf_act_missense_mutation_sift_correlation1.png")


#' boxplot
missense_impact <- missense_nci60_crc65 %>%
  inner_join(tf_activity, by = c("sample", "gene" = "tf")) %>%
  inner_join(sift, by = c("gene", "prot_pos_i"="pos","aa_wt"="ref","aa_mut"="alt")) %>%
  group_by(sample, gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  mutate(sift = -log10(score+1e-10)) %>%
  #filter(sift < 5) %>%
  mutate(sift_factor = map_chr(.x=sift, .f=sift_to_factor)) %>%
  ggplot(mapping = aes(x = sift_factor, y = tf_activity)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, alpha=0.1)

#+ fig.width=6, fig.height=6
missense_impact

ggsave(filename = "tf_act_missense_mutation_sift_boxplot1.png", plot = missense_impact, path = "./output/plots/nci60_crc65/", width = 6, height = 6)
unlink("tf_act_missense_mutation_sift_boxplot1.png")


#' scatterplot
missense_impact <- missense_nci60_crc65 %>%
  bind_rows(missense_roum) %>%
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

ggsave(filename = "tf_act_missense_mutation_sift_correlation2.png", plot = missense_impact, path = "./output/plots/nci60_crc65/", width = 6, height = 6)
unlink("tf_act_missense_mutation_sift_correlation2.png")


#' boxplot
missense_impact <- missense_nci60_crc65 %>%
  bind_rows(missense_roum) %>%
  inner_join(tf_activity, by = c("sample", "gene" = "tf")) %>%
  inner_join(sift, by = c("gene", "prot_pos_i"="pos","aa_wt"="ref","aa_mut"="alt")) %>%
  group_by(sample, gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  mutate(sift = -log10(score+1e-10)) %>%
  #filter(sift < 5) %>%
  mutate(sift_factor = map_chr(.x=sift, .f=sift_to_factor)) %>%
  ggplot(mapping = aes(x = sift_factor, y = tf_activity)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width=0.1, alpha=0.1)

#+ fig.width=6, fig.height=6
missense_impact

ggsave(filename = "tf_act_missense_mutation_sift_boxplot2.png", plot = missense_impact, path = "./output/plots/nci60_crc65/", width = 6, height = 6)
unlink("tf_act_missense_mutation_sift_boxplot2.png")
