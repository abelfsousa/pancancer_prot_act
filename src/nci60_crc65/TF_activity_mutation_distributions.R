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


#' # -- assess the TF activity distributions between different mutation types
mut1 <- nci60_crc65_mutation %>%
  select(sample, gene, mutation_type) %>%
  distinct()

mut2 <- roum_mutation %>%
  select(sample, gene=gene_symbol, mutation_type=variant_class) %>%
  filter(!mutation_type == "Splice_Site") %>%
  mutate(mutation_type = str_replace(mutation_type, "Missense_Mutation", "Missense")) %>%
  mutate(mutation_type = str_replace(mutation_type, "Nonsense_Mutation", "Nonsense")) %>%
  mutate(mutation_type = str_replace(mutation_type, "Frame_Shift_Ins|Frame_Shift_Del", "Frameshift")) %>%
  mutate(mutation_type = str_replace(mutation_type, "In_Frame_Del|In_Frame_Ins", "Nonframeshift")) %>%
  mutate(mutation_type = str_replace(mutation_type, "Nonstop_Mutation", "Read_through")) %>%
  distinct()


#' ## using only NCI60/CRC65 datasets
mutations <- mut1 %>%
  semi_join(tf_activity, by = c("gene" = "tf", "sample")) %>%
  group_by(mutation_type) %>%
  nest() %>%
  ungroup() #%>%
  #mutate(data2 = (function(x) map2(.x=1:length(x),.y=x,.f=function(a,b)reduce(.x=c(list(b),x[-a]),.f=setdiff))) (data))

mutations2 <- mutations %>%
  mutate(data2 = map2(.x = mutation_type, .y = data, .f = function(x,y) reduce(.x=c(list(y), mutations[mutations$mutation_type != x, "data", drop = T]), .f=setdiff))) %>%
  select(-data) %>%
  unnest()

mutations <- mutations %>%
  unnest()

distribution <- tf_activity %>%
  inner_join(mutations2, by = c("sample", "tf" = "gene")) %>%
  ggplot(mapping = aes(x=mutation_type, y = tf_activity, fill = mutation_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, show.legend = F) +
  geom_jitter(size = 1, width = 0.1, alpha = 0.1, show.legend = F) +
  ggpubr::stat_compare_means(comparisons = list(c(1,6),c(1,3))) +
  coord_flip()

#+ fig.width=7, fig.height=3
distribution

distribution <- tf_activity %>%
  inner_join(mutations2, by = c("sample", "tf" = "gene")) %>%
  ggplot(mapping = aes(x=tf_activity, fill = mutation_type, color = mutation_type)) +
  geom_density(alpha = 0.5)

#+ fig.width=7, fig.height=2
distribution

frameshift <- tf_activity %>%
  inner_join(mutations2, by = c("sample", "tf" = "gene")) %>%
  filter(mutation_type == "Frameshift") %>%
  arrange((tf_activity))


#' ## using NCI60/CRC65 datasets + Roumeliotis
mutations <- bind_rows(mut1, mut2) %>%
  semi_join(tf_activity, by = c("gene" = "tf", "sample")) %>%
  group_by(mutation_type) %>%
  nest() %>%
  ungroup()

mutations2 <- mutations %>%
  mutate(data2 = map2(.x = mutation_type, .y = data, .f = function(x,y) reduce(.x=c(list(y), mutations[mutations$mutation_type != x, "data", drop = T]), .f=setdiff))) %>%
  select(-data) %>%
  unnest()

mutations <- mutations %>%
  unnest()

distribution <- tf_activity %>%
  inner_join(mutations2, by = c("sample", "tf" = "gene")) %>%
  ggplot(mapping = aes(x=mutation_type, y = tf_activity, fill = mutation_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, show.legend = F) +
  geom_jitter(size = 2, width = 0.1, alpha = 0.1, show.legend = F) +
  ggpubr::stat_compare_means(comparisons = list(c(1,6),c(1,3))) +
  coord_flip()

#+ fig.width=7, fig.height=3
distribution

distribution <- tf_activity %>%
  inner_join(mutations2, by = c("sample", "tf" = "gene")) %>%
  ggplot(mapping = aes(x=tf_activity, fill = mutation_type, color = mutation_type)) +
  geom_density(alpha = 0.5)

#+ fig.width=7, fig.height=2
distribution
