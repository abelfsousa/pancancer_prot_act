# assess the correlation between kinase/TF activities and the SIFT score of missense mutations

# load R packages
library(ggpubr)
library(tidyverse)

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load SIFT data
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
  select(gene, pos, ref, alt, score)


# -- for CPTAC tumours


# load CPTAC samples
source("./src/utils/getSamples.R")
cptac_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
  getSamples("protein") %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue)


# load CPTAC mutation data
cptac_mutations <- read_tsv(file = "./output/files/mutations_protpos.txt.gz")


# load CPTAC TF activities
tf_activity <- read_csv(file = "./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "activity") %>%
  rename(tf=X1) %>%
  select(sample, gene = tf, activity)


# load CPTAC kinase activities (with quantile-normalized protein regressed-out phosphorylation data)
k_subN <- 3
kin_activity <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  select(sample, gene = kinase, activity = log10P)


# select only missense mutations
# select only the sample, gene, protein position, wt aa and mutated aa
# remove duplicates (same gene in the same sample can have multiple SNPs in different genomic positions leading to the same protein alteration (same missense mutation))
missense <- cptac_mutations %>%
  filter(variant_class == "Missense_Mutation") %>%
  select(sample, gene=gene_symbol, prot_pos_i = prot_pos, aa_wt, aa_mut) %>%
  distinct()


# map missense mutations to kinase activities and SIFT scores
# remove cases where the same sample has multiple missense mutations
# in the same gene to prevent the assignment of the same activity the different SIFT scores
sift_kin <- sift %>%
  filter(gene %in% kin_activity$gene)

missense_kinases <- missense %>%
  inner_join(kin_activity, by = c("sample", "gene")) %>%
  inner_join(sift_kin, by = c("gene", "prot_pos_i"="pos","aa_wt"="ref","aa_mut"="alt")) %>%
  group_by(sample, gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  mutate(sift = -log10(score+1e-10)) %>%
  mutate(protein = "Kinases")


# map missense mutations to TF activities
# remove cases where the same sample has multiple missense mutations
# in the same gene to prevent the assignment of the same activity the different SIFT scores
sift_tf <- sift %>%
  filter(gene %in% tf_activity$gene)

missense_tfs <- missense %>%
  inner_join(tf_activity, by = c("sample", "gene")) %>%
  inner_join(sift_tf, by = c("gene", "prot_pos_i"="pos","aa_wt"="ref","aa_mut"="alt")) %>%
  group_by(sample, gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  mutate(sift = -log10(score+1e-10)) %>%
  mutate(protein = "TFs")


missense_activities <- bind_rows(missense_kinases, missense_tfs)


correlation_plot <- missense_activities %>%
  ggplot(mapping = aes(x = sift, y = activity)) +
  geom_point(size = 1.5, alpha = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~ protein, scales = "free") +
  stat_cor(method = "pearson", label.y.npc = 0.9, label.x.npc = 0.45, size = 7) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, colour = "black", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 22, colour = "black"),
    axis.title = element_text(size = 22, colour = "black"),
    axis.text = element_text(size = 20, colour = "black")) +
  labs(x = "SIFT score (-log10)", y = "Activity score", title = "CPTAC tumours")

ggsave(filename = "cptac_kin_tf_act_missense_sift_corr.png", plot = correlation_plot, path = "./output/plots/mutation_impact_kin_tf_activity/", width = 12, height = 5)
ggsave(filename = "cptac_kin_tf_act_missense_sift_corr.pdf", plot = correlation_plot, path = "./output/plots/mutation_impact_kin_tf_activity/", width = 12, height = 5)


# load CPTAC tumours purity scores
purity_scores <- read_tsv(file = "./output/files/samples_purity.txt")
pure_samples <- purity_scores %>%
  filter(score > 0.8)

correlation_plot <- missense_activities %>%
  filter(sample %in% pure_samples$sample) %>%
  ggplot(mapping = aes(x = sift, y = activity)) +
  geom_point(size = 1.5, alpha = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~ protein, scales = "free") +
  stat_cor(method = "pearson", label.y.npc = 0.9, label.x.npc = 0.45, size = 7) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, colour = "black", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 22, colour = "black"),
    axis.title = element_text(size = 22, colour = "black"),
    axis.text = element_text(size = 20, colour = "black")) +
  labs(x = "SIFT score (-log10)", y = "Activity score", title = "CPTAC tumours (purity > 0.8)")

ggsave(filename = "cptac_kin_tf_act_missense_sift_corr_pure_samples.png", plot = correlation_plot, path = "./output/plots/mutation_impact_kin_tf_activity/", width = 12, height = 5)
ggsave(filename = "cptac_kin_tf_act_missense_sift_corr_pure_samples.pdf", plot = correlation_plot, path = "./output/plots/mutation_impact_kin_tf_activity/", width = 12, height = 5)




# -- for NCI60/CRC65 cell lines


# load NCI60/CRC65 cell lines
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


# load NCI60/CRC65 mutations
cells_mutations <- read_tsv(file = "./output/files/nci60_crc65_mutations.txt.gz")


# load TF activity inference data from NCI60/CRC65 dataset
tf_activity_crc65 <- read_csv(file = "./data/dorothea/nci60_crc65/TF_activity_log2FC_allsamples_CRC65.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "activity") %>%
  rename(gene=X1) %>%
  filter(!sample %in% shared_cell_lines) %>%
  select(sample, gene, activity)

tf_activity_nci60 <- read_csv(file = "./data/dorothea/nci60_crc65/TF_activity_log2FC_allsamples_NCI60.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "activity") %>%
  rename(gene=X1) %>%
  select(sample, gene, activity)

tf_activity_cells <- bind_rows(tf_activity_crc65, tf_activity_nci60)


# load kinase activity inference data from NCI60/CRC65 dataset
# select quantifications with more than 3 substrates
kin_activity_cells <- read_tsv(file = "./output/files/nci60_crc65_kinase_activities.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= 3) %>%
  select(-source_type, -n) %>%
  filter(!(batch == "CRC65" & sample %in% shared_cell_lines)) %>%
  select(sample, gene = kinase, activity = log10P)


# select only missense mutations
# select only the sample, gene, protein position, wt aa and mutated aa
# remove duplicates (same gene in the same sample can have multiple SNPs in different genomic positions leading to the same protein alteration (same missense mutation))
missense_cells <- cells_mutations %>%
  filter(mutation_type == "Missense") %>%
  select(sample, gene, prot_pos_i, aa_wt, aa_mut) %>%
  distinct()


# map missense mutations to kinase activities and SIFT scores
# remove cases where the same sample has multiple missense mutations
# in the same gene to prevent the assignment of the same activity the different SIFT scores
sift_kin_cells <- sift %>%
  filter(gene %in% kin_activity_cells$gene)

missense_kinases_cells <- missense_cells %>%
  inner_join(kin_activity_cells, by = c("sample", "gene")) %>%
  inner_join(sift_kin_cells, by = c("gene", "prot_pos_i"="pos","aa_wt"="ref","aa_mut"="alt")) %>%
  group_by(sample, gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  mutate(sift = -log10(score+1e-10)) %>%
  mutate(protein = "Kinases")


# map missense mutations to TF activities
# remove cases where the same sample has multiple missense mutations
# in the same gene to prevent the assignment of the same activity the different SIFT scores
sift_tf_cells <- sift %>%
  filter(gene %in% tf_activity_cells$gene)

missense_tfs_cells <- missense_cells %>%
  inner_join(tf_activity_cells, by = c("sample", "gene")) %>%
  inner_join(sift_tf_cells, by = c("gene", "prot_pos_i"="pos","aa_wt"="ref","aa_mut"="alt")) %>%
  group_by(sample, gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  mutate(sift = -log10(score+1e-10)) %>%
  mutate(protein = "TFs")


missense_activities_cells <- bind_rows(missense_kinases_cells, missense_tfs_cells)


correlation_plot <- missense_activities_cells %>%
  ggplot(mapping = aes(x = sift, y = activity)) +
  geom_point(size = 1.5, alpha = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~ protein, scales = "free") +
  stat_cor(method = "pearson", label.y.npc = 0.9, label.x.npc = 0.45, size = 7) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, colour = "black", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 22, colour = "black"),
    axis.title = element_text(size = 22, colour = "black"),
    axis.text = element_text(size = 20, colour = "black")) +
  labs(x = "SIFT score (-log10)", y = "Activity score", title = "NCI60/CRC65 cancer cell lines")

ggsave(filename = "nci60_crc65_kin_tf_act_missense_sift_corr.png", plot = correlation_plot, path = "./output/plots/mutation_impact_kin_tf_activity/", width = 12, height = 5)
ggsave(filename = "nci60_crc65_kin_tf_act_missense_sift_corr.pdf", plot = correlation_plot, path = "./output/plots/mutation_impact_kin_tf_activity/", width = 12, height = 5)
