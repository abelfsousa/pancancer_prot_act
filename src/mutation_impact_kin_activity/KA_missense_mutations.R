# load R packages
library(tidyverse)

source("./src/utils/getSamples.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load samples annotation
samples_annotation <- read_tsv(file = "./output/files/all_samples_annotation.txt")
samples_annotation <- getSamples(samples_annotation, c("phosphorylation"))


# load mutation data
mut <- read_tsv(file = "./output/files/mutations_protpos.txt.gz")

# select missense mutations
missense <- mut %>%
  filter(variant_class == "Missense_Mutation" & variant_type == "SNP") %>%
  select(sample, gene = gene_symbol, pos = prot_pos, ref = aa_wt, alt = aa_mut) %>%
  distinct() %>%
  # remove cases of samples with multiple mutations on the same gene
  group_by(sample, gene) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)


# load kinase-activity inference data without imputations
# (quantile-normalized protein regressed-out phosphorylation data)
# select quantifications with more than 3 substrates
k_subN <- 3
ka1 <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  select(sample, kinase, kin_activity=log10P)


# load kinase activities and imputed values
# matrix used for PCA analysis and UMAP
ka2 <- data.table::fread("./data/Danish/kinaseActMatImputed.tsv") %>%
  as_tibble() %>%
  rename(kinase = V1) %>%
  pivot_longer(-kinase, names_to = "sample", values_to = "kin_activity")


# regress out batch covariate from kinase activities

# function to fit a linear model and extract residuals
reg <- function(df){
  if(length(unique(df$batch)) == 1){
    res <- df$kin_activity
  } else {
    # fit the model
    m <- lm(kin_activity ~ batch, data = df)
    
    # extract residuals
    res <- residuals(m)
  }
  res
}

ka1_regout <- ka1 %>%
  inner_join(samples_annotation[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, ~ nrow(.x) > 10)) %>%
  mutate(kin_activity_resd = map(.x = data, .f = reg)) %>%
  unnest() %>%
  select(kinase, sample, kin_activity = kin_activity_resd)

ka2_regout <- ka2 %>%
  inner_join(samples_annotation[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(kin_activity_resd = map(.x = data, .f = reg)) %>%
  unnest() %>%
  select(kinase, sample, kin_activity = kin_activity_resd)


# load SIFT data

# uniprot to gene name mapping
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
  filter(!gene %in% exclude$gene) %>%
  mutate(score = -log10(score+1e-10)) %>%
  select(gene, pos, ref, alt, score)


# correlate SIFT score to kinase activity
# imputed values and batch regressed-out
sift2 <- sift %>%
  filter(gene %in% ka2_regout$kinase)

ka_sift <- ka2_regout %>%
  inner_join(missense, by = c("sample", "kinase" = "gene")) %>%
  inner_join(sift2, by = c("kinase" = "gene", "pos", "ref", "alt")) %>%
  filter(score <= 5) %>%
  ggplot(mapping = aes(x = score, y = kin_activity)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  theme_classic() +
  ggpubr::stat_cor(size = 3) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10),
    aspect.ratio = 1) +
  labs(x = "SIFT score (-log10)", y = "Kinase activity")

ggsave(filename = "kinase_activity_imputed_batch_reg_out_missense_sift_cor.png", plot = ka_sift, path = "./output/plots/mutation_impact_kin_activity/", height = 4, width = 4)
unlink("kinase_activity_imputed_batch_reg_out_missense_sift_cor.png")


# correlate SIFT score to kinase activity
# without imputed values and batch regressed-out
sift2 <- sift %>%
  filter(gene %in% ka1_regout$kinase)

ka_sift <- ka1_regout %>%
  inner_join(missense, by = c("sample", "kinase" = "gene")) %>%
  inner_join(sift2, by = c("kinase" = "gene", "pos", "ref", "alt")) %>%
  filter(score <= 5) %>%
  ggplot(mapping = aes(x = score, y = kin_activity)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  theme_classic() +
  ggpubr::stat_cor(size = 3) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10),
    aspect.ratio = 1) +
  labs(x = "SIFT score (-log10)", y = "Kinase activity")

ggsave(filename = "kinase_activity_not_imputed_batch_reg_out_missense_sift_cor.png", plot = ka_sift, path = "./output/plots/mutation_impact_kin_activity/", height = 4, width = 4)
unlink("kinase_activity_not_imputed_batch_reg_out_missense_sift_cor.png")
