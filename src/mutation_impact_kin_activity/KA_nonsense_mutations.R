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


# correlate kinase activity with the position of nonsense mutations
# kinase activity with imputations

# regress out batch covariate from kinase activities
ka_regout <- ka2 %>%
  inner_join(samples_annotation[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(kin_activity_resd = map(.x = data, .f = ~ residuals(lm(kin_activity ~ batch, data = .x)))) %>%
  unnest() %>%
  select(kinase, sample, kin_activity = kin_activity_resd)

# select nonsense mutations
nonsense <- mut %>%
  filter(variant_class == "Nonsense_Mutation") %>%
  mutate(prot_length = nchar(seq)) %>%
  select(sample, gene = gene_symbol, prot_pos, prot_length) %>%
  distinct() %>%
  # remove cases of samples with multiple mutations on the same gene
  group_by(sample, gene) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)

ka_mut <- ka_regout %>%
  inner_join(nonsense, by = c("kinase" = "gene", "sample")) %>%
  mutate(prot_pos_p = prot_pos/prot_length) %>%
  select(-prot_length) %>%
  pivot_longer(c(prot_pos, prot_pos_p), names_to = "position", values_to = "value") %>%
  ggplot(mapping = aes(x = value, y = kin_activity)) +
  facet_wrap( ~ position, scales = "free_x", labeller = labeller(position = c("prot_pos" = "Protein position (aa)", "prot_pos_p" = "Protein position (%)"))) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  theme_classic() +
  ggpubr::stat_cor(size = 3) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10),
    panel.spacing = unit(1, "lines")) +
  labs(x = "Nonsense mutation position", y = "Kinase activity")


# correlate kinase activity with the position of nonsense mutations
# kinase activity without imputations

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

# regress out batch covariate from kinase activities
ka_regout <- ka1 %>%
  inner_join(samples_annotation[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, ~ nrow(.x) > 10)) %>%
  mutate(kin_activity_resd = map(.x = data, .f = reg)) %>%
  unnest() %>%
  select(kinase, sample, kin_activity = kin_activity_resd)

# select nonsense mutations
nonsense <- mut %>%
  filter(variant_class == "Nonsense_Mutation") %>%
  mutate(prot_length = nchar(seq)) %>%
  select(sample, gene = gene_symbol, prot_pos, prot_length) %>%
  distinct() %>%
  # remove cases of samples with multiple mutations on the same gene
  group_by(sample, gene) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)

ka_mut <- ka_regout %>%
  inner_join(nonsense, by = c("kinase" = "gene", "sample")) %>%
  mutate(prot_pos_p = prot_pos/prot_length) %>%
  select(-prot_length) %>%
  pivot_longer(c(prot_pos, prot_pos_p), names_to = "position", values_to = "value") %>%
  ggplot(mapping = aes(x = value, y = kin_activity)) +
  facet_wrap( ~ position, scales = "free_x", labeller = labeller(position = c("prot_pos" = "Protein position (aa)", "prot_pos_p" = "Protein position (%)"))) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  theme_classic() +
  ggpubr::stat_cor(size = 3) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10),
    panel.spacing = unit(1, "lines")) +
  labs(x = "Nonsense mutation position", y = "Kinase activity")
