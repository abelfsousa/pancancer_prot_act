# Correlation between TF activities and TF phosphosites

library(tidyverse)
library(ggpubr)
library(viridis)
library(Hmisc)
library(ggbeeswarm)
source(file = "./src/utils/getSamples.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy

set.seed(123)


# load TF activities
tf_act <- data.table::fread("./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  as_tibble() %>%
  rename(tf=V1) %>%
  pivot_longer(-tf, names_to = "sample", values_to = "tf_activity")


# load CPTAC phosphorylation data
phospho <- read_tsv(file = "./output/files/phosphoproteomicsQ_Protreg_allSamp_withZtransf.txt.gz") %>%
#phospho <- read_tsv(file = "./output/files/phosphoproteomicsQ.txt.gz") %>%
  pivot_longer(-c(gene, psite, psites), names_to = "sample", values_to = "log2fc") %>%
  filter(!is.na(log2fc)) %>%
  filter(psites == 1) %>%
  select(-gene, -psites) %>%
  separate(col = "psite", into = c("gene", "psite"), sep = "_") %>%
  mutate(residue = str_extract(psite, "^[STY]{1}")) %>%
  mutate(position = as.numeric(str_extract(psite, "[0-9]+"))) %>%
  select(-psite) %>%
  select(gene, position, residue, everything())


# load TF phosphosites with respective functional score
regPsite <- read_tsv(file = "./output/files/cptac_tf_psites_funscoR.txt.gz")


# load PhosphoSitePlus regulatory sites annotation

# function to evaluate the phosphosites regulatory status
status <- function(strg){
  if (is.na(strg)){
    "unknown"
  }
  else if(str_detect(string = strg, pattern = "inhibited") & str_detect(string = strg, pattern = "induced")){
    "both"
  } else if(str_detect(string = strg, pattern = "inhibited")){
    "inhibition"
  } else if(str_detect(string = strg, pattern = "induced")){
    "activation"
  } else{
    "unknown"
  }
}

psp_reg <- read_tsv(file = "./data/PhosphoSitePlus/Regulatory_sites.gz", skip = 3) %>%
  filter(ORGANISM == "human") %>%
  select(gene=GENE, psite=MOD_RSD, reg_fun = ON_FUNCTION) %>%
  na.exclude() %>%
  separate(col = psite, into = c("psite", "modification")) %>%
  filter(modification == "p") %>%
  select(-modification) %>%
  mutate(residue = str_extract(psite, "^[A-Z]{1}"), position = str_extract(psite, "[0-9]+")) %>%
  select(-psite) %>%
  select(gene, residue, position, reg_fun) %>%
  distinct() %>%
  mutate(position = as.numeric(position)) %>%
  mutate(reg_status = map_chr(.x = reg_fun, .f = status)) %>%
  select(-reg_fun)


# load metadata
samples_metadata <- read_tsv("./output/files/all_samples_annotation.txt")
samples_metadata <- getSamples(samples_metadata, data_types = c("mRNA", "protein")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue)


# correlate TF activities with CPTAC TF phosphosites on the same TF

# set a correlation function
corr <- function(df, x, y, method){
  a = df[[x]]
  b = df[[y]]
  cor <- broom::tidy(cor.test(a, b, method = method))
  
  res <- cor %>% select(estimate, p.value)
  
  res
}


TFact_RegPsite_cor <- tf_act %>%
  inner_join(regPsite, by = c("tf" = "gene")) %>%
  select(sample, tf, activity=tf_activity, position, residue, funscoR=probabilities) %>%
  inner_join(phospho, by = c("sample", "tf" = "gene", "position", "residue")) %>%
  select(sample, tf, activity, position, residue, log2fc, everything()) %>%
  group_by(tf, position, residue, funscoR) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(data = map(.x = data, .f = ~ .x %>% mutate(activity = scale(activity)[,1], log2fc = scale(log2fc)[,1]))) %>%
  mutate(cor = map(.x = data, .f = corr, x="activity", y = "log2fc", method = "pearson")) %>%
  select(-data) %>%
  unnest() %>%
  rename(pearson_r=estimate) %>%
  mutate(class1 = "TF-Phosphosite", class2="Substrate-based") %>%
  rename(TFA=tf) %>%
  mutate(TFB=TFA)


TFact_RegPsite_cor_batch <- tf_act %>%
  inner_join(regPsite, by = c("tf" = "gene")) %>%
  select(sample, tf, activity=tf_activity, position, residue, funscoR=probabilities) %>%
  inner_join(phospho, by = c("sample", "tf" = "gene", "position", "residue")) %>%
  inner_join(samples_metadata, by = "sample") %>%
  select(batch, sample, tf, activity, position, residue, log2fc, everything()) %>%
  group_by(batch, tf, position, residue, funscoR) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 5)) %>%
  mutate(data = map(.x = data, .f = ~ .x %>% mutate(activity = scale(activity)[,1], log2fc = scale(log2fc)[,1]))) %>%
  mutate(cor = map(.x = data, .f = corr, x="activity", y = "log2fc", method = "pearson")) %>%
  select(-data) %>%
  unnest() %>%
  rename(pearson_r=estimate) %>%
  mutate(class1 = "TF-Phosphosite", class2="Substrate-based") %>%
  rename(TFA=tf) %>%
  mutate(TFB=TFA)


# correlate TF activities with all CPTAC TF phosphosites (all pairs)

# set up a function to convert a squared matrix into a data frame (tibble)
matrix_to_tb <- function(x, y, k, p){
  
  mat <- x
  mat_data <- y
  KIN <- k
  PHO <- p
  
  mat <- mat[p,k]
  df <- mat %>%
    as.data.frame() %>%
    rownames_to_column(var="kinaseA") %>%
    as_tibble() %>%
    pivot_longer(-c("kinaseA"), names_to="kinaseB", values_to=mat_data) %>%
    separate(col = "kinaseA", into = c("kinaseA", "position", "residue"), sep = "_") %>%
    mutate(position = as.numeric(position))
  
  df
}

# set up a function to join the rcorr output into a single tibble
joinMatrices <- function(x, kin, pho){
  
  rcorr_matrices <- x
  rcorr_names <- names(rcorr_matrices)
  
  KIN <- kin
  PHO <- pho
  
  dfs <- map2(.x = rcorr_matrices, .y = rcorr_names, .f = matrix_to_tb, k = KIN, p = PHO)
  
  tb <- reduce(.x = dfs, .f = ~ inner_join(.x, .y, by = c("kinaseA", "position", "residue", "kinaseB")))
  
  tb
}


tf_mat <- tf_act %>%
  semi_join(regPsite, by = c("tf" = "gene")) %>%
  semi_join(phospho, by = c("sample")) %>%
  arrange(sample) %>%
  pivot_wider(names_from = "sample", values_from = "tf_activity") %>%
  as.data.frame() %>%
  column_to_rownames(var = "tf") %>%
  as.matrix() %>%
  t()

tfs <- colnames(tf_mat)

phospho_mat <- phospho %>%
  semi_join(regPsite, by = c("gene", "position", "residue")) %>%
  semi_join(tf_act, by = c("sample")) %>%
  arrange(sample) %>%
  pivot_wider(names_from = "sample", values_from = "log2fc") %>%
  unite(gene, position, residue, col = psite) %>%
  as.data.frame() %>%
  column_to_rownames(var = "psite") %>%
  as.matrix() %>%
  t()

psites <- colnames(phospho_mat)

phospho_mat <- scale(phospho_mat)
tf_mat <- scale(tf_mat)

TFact_RegPsite_cor_all <- rcorr(phospho_mat, tf_mat) %>%
  joinMatrices(tfs, psites) %>%
  filter(n > 10) %>%
  select(-n, pearson_r=r, p.value=P) %>%
  inner_join(regPsite, by = c("kinaseA"="gene", "position", "residue")) %>%
  rename(funscoR=probabilities) %>%
  mutate(class1 = "All pairs", class2="Substrate-based") %>%
  filter(kinaseA != kinaseB) %>%
  rename(TFA = kinaseA, TFB = kinaseB)


# join all correlations
all_cor <- bind_rows(TFact_RegPsite_cor, TFact_RegPsite_cor_all) %>%
  select(TFB, TFA, position, residue, everything())


# plot the distributions

TFact_RegPsite_cor_distribution <- all_cor %>%
  filter(class1 == "TF-Phosphosite", class2 == "Substrate-based") %>%
  left_join(psp_reg, by = c("TFA" = "gene", "residue", "position")) %>%
  mutate(reg_status = replace_na(reg_status, "unknown")) %>%
  #filter(reg_status %in% c("inhibition", "activation")) %>%
  ggplot(mapping = aes(x = reg_status, y = pearson_r)) +
  geom_boxplot(notch = F, color = "black", fill = "#f0f0f0", show.legend = F) +
  geom_jitter(mapping = aes(color = reg_status), width = 0.1, alpha = 0.4, show.legend = F) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c(1,4), c(1,3), c(2,4), c(3,4))) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 18),
    axis.text = element_text(colour = "black", size = 16)) +
  labs(x = "PhosphositePlus annotation", y = "Pearson's r (TF activity vs TF phosphosite)")

ggsave(filename = "tf_phosphosite_correlation_psp_status_box.pdf", plot = TFact_RegPsite_cor_distribution, path = "./output/plots/tf_reg_psites_correlation/", height = 6, width = 6)
ggsave(filename = "tf_phosphosite_correlation_psp_status_box.png", plot = TFact_RegPsite_cor_distribution, path = "./output/plots/tf_reg_psites_correlation/", height = 6, width = 6)


TFact_RegPsite_cor_distribution <- all_cor %>%
  filter(class1 == "TF-Phosphosite", class2 == "Substrate-based") %>%
  mutate(group1 = cut_number(funscoR, 5)) %>%
  mutate(group2 = cut(funscoR, 5)) %>%
  ggplot(mapping = aes(x = group2, y = pearson_r)) +
  geom_boxplot(notch = T, color = "black", fill = "#f0f0f0", show.legend = F) +
  geom_jitter(mapping = aes(color = group2), width = 0.1, alpha = 0.4, show.legend = F) +
  #stat_compare_means(method = "wilcox.test", comparisons = list(c(1,4), c(1,3), c(2,4), c(3,4))) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)) +
  labs(x = "Phosphosite functional score", y = "Pearson's r (Activity vs Phosphorylation)")

ggsave(filename = "tf_phosphosite_correlation_box_funscore.pdf", plot = TFact_RegPsite_cor_distribution, path = "./output/plots/tf_reg_psites_correlation/", height = 4, width = 6)
ggsave(filename = "tf_phosphosite_correlation_box_funscore.png", plot = TFact_RegPsite_cor_distribution, path = "./output/plots/tf_reg_psites_correlation/", height = 4, width = 6)


# same
# TFact_RegPsite_cor_distribution <- all_cor %>%
#   filter(class1 == "TF-Phosphosite", class2 == "Substrate-based") %>%
#   mutate(group1 = cut(funscoR, 5, dig.lab = 1)) %>%
#   #filter(abs(pearson_r) > 0.4) %>%
#   mutate(signf = if_else(p.value < 0.05, "sig", "nonsig")) %>%
#   group_by(group1, signf) %>%
#   summarise(count = n()) %>%
#   ungroup() %>%
#   ggplot(mapping = aes(x = group1, y = count, fill = signf)) +
#   geom_col(position = "fill", width = 1, color = "white")
# TFact_RegPsite_cor_distribution
# 
# TFact_RegPsite_cor_distribution <- all_cor %>%
#   filter(class1 == "TF-Phosphosite", class2 == "Substrate-based") %>%
#   mutate(group1 = cut(funscoR, 5, dig.lab = 1)) %>%
#   #filter(abs(pearson_r) > 0.4) %>%
#   mutate(signf = if_else(p.value < 0.05, "sig", "nonsig")) %>%
#   ggplot(mapping = aes(x = group1, fill = signf)) +
#   geom_bar(position = "fill", width = 1, color = "white")
# TFact_RegPsite_cor_distribution
# 
# # same
# TFact_RegPsite_cor_distribution <- all_cor %>%
#   filter(class1 == "TF-Phosphosite", class2 == "Substrate-based") %>%
#   mutate(group1 = cut(funscoR, 5, dig.lab = 1)) %>%
#   #filter(abs(pearson_r) > 0.4) %>%
#   mutate(signf = if_else(p.value < 0.05, "sig", "nonsig")) %>%
#   group_by(group1, signf) %>%
#   summarise(count = n()) %>%
#   ungroup() %>%
#   group_by(group1) %>%
#   mutate(perc = count/sum(count)) %>%
#   ungroup() %>%
#   ggplot(mapping = aes(x = group1, y = perc, fill = signf)) +
#   geom_col(width = 1, color = "white")

TFact_RegPsite_cor_distribution <- all_cor %>%
  filter(class1 == "TF-Phosphosite", class2 == "Substrate-based") %>%
  mutate(group1 = cut(funscoR, 5, dig.lab = 1)) %>%
  #filter(abs(pearson_r) > 0.4) %>%
  group_by(group1) %>%
  summarise(sig = sum(p.value < 0.05), nonsig = sum(p.value > 0.05)) %>%
  ungroup() %>%
  pivot_longer(-group1, names_to = "signf", values_to = "count") %>%
  group_by(group1) %>%
  mutate(perc = count/sum(count)) %>%
  ungroup()

counts <- TFact_RegPsite_cor_distribution %>%
  group_by(group1) %>%
  summarise(count = sum(count)) %>%
  ungroup()

TFact_RegPsite_cor_distribution <- TFact_RegPsite_cor_distribution %>%
  ggplot(mapping = aes(x = group1, y = perc, fill = signf)) +
  geom_col(color = "white", width = 1) +
  geom_text(data = counts, mapping = aes(x = group1, y = 0.95, label = count), size = 5, inherit.aes = F) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 18, color = "black"),
    legend.text = element_text(size = 16, color = "black"),
    legend.position = "bottom") +
  scale_fill_manual(values = c("nonsig" = "#deebf7", "sig" = "#3182bd"), labels = c("nonsig" = "P > 5%", "sig" = "P < 5%")) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  #scale_y_continuous(labels = scales::percent) +
  labs(x = "Functional score (binned)", y = "Percentage", fill = "Pearson's r") +
  guides(fill = guide_legend(nrow = 1))

ggsave(filename = "tf_phosphosite_correlation_barplot_funscore.pdf", plot = TFact_RegPsite_cor_distribution, path = "./output/plots/tf_reg_psites_correlation/", height = 4, width = 5)
ggsave(filename = "tf_phosphosite_correlation_barplot_funscore.png", plot = TFact_RegPsite_cor_distribution, path = "./output/plots/tf_reg_psites_correlation/", height = 4, width = 5)


TFact_RegPsite_cor_distribution_batch <- TFact_RegPsite_cor_batch %>%
  mutate(group1 = cut(funscoR, 5, dig.lab = 1)) %>%
  #filter(abs(pearson_r) > 0.4) %>%
  group_by(batch, group1) %>%
  summarise(sig = sum(p.value < 0.05), nonsig = sum(p.value > 0.05)) %>%
  ungroup() %>%
  pivot_longer(-c(group1, batch), names_to = "signf", values_to = "count") %>%
  group_by(batch, group1) %>%
  mutate(perc = count/sum(count)) %>%
  ungroup()

counts <- TFact_RegPsite_cor_distribution_batch %>%
  group_by(batch, group1) %>%
  summarise(count = sum(count)) %>%
  ungroup()

TFact_RegPsite_cor_distribution_batch <- TFact_RegPsite_cor_distribution_batch %>%
  ggplot(mapping = aes(x = group1, y = perc, fill = signf)) +
  geom_col(color = "white", width = 1) +
  geom_text(data = counts, mapping = aes(x = group1, y = 0.9, label = count), size = 3, inherit.aes = F) +
  facet_wrap(~ batch, ncol = 5) +
  coord_flip() +
  theme_classic() +
  theme(
    panel.spacing = unit(0.4, "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size = 7),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "bottom") +
  scale_fill_manual(values = c("nonsig" = "#deebf7", "sig" = "#3182bd"), labels = c("nonsig" = "P > 5%", "sig" = "P < 5%")) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%")) +
  #scale_y_continuous(labels = scales::percent) +
  labs(x = "Functional score (binned)", y = "Percentage", fill = "Pearson's r") +
  guides(fill = guide_legend(nrow = 1))

ggsave(filename = "tf_phosphosite_correlation_barplot_funscore_batch.pdf", plot = TFact_RegPsite_cor_distribution_batch, path = "./output/plots/tf_reg_psites_correlation/", height = 4, width = 7)
ggsave(filename = "tf_phosphosite_correlation_barplot_funscore_batch.png", plot = TFact_RegPsite_cor_distribution_batch, path = "./output/plots/tf_reg_psites_correlation/", height = 4, width = 7)


TFact_RegPsite_cor_distribution <- all_cor %>%
  left_join(psp_reg, by = c("TFA" = "gene", "residue", "position")) %>%
  mutate(reg_status = replace_na(reg_status, "unknown")) %>%
  filter(reg_status == "activation") %>%
  mutate(class1 = fct_recode(class1, `Auto-psite` = "TF-Phosphosite", Others = "All pairs")) %>%
  mutate(class1 = fct_relevel(class1, "Others", "Auto-psite")) %>%
  ggplot(mapping = aes(x = class1, y = pearson_r)) +
  geom_boxplot(notch = T, color = "black", fill = "#f0f0f0", show.legend = F) +
  geom_jitter(mapping = aes(color = class1), width = 0.1, alpha = 0.4, show.legend = F) +
  stat_compare_means(method = "wilcox.test") +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 18),
    axis.text = element_text(colour = "black", size = 16)) +
  labs(x = "TF-phosphosite pair", y = "Pearson's r\n(TF activity vs TF activating phosphosites)")

ggsave(filename = "tf_phosphosite_correlation_auto_all_pairs_box.png", plot = TFact_RegPsite_cor_distribution, path = "./output/plots/tf_reg_psites_correlation/", height = 6, width = 4)
ggsave(filename = "tf_phosphosite_correlation_auto_all_pairs_box.pdf", plot = TFact_RegPsite_cor_distribution, path = "./output/plots/tf_reg_psites_correlation/", height = 6, width = 4)
