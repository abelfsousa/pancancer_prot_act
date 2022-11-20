#' ---
#' title: "Protein/mRNA abundance to kinase/TF activity correlation"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=6, fig.height=6)


#' load R packages
library(tidyverse)
library(ggpubr)
library(viridis)

source(file = "./src/utils/getSamples.R")

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load protein data
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "prot_log2fc") %>%
  filter(!is.na(prot_log2fc))


#' load RNA data
rna <- read_tsv(file = "./output/files/transcriptomics_log2fc.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "rna_log2fc") %>%
  filter(!is.na(rna_log2fc))


#' load kinase-activity inference data\
#' (quantile-normalized protein regressed-out phosphorylation data)
ka <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz")

k_subN <- 3
k_sampN <- 0
kin_activity <- ka %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  group_by(kinase) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > k_sampN) %>%
  select(-n) %>%
  rename(activity = log10P)


#' load CPTAC TF activities
tf_activity <- read_csv(file = "./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "activity") %>%
  rename(tf=X1) %>%
  select(sample, tf, activity)


#' load metadata
cptac_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")
cptac_samples <- getSamples(cptac_samples, c("phosphorylation", "mRNA")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue)


#' remove possible batch effects from the protein abundance

# function to compute the residuals from a linear model
get_residuals <- function(df){
  dfM <- df %>%
    rename_at(.vars = 2, ~ "log2fc") %>%
    as.data.frame() %>%
    column_to_rownames(var = "sample")
  
  # fit the models
  if(length(unique(dfM$batch)) == 1){
    resd <- tibble(sample = rownames(dfM), residuals = dfM[["log2fc"]])
  } else {
    m <- lm(log2fc ~ batch, data = dfM)
    resd <- residuals(m)
    resd <- tibble(sample = names(resd), residuals = unname(resd))
  }
  
  res <- df %>%
    inner_join(resd, by = "sample")
  
  res
}

protein <- protein %>%
  inner_join(cptac_samples, by = "sample") %>%
  group_by(gene) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(data2 = map(data, get_residuals)) %>%
  select(-data) %>%
  unnest() %>%
  select(-batch) %>%
  rename(prot_log2fc_resd = residuals)

rna <- rna %>%
  inner_join(cptac_samples, by = "sample") %>%
  group_by(gene) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(data2 = map(data, get_residuals)) %>%
  select(-data) %>%
  unnest() %>%
  select(-batch) %>%
  rename(rna_log2fc_resd = residuals)


#' kinase/tf to protein/mRNA abundance correlation
# k <- kin_activity %>%
#   rename(id = kinase) %>%
#   mutate(protein_type = "kinase")
# 
# t <- tf_activity %>%
#   rename(id = tf) %>%
#   mutate(protein_type = "TF")
# 
# x <- bind_rows(k, t) %>%
#   inner_join(protein, by = c("sample", "id" = "gene")) %>%
#   inner_join(rna, by = c("sample", "id" = "gene")) %>%
#   inner_join(cptac_samples, by = "sample") %>%
#   group_by(protein_type, id) %>%
#   nest() %>%
#   ungroup() %>%
#   filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
#   mutate(cor1 = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$prot_log2fc_resd, method = "pearson")) %>% select(estimate))) %>%
#   mutate(cor2 = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$rna_log2fc_resd, method = "pearson")) %>% select(estimate))) %>%
#   select(-data) %>%
#   unnest(cor1) %>%
#   rename(protein = estimate) %>%
#   unnest(cor2) %>%
#   rename(rna = estimate) %>%
#   pivot_longer(-c(protein_type, id), names_to = "cor_type", values_to = "cor_value")


#' kinase to protein abundance correlation
kin_prot <- kin_activity %>%
  inner_join(protein, by = c("sample", "kinase" = "gene")) %>%
  inner_join(cptac_samples, by = "sample")

#' all samples
kin_prot_corr_all <- kin_prot %>%
  group_by(kinase) %>%
  nest() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$prot_log2fc_resd, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(id = kinase, cor_value = estimate) %>%
  mutate(protein_type = "kinase", cor_type = "protein")

#' by batch
kin_prot_corr_by_batch <- kin_prot %>%
  group_by(batch, kinase) %>%
  nest() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$prot_log2fc_resd, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(batch, id = kinase, cor_value = estimate) %>%
  mutate(protein_type = "kinase", cor_type = "protein")


#' kinase to mRNA expression correlation
kin_rna <- kin_activity %>%
  inner_join(rna, by = c("sample", "kinase" = "gene")) %>%
  inner_join(cptac_samples, by = "sample")

#' all samples
kin_rna_corr_all <- kin_rna %>%
  group_by(kinase) %>%
  nest() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$rna_log2fc_resd, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(id = kinase, cor_value = estimate) %>%
  mutate(protein_type = "kinase", cor_type = "rna")

#' by batch
kin_rna_corr_by_batch <- kin_rna %>%
  group_by(batch, kinase) %>%
  nest() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$rna_log2fc_resd, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(batch, id = kinase, cor_value = estimate) %>%
  mutate(protein_type = "kinase", cor_type = "rna")


#' TF to protein abundance correlation
tf_prot <- tf_activity %>%
  inner_join(protein, by = c("sample", "tf" = "gene")) %>%
  inner_join(cptac_samples, by = "sample")

#' all samples
tf_prot_corr_all <- tf_prot %>%
  group_by(tf) %>%
  nest() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$prot_log2fc_resd, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(id = tf, cor_value = estimate) %>%
  mutate(protein_type = "TF", cor_type = "protein")

#' by batch
tf_prot_corr_by_batch <- tf_prot %>%
  group_by(batch, tf) %>%
  nest() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$prot_log2fc_resd, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(batch, id = tf, cor_value = estimate) %>%
  mutate(protein_type = "TF", cor_type = "protein")


#' TF to mRNA expression correlation
tf_rna <- tf_activity %>%
  inner_join(rna, by = c("sample", "tf" = "gene")) %>%
  inner_join(cptac_samples, by = "sample")

#' all samples
tf_rna_corr_all <- tf_rna %>%
  group_by(tf) %>%
  nest() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$rna_log2fc_resd, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(id = tf, cor_value = estimate) %>%
  mutate(protein_type = "TF", cor_type = "rna")

#' by batch
tf_rna_corr_by_batch <- tf_rna %>%
  group_by(batch, tf) %>%
  nest() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$rna_log2fc_resd, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(batch, id = tf, cor_value = estimate) %>%
  mutate(protein_type = "TF", cor_type = "rna")


#' join kinase and TF correlations (all samples)
kin_tf_corr_all <- kin_prot_corr_all %>%
  bind_rows(kin_rna_corr_all) %>%
  bind_rows(tf_prot_corr_all) %>%
  bind_rows(tf_rna_corr_all)

#' plot correlation distributions
kin_tf_prot_cor_all_plot <- kin_tf_corr_all %>%
  filter(cor_type == "protein") %>%
  ggplot(mapping = aes(x = protein_type, y = cor_value, fill = protein_type)) +
  geom_boxplot(size = 2, outlier.size = 5, show.legend = F, alpha = 0.7) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 1),
    plot.title = element_text(size = 28, colour = "black", hjust = 0.5),
    axis.title.x = element_text(size = 28, colour = "black"),
    axis.text.x = element_text(size = 26, colour = "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    legend.position = "none") +
  scale_x_discrete(labels = c(kinase = "Kinases", TF = "TFs")) +
  scale_fill_viridis(discrete = T) +
  labs(y = "Pearson's r", title =  "Protein abundance vs protein activity")

#+ fig.width=14, fig.height=5
kin_tf_prot_cor_all_plot

ggsave(filename = "kinase_tf_prot_correlation_boxplot.png", plot = kin_tf_prot_cor_all_plot, path = "./output/plots/prot_rna_kinase_tf_act_correlation/", height = 5, width = 14)
ggsave(filename = "kinase_tf_prot_correlation_boxplot.pdf", plot = kin_tf_prot_cor_all_plot, path = "./output/plots/prot_rna_kinase_tf_act_correlation/", height = 5, width = 14)


#' join kinase and TF correlations (by batch)
kin_tf_corr_by_batch <- kin_prot_corr_by_batch %>%
  bind_rows(kin_rna_corr_by_batch) %>%
  bind_rows(tf_prot_corr_by_batch) %>%
  bind_rows(tf_rna_corr_by_batch)

#' plot correlation distributions
N <- kin_tf_corr_by_batch %>%
  filter(cor_type == "protein") %>%
  group_by(batch, protein_type) %>%
  tally() %>%
  ungroup() %>%
  mutate(batch = fct_recode(batch, Brain = "cbttc", `Colorectal (ccle)` = "ccle-colorectal", `Colon`= "colon-opportunities", `Kidney` = "discovery-ccrcc", `Uterus` = "discovery-ucec", `Lung` = "discovery-luad", `Stomach` = "eogc-proteogenomics", `Liver` = "hcc-proteogenomics", `Breast (tcga)` = "tcga-brca", `Ovary` = "tcga-ov")) %>%
  mutate(batch = factor(x = as.character(batch), levels = sort(unique(as.character(batch)))))

kin_tf_prot_cor_by_batch_plot <- kin_tf_corr_by_batch %>%
  filter(cor_type == "protein") %>%
  mutate(batch = fct_recode(batch, Brain = "cbttc", `Colorectal (ccle)` = "ccle-colorectal", `Colon`= "colon-opportunities", `Kidney` = "discovery-ccrcc", `Uterus` = "discovery-ucec", `Lung` = "discovery-luad", `Stomach` = "eogc-proteogenomics", `Liver` = "hcc-proteogenomics", `Breast (tcga)` = "tcga-brca", `Ovary` = "tcga-ov")) %>%
  mutate(batch = factor(x = as.character(batch), levels = sort(unique(as.character(batch))))) %>%
  ggplot(mapping = aes(x = protein_type, y = cor_value, fill = protein_type)) +
  geom_boxplot(size = 1, outlier.size = 2, show.legend = F, alpha = 0.7) +
  geom_text(data = N, mapping = aes(x = rep(c(1.3, 2.3),10), y = -0.6, label = n), size = 6) +
  theme_classic() +
  coord_flip() +
  facet_wrap(facets = vars(batch), scales = "fixed", nrow = 2) +
  theme(
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 1),
    plot.title = element_text(size = 28, colour = "black", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 24, colour = "black"),
    axis.title.x = element_text(size = 28, colour = "black"),
    axis.text.x = element_text(size = 26, colour = "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    legend.position = "none") +
  scale_x_discrete(labels = c(kinase = "Kinases", TF = "TFs")) +
  scale_fill_viridis(discrete = T) +
  labs(y = "Pearson's r", title =  "Protein abundance vs protein activity")

#+ fig.width=14, fig.height=5
kin_tf_prot_cor_by_batch_plot

ggsave(filename = "kinase_tf_prot_correlation_by_batch_boxplot.png", plot = kin_tf_prot_cor_by_batch_plot, path = "./output/plots/prot_rna_kinase_tf_act_correlation/", height = 5, width = 14)
ggsave(filename = "kinase_tf_prot_correlation_by_batch_boxplot.pdf", plot = kin_tf_prot_cor_by_batch_plot, path = "./output/plots/prot_rna_kinase_tf_act_correlation/", height = 5, width = 14)
