#' ---
#' title: "Correlation between phosphorylation quantifications and protein abundance"
#' author: "abelsousa"
#' ---

#+ setup, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)



#' load R packages
library(RColorBrewer)
library(tidyverse)
library(viridis)

source("./src/utils/getSamples.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load all samples with phosphorylation data gathered in this study
phospho_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
  getSamples("phosphorylation") %>%
  mutate(batch = map2_chr(batch, tissue, ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(sample, batch) %>%
  mutate(batch = as.factor(batch))


#' load quantile normalized phosphorylation data (log2FC)\
#' remove some phosphosites with NAs in all samples
phospho1 <- read_tsv(file = "./output/files/phosphoproteomicsQ.txt.gz") %>%
  pivot_longer(names_to = "sample", values_to = "phos_log2fc", -c(gene, psite, psites)) %>%
  select(sample, protein=gene, psite, psites, phos_log2fc) %>%
  group_by(protein, psite, psites) %>%
  filter( !(sum(is.na(phos_log2fc)) == n()) ) %>%
  ungroup()


#' load protein regressed-out quantile normalized phosphorylation data (log2FC)\
#' remove some phosphosites with NAs in all samples
phospho2 <- read_tsv(file = "./output/files/phosphoproteomicsQ_Protreg_allSamp_withZtransf.txt.gz") %>%
  pivot_longer(names_to = "sample", values_to = "phos_log2fc", -c(gene, psite, psites)) %>%
  select(sample, protein=gene, psite, psites, phos_log2fc) %>%
  group_by(protein, psite, psites) %>%
  filter( !(sum(is.na(phos_log2fc)) == n()) ) %>%
  ungroup()


#' load protein abundance data (log2FC)
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(names_to = "sample", values_to = "prot_log2fc", -gene) %>%
  select(sample, protein=gene, prot_log2fc)


#' regress-out batch from protein/phosphosite abundance?
# source("./src/utils/regress_outCovs2.R")
# protein_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
#   getSamples("protein") %>%
#   mutate(batch = map2_chr(batch, tissue, ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
#   select(sample, batch) %>%
#   mutate(batch = as.factor(batch))
# 
# protein <- protein %>%
#   inner_join(protein_samples, by = "sample") %>%
#   select(protein, sample, prot_log2fc, batch) %>%
#   group_by(protein) %>%
#   nest() %>%
#   ungroup() %>%
#   mutate(residuals = map(.x = data, .f = regress_outCovs2)) %>%
#   unnest() %>%
#   select(-prot_log2fc, -batch) %>%
#   rename(prot_log2fc = residual)
# 
# phospho1 <- phospho1 %>%
#   inner_join(phospho_samples, by = "sample") %>%
#   select(psite, psites, protein, sample, phos_log2fc, batch) %>%
#   group_by(psite, psites, protein) %>%
#   nest() %>%
#   ungroup() %>%
#   mutate(residuals = map(.x = data, .f = regress_outCovs2)) %>%
#   unnest() %>%
#   select(-phos_log2fc, -batch) %>%
#   rename(phos_log2fc = residual)
# 
# phospho2 <- phospho2 %>%
#   inner_join(phospho_samples, by = "sample") %>%
#   select(psite, psites, protein, sample, phos_log2fc, batch) %>%
#   group_by(psite, psites, protein) %>%
#   nest() %>%
#   ungroup() %>%
#   mutate(residuals = map(.x = data, .f = regress_outCovs2)) %>%
#   unnest() %>%
#   select(-phos_log2fc, -batch) %>%
#   rename(phos_log2fc = residual)


#' set up a correlation function
cor_test <- function(df, x, y, met = "pearson"){
  a <- df[[x]]
  b <- df[[y]]
  
  corr <- broom::tidy(cor.test(a, b, method = met))
  
  res <- corr %>%
    select(estimate, p.value) %>%
    rename_all(.funs = ~ str_c(met, .x, sep = "_"))
  
  return(res)
}


#' # phosphosite to protein correlation (all studies together)\
#' ### select phosphosites with expression in at least 1% of samples (n > ~9 samples)
pho_prot1 <- phospho1 %>%
  #filter(psite %in% unique(phospho1$psite)[1:1000]) %>%
  inner_join(protein, by = c("sample", "protein")) %>%
  group_by(protein, psite, psites) %>%
  filter(sum(!is.na(phos_log2fc) & !is.na(prot_log2fc)) > 0.01*n()) %>%
  #filter(sum(is.na(phos_log2fc) | is.na(prot_log2fc)) > 0.01*n()) %>%
  nest() %>%
  ungroup() %>%
  mutate(pearson = map(.x=data, .f = cor_test, x = "phos_log2fc", y = "prot_log2fc")) %>%
  #mutate(spearman = map(.x=data, .f = cor_test, x = "phos_log2fc", y = "prot_log2fc", met = "spearman")) %>%
  select(-data) %>%
  unnest() %>%
  mutate(class1 = "Phosphosite vs Protein", class2 = "Non-regressed out")

pho_prot2 <- phospho2 %>%
  #filter(psite %in% unique(phospho2$psite)[1:1000]) %>%
  inner_join(protein, by = c("sample", "protein")) %>%
  group_by(protein, psite, psites) %>%
  filter(sum(!is.na(phos_log2fc) & !is.na(prot_log2fc)) > 0.01*n()) %>%
  #filter(sum(is.na(phos_log2fc) | is.na(prot_log2fc)) > 0.01*n()) %>%
  nest() %>%
  ungroup() %>%
  mutate(pearson = map(.x=data, .f = cor_test, x = "phos_log2fc", y = "prot_log2fc")) %>%
  #mutate(spearman = map(.x=data, .f = cor_test, x = "phos_log2fc", y = "prot_log2fc", met = "spearman")) %>%
  select(-data) %>%
  unnest() %>%
  mutate(class1 = "Phosphosite vs Protein", class2 = "Regressed out")

common_psites <- intersect(pho_prot1$psite, pho_prot2$psite)

#pho_prot1 <- pho_prot1 %>%
#  filter(psite %in% common_psites)

#pho_prot2 <- pho_prot2 %>%
#  filter(psite %in% common_psites)

pho_prot <- bind_rows(pho_prot1, pho_prot2)

N <- pho_prot %>%
  group_by(class2) %>%
  tally() %>%
  ungroup() %>%
  filter(class2 == "Regressed out") %>%
  mutate(x = 1.5, y = -0.6)

pho_prot_plot <- pho_prot %>%
  ggplot(mapping = aes(x = class2, y = pearson_estimate, fill = class2)) +
  geom_boxplot(outlier.shape = NA, show.legend = F, alpha = 0.8, size = 1, fatten = 1.5) +
  geom_text(data = N, mapping = aes(x = x, y = y, label = n), color = "black", size = 5) +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14)) +
  #scale_y_continuous(labels = function(x) formatC(x, 2)) +
  scale_fill_manual(values = c("#bdbdbd", "#67a9cf")) +
  scale_x_discrete(labels = c("Non-regressed out" = "Non-regressed out", "Regressed out" = "Regressed out")) +
  labs(x = "Phosphosites", y = "Pearson's r (phosphosite vs protein)")

#+ fig.width=3, fig.height=4
#pho_prot_plot

ggsave(filename="phospho_protein_correlation.png", plot = pho_prot_plot, path = "./output/plots/phospho_ProtRNA_correlation/", width=4, height=4)
ggsave(filename="phospho_protein_correlation.pdf", plot = pho_prot_plot, path = "./output/plots/phospho_ProtRNA_correlation/", width=4, height=4)


#' # phosphosite to protein correlation by batch\
#' ### select phosphosites with expression in at least 10% of samples of each batch
pho_prot1 <- phospho1 %>%
  #filter(psite %in% unique(phospho1$psite)[1:1000]) %>%
  inner_join(phospho_samples, by = "sample") %>%
  inner_join(protein, by = c("sample", "protein")) %>%
  group_by(batch, protein, psite, psites) %>%
  filter(sum(!is.na(phos_log2fc) & !is.na(prot_log2fc)) > 0.1*n()) %>%
  #filter(sum(is.na(phos_log2fc) | is.na(prot_log2fc)) > 0.1*n()) %>%
  nest() %>%
  ungroup() %>%
  mutate(pearson = map(.x=data, .f = cor_test, x = "phos_log2fc", y = "prot_log2fc")) %>%
  #mutate(spearman = map(.x=data, .f = cor_test, x = "phos_log2fc", y = "prot_log2fc", met = "spearman")) %>%
  select(-data) %>%
  unnest() %>%
  mutate(class1 = "Phosphosite vs Protein", class2 = "Non-regressed out")

pho_prot2 <- phospho2 %>%
  #filter(psite %in% unique(phospho2$psite)[1:1000]) %>%
  inner_join(phospho_samples, by = "sample") %>%
  inner_join(protein, by = c("sample", "protein")) %>%
  group_by(batch, protein, psite, psites) %>%
  filter(sum(!is.na(phos_log2fc) & !is.na(prot_log2fc)) > 0.1*n()) %>%
  #filter(sum(is.na(phos_log2fc) | is.na(prot_log2fc)) > 0.1*n()) %>%
  nest() %>%
  ungroup() %>%
  mutate(pearson = map(.x=data, .f = cor_test, x = "phos_log2fc", y = "prot_log2fc")) %>%
  #mutate(spearman = map(.x=data, .f = cor_test, x = "phos_log2fc", y = "prot_log2fc", met = "spearman")) %>%
  select(-data) %>%
  unnest() %>%
  mutate(class1 = "Phosphosite vs Protein", class2 = "Regressed out")

common_psites <- intersect(pho_prot1$psite, pho_prot2$psite)

#pho_prot1 <- pho_prot1 %>%
#  filter(psite %in% common_psites)

#pho_prot2 <- pho_prot2 %>%
#  filter(psite %in% common_psites)

pho_prot <- bind_rows(pho_prot1, pho_prot2)

N <- pho_prot %>%
  group_by(batch, class2) %>%
  tally() %>%
  ungroup() %>%
  filter(class2 == "Regressed out") %>%
  mutate(x = 1.5, y = -0.9) %>%
  mutate(batch = as.character(batch)) %>%
  mutate(batch = fct_recode(batch, Brain = "cbttc", `Colorectal (ccle)` = "ccle-colorectal", `Colon`= "colon-opportunities", `Kidney` = "discovery-ccrcc", `Lung` = "discovery-luad", `Uterus` = "discovery-ucec", `Stomach` = "eogc-proteogenomics", `Liver` = "hcc-proteogenomics", `Breast (tcga)` = "tcga-brca", `Ovary` = "tcga-ov")) %>%
  mutate(batch = factor(x = as.character(batch), levels = sort(unique(as.character(batch)))))

pho_prot_plot <- pho_prot %>%
  #mutate(batch = fct_reorder(batch, .x=pearson_estimate, .desc = TRUE)) %>%
  mutate(batch = fct_recode(batch, Brain = "cbttc", `Colorectal (ccle)` = "ccle-colorectal", `Colon`= "colon-opportunities", `Kidney` = "discovery-ccrcc", `Lung` = "discovery-luad", `Uterus` = "discovery-ucec", `Stomach` = "eogc-proteogenomics", `Liver` = "hcc-proteogenomics", `Breast (tcga)` = "tcga-brca", `Ovary` = "tcga-ov")) %>%
  mutate(batch = factor(x = as.character(batch), levels = sort(unique(as.character(batch))))) %>%
  ggplot(mapping = aes(x = class2, y = pearson_estimate, fill = class2)) +
  geom_boxplot(outlier.shape = NA, size = 0.8, alpha = 0.8, fatten = 1) +
  geom_text(data = N, mapping = aes(x = x, y = y, label = n), color = "black", size = 3.5) +
  facet_wrap(~ batch, nrow=2, scales = "fixed") +
  theme_classic() +
  theme(
    #panel.spacing = unit(1, "cm"),
    legend.position = "bottom",
    legend.text = element_text(color = "black", size = 12),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 11),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14)) +
  #scale_x_discrete(labels = c("Non-regressed out" = "Non-regr. out", "Regressed out" = "Regr. out")) +
  scale_fill_manual(values = c("#bdbdbd", "#67a9cf")) +
  labs(x = "Phosphosites", y = "Pearson's r (phosphosite vs protein)")

#+ fig.width=6, fig.height=4
#pho_prot_plot

ggsave(filename="phospho_protein_correlation_batch.png", plot = pho_prot_plot, path = "./output/plots/phospho_ProtRNA_correlation/", width=6, height=4)
ggsave(filename="phospho_protein_correlation_batch.pdf", plot = pho_prot_plot, path = "./output/plots/phospho_ProtRNA_correlation/", width=6, height=4)
