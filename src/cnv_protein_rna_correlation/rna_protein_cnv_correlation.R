#' ---
#' title: "CNV to protein/rna correlation"
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



#' load CNV data
cnv <- read_tsv(file = "./output/files/cnv.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "cnv")


#' load protein data
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "prot_log2fc") %>%
  filter(!is.na(prot_log2fc))


#' load RNA data
rna <- read_tsv(file = "./output/files/transcriptomics_log2fc.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "rna_log2fc") %>%
  filter(!is.na(rna_log2fc))


#' load metadata
samples_metadata <- read_tsv("./output/files/all_samples_annotation.txt")

cnv_rna_prot_samples <- getSamples(samples_metadata, data_types = c("cnv", "mRNA", "protein")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue)


#' remove possible batch effects from mRNA and protein abundance

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
  inner_join(cnv_rna_prot_samples, by = "sample") %>%
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
  inner_join(cnv_rna_prot_samples, by = "sample") %>%
  group_by(gene) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(data2 = map(data, get_residuals)) %>%
  select(-data) %>%
  unnest() %>%
  select(-batch) %>%
  rename(rna_log2fc_resd = residuals)

multi_data <- cnv %>%
  inner_join(protein, by = c("gene", "sample")) %>%
  inner_join(rna, by = c("gene", "sample")) %>%
  inner_join(cnv_rna_prot_samples, by = "sample")


cnv_protein_plot1 <- multi_data %>%
  #slice(1:10000) %>%
  mutate(cnv2 = as.character(cnv)) %>%
  mutate(cnv2 = fct_relevel(cnv2, "-2", "-1", "0", "1", "2")) %>%
  ggplot(mapping = aes(x = cnv, y = prot_log2fc_resd)) +
  geom_boxplot(mapping = aes(fill = cnv, group = cnv2), outlier.shape = NA, show.legend = F) +
  geom_smooth(method = "lm", show.legend = F) +
  stat_cor(size = 5, label.y.npc = 1) +
  scale_y_continuous(limits = c(-2.5,2.5)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 16, colour = "black", hjust = 0.5)) +
  scale_fill_gradient2(low = "#a6611a", mid = "#f5f5f5", high = "#018571") +
  labs(x = "CNV (GISTIC2)", y = "Protein abundance (log2 fold-change)", title = "Protein vs CNV")

#+ fig.width=3, fig.height=5
cnv_protein_plot1

ggsave(filename = "protein_cnv_distribution_all_studies.png", plot = cnv_protein_plot1, path = "./output/plots/cnv_protein_rna_correlation/", height = 5, width = 3)
ggsave(filename = "protein_cnv_distribution_all_studies.pdf", plot = cnv_protein_plot1, path = "./output/plots/cnv_protein_rna_correlation/", height = 5, width = 3)
unlink("protein_cnv_distribution_all_studies.png")
unlink("protein_cnv_distribution_all_studies.pdf")


cnv_protein_plot2 <- multi_data %>%
  #slice(1:10000) %>%
  mutate(cnv2 = as.character(cnv)) %>%
  mutate(cnv2 = fct_relevel(cnv2, "-2", "-1", "0", "1", "2")) %>%
  mutate(batch = fct_recode(batch, Brain = "cbttc", `Breast (ccle)` = "ccle-breast", `Colorectal (ccle)` = "ccle-colorectal", `Colon`= "colon-opportunities", `Kidney` = "discovery-ccrcc", `Uterus` = "discovery-ucec", `Liver` = "hcc-proteogenomics", `Breast (tcga)` = "tcga-brca", `Colorectal (tcga)` = "tcga-coread", `Ovary` = "tcga-ov")) %>%
  ggplot(mapping = aes(x = cnv, y = prot_log2fc_resd)) +
  geom_boxplot(mapping = aes(fill = cnv, group = cnv2), outlier.shape = NA, show.legend = F) +
  geom_smooth(method = "lm", show.legend = F) +
  stat_cor(aes(label = ..r.label..), size = 4, label.y.npc = 1, label.x.npc = 0) +
  facet_wrap(~ batch, nrow = 2) +
  scale_y_continuous(limits = c(-2,2)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    plot.title = element_text(size = 16, colour = "black", hjust = 0.5)) +
  scale_fill_gradient2(low = "#a6611a", mid = "#f5f5f5", high = "#018571") +
  labs(x = "CNV (GISTIC2)", y = "Protein abundance (log2 fold-change)", title = "Protein vs CNV by study")

#+ fig.width=7, fig.height=5
cnv_protein_plot2

ggsave(filename = "protein_cnv_distribution_by_study.png", plot = cnv_protein_plot2, path = "./output/plots/cnv_protein_rna_correlation/", height = 5, width = 7)
ggsave(filename = "protein_cnv_distribution_by_study.pdf", plot = cnv_protein_plot2, path = "./output/plots/cnv_protein_rna_correlation/", height = 5, width = 7)
unlink("protein_cnv_distribution_by_study.png")
unlink("protein_cnv_distribution_by_study.pdf")


cnv_rna_plot1 <- multi_data %>%
  #slice(1:10000) %>%
  mutate(cnv2 = as.character(cnv)) %>%
  mutate(cnv2 = fct_relevel(cnv2, "-2", "-1", "0", "1", "2")) %>%
  ggplot(mapping = aes(x = cnv, y = rna_log2fc_resd)) +
  geom_boxplot(mapping = aes(fill = cnv, group = cnv2), outlier.shape = NA, show.legend = F) +
  geom_smooth(method = "lm", show.legend = F) +
  stat_cor(size = 5, label.y.npc = 1) +
  scale_y_continuous(limits = c(-2.5,2.5)) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 16, colour = "black", hjust = 0.5)) +
  scale_fill_gradient2(low = "#ca0020", mid = "#f7f7f7", high = "#0571b0") +
  labs(x = "CNV (GISTIC2)", y = "mRNA abundance (log2 fold-change)", title = "mRNA vs CNV")

#+ fig.width=3, fig.height=5
cnv_rna_plot1

ggsave(filename = "rna_cnv_distribution_all_studies.pdf", plot = cnv_rna_plot1, path = "./output/plots/cnv_protein_rna_correlation/", height = 5, width = 3)
ggsave(filename = "rna_cnv_distribution_all_studies.png", plot = cnv_rna_plot1, path = "./output/plots/cnv_protein_rna_correlation/", height = 5, width = 3)
unlink("rna_cnv_distribution_all_studies.pdf")
unlink("rna_cnv_distribution_all_studies.png")


cnv_rna_plot2 <- multi_data %>%
  #slice(1:10000) %>%
  mutate(cnv2 = as.character(cnv)) %>%
  mutate(cnv2 = fct_relevel(cnv2, "-2", "-1", "0", "1", "2")) %>%
  mutate(batch = fct_recode(batch, Brain = "cbttc", `Breast (ccle)` = "ccle-breast", `Colorectal (ccle)` = "ccle-colorectal", `Colon`= "colon-opportunities", `Kidney` = "discovery-ccrcc", `Uterus` = "discovery-ucec", `Liver` = "hcc-proteogenomics", `Breast (tcga)` = "tcga-brca", `Colorectal (tcga)` = "tcga-coread", `Ovary` = "tcga-ov")) %>%
  ggplot(mapping = aes(x = cnv, y = rna_log2fc_resd)) +
  geom_boxplot(mapping = aes(fill = cnv, group = cnv2), outlier.shape = NA, show.legend = F) +
  geom_smooth(method = "lm", show.legend = F) +
  stat_cor(aes(label = ..r.label..), size = 4, label.y.npc = 1, label.x.npc = 0) +
  facet_wrap(~ batch, nrow = 2) +
  scale_y_continuous(limits = c(-2,2)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    plot.title = element_text(size = 16, colour = "black", hjust = 0.5)) +
  scale_fill_gradient2(low = "#ca0020", mid = "#f7f7f7", high = "#0571b0") +
  labs(x = "CNV (GISTIC2)", y = "mRNA abundance (log2 fold-change)", title = "mRNA vs CNV by study")

#+ fig.width=7, fig.height=5
cnv_rna_plot2

ggsave(filename = "rna_cnv_distribution_by_study.pdf", plot = cnv_rna_plot2, path = "./output/plots/cnv_protein_rna_correlation/", height = 5, width = 7)
ggsave(filename = "rna_cnv_distribution_by_study.png", plot = cnv_rna_plot2, path = "./output/plots/cnv_protein_rna_correlation/", height = 5, width = 7)
unlink("rna_cnv_distribution_by_study.pdf")
unlink("rna_cnv_distribution_by_study.png")


#CNV to mRNA/protein correlation

#set up a correlation function
cor2 <- function(df, x, y, method){
  a = df[[x]]
  b = df[[y]]
  cor <- broom::tidy(cor.test(a, b, method = method))
  
  res <- cor %>%
    select(estimate, p.value)
  
  res
}

cnv_prot_rna_cor <- multi_data %>%
  filter(!(is.na(prot_log2fc) | is.na(rna_log2fc))) %>%
  group_by(gene) %>%
  nest() %>%
  ungroup() %>%
  mutate(f = map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  filter(f) %>%
  select(-f) %>%
  mutate(cor_prot = map(.x = data, .f = cor2, x = "prot_log2fc_resd", y = "cnv", "pearson")) %>%
  mutate(cor_rna = map(.x = data, .f = cor2, x = "rna_log2fc_resd", y = "cnv", "pearson")) %>%
  select(-data) %>%
  unnest(cor_prot) %>%
  rename(corr_cnv_prot = estimate, pval_cnv_prot = p.value) %>%
  unnest(cor_rna) %>%
  rename(corr_cnv_rna = estimate, pval_cnv_rna = p.value)

cnv_prot_rna_cor_plot <- cnv_prot_rna_cor %>%
  select(-pval_cnv_prot, -pval_cnv_rna) %>%
  pivot_longer(-gene, names_to = "corr", values_to = "pearson_r") %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = corr, y = pearson_r, fill=corr), outlier.shape = NA) +
  scale_fill_viridis(discrete = T) +
  scale_y_continuous(limits = c(-0.25,0.75)) +
  labs(title = "All studies")

#+ fig.width=4, fig.height=4
cnv_prot_rna_cor_plot

ggsave(filename = "prot_rna_cnv_cor_bygene_all_studies.png", plot = cnv_prot_rna_cor_plot, path = "./output/plots/cnv_protein_rna_correlation/", height = 4, width = 4)
ggsave(filename = "prot_rna_cnv_cor_bygene_all_studies.pdf", plot = cnv_prot_rna_cor_plot, path = "./output/plots/cnv_protein_rna_correlation/", height = 4, width = 4)


pvalues <- cnv_prot_rna_cor %>%
  select(-pval_cnv_prot, -pval_cnv_rna) %>%
  pivot_longer(-gene, names_to = "corr", values_to = "pearson_r") %>%
  nest() %>%
  mutate(pvalue = map(.x = data, .f = ~ broom::tidy(wilcox.test(formula = pearson_r ~ corr, data = .x))))  %>%
  select(-data) %>%
  unnest() %>%
  select(p.value) %>%
  mutate(label = if_else(p.value < 2.2e-16, "2.2e-16", as.character(p.value))) %>%
  mutate(label = str_c("Wilcoxon, p < ", label))

cnv_prot_rna_cor_plot <- cnv_prot_rna_cor %>%
  select(-pval_cnv_prot, -pval_cnv_rna) %>%
  pivot_longer(-gene, names_to = "corr", values_to = "pearson_r") %>%
  ggplot(mapping = aes(x = pearson_r, color = corr, fill = corr)) +
  geom_histogram(mapping = aes(y = stat(density)), position = "identity", alpha = 0.8, bins = 30) +
  geom_density(alpha = 0.4) +
  geom_text(data = pvalues, mapping = aes(x = 0.5, y = 3.5, label = label), inherit.aes = F, size = 5) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16, colour = "black"),
    axis.title = element_text(size = 18, colour = "black"),
    legend.text = element_text(size = 16, colour = "black"),
    legend.title = element_text(size = 18, colour = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal") +
  scale_fill_viridis(discrete = T, label = c("corr_cnv_prot" = "CNV vs Protein", "corr_cnv_rna" = "CNV vs mRNA")) +
  scale_color_viridis(discrete = T, label = c("corr_cnv_prot" = "CNV vs Protein", "corr_cnv_rna" = "CNV vs mRNA")) +
  labs(x = "Pearson's r", y = "Density", fill = "Correlations", color = "Correlations")

#+ fig.width=6, fig.height=3
cnv_prot_rna_cor_plot

ggsave(filename = "prot_rna_cnv_cor_bygene_all_studies2.png", plot = cnv_prot_rna_cor_plot, path = "./output/plots/cnv_protein_rna_correlation/", height = 3, width = 6)
ggsave(filename = "prot_rna_cnv_cor_bygene_all_studies2.pdf", plot = cnv_prot_rna_cor_plot, path = "./output/plots/cnv_protein_rna_correlation/", height = 3, width = 6)


cnv_prot_rna_cor_batch <- multi_data %>%
  filter(!(is.na(prot_log2fc) | is.na(rna_log2fc))) %>%
  group_by(batch, gene) %>%
  nest() %>%
  ungroup() %>%
  mutate(f = map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  filter(f) %>%
  select(-f) %>%
  mutate(cor_prot = map(.x = data, .f = cor2, x = "prot_log2fc", y = "cnv", "pearson")) %>%
  mutate(cor_rna = map(.x = data, .f = cor2, x = "rna_log2fc", y = "cnv", "pearson")) %>%
  select(-data) %>%
  unnest(cor_prot) %>%
  rename(corr_cnv_prot = estimate, pval_cnv_prot = p.value) %>%
  unnest(cor_rna) %>%
  rename(corr_cnv_rna = estimate, pval_cnv_rna = p.value)

cnv_prot_rna_cor_batch_plot <- cnv_prot_rna_cor_batch %>%
  select(-pval_cnv_prot, -pval_cnv_rna) %>%
  pivot_longer(-c(gene, batch), names_to = "corr", values_to = "pearson_r") %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = corr, y = pearson_r, fill=corr), outlier.shape = NA) +
  facet_wrap(~ batch, nrow = 2) +
  scale_fill_viridis(discrete = T) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.text = element_text(size=7)) +
  labs(title = "By study")

#+ fig.width=6, fig.height=4
cnv_prot_rna_cor_batch_plot

ggsave(filename = "prot_rna_cnv_cor_bygene_by_study.png", plot = cnv_prot_rna_cor_batch_plot, path = "./output/plots/cnv_protein_rna_correlation/", height = 4, width = 6)
ggsave(filename = "prot_rna_cnv_cor_bygene_by_study.pdf", plot = cnv_prot_rna_cor_batch_plot, path = "./output/plots/cnv_protein_rna_correlation/", height = 4, width = 6)


pvalues <- cnv_prot_rna_cor_batch %>%
  select(-pval_cnv_prot, -pval_cnv_rna) %>%
  pivot_longer(-c(gene, batch), names_to = "corr", values_to = "pearson_r") %>%
  group_by(batch) %>%
  nest() %>%
  mutate(pvalue = map(.x = data, .f = ~ broom::tidy(wilcox.test(formula = pearson_r ~ corr, data = .x))))  %>%
  select(-data) %>%
  unnest() %>%
  select(batch, p.value) %>%
  mutate(label = if_else(p.value < 2.2e-16, "2.2e-16", as.character(p.value))) %>%
  mutate(label = str_c("Wilcoxon, p < ", label)) %>%
  mutate(batch = fct_recode(batch, Brain = "cbttc", `Breast (ccle)` = "ccle-breast", `Colorectal (ccle)` = "ccle-colorectal", `Colon`= "colon-opportunities", `Kidney` = "discovery-ccrcc", `Uterus` = "discovery-ucec", `Liver` = "hcc-proteogenomics", `Breast (tcga)` = "tcga-brca", `Colorectal (tcga)` = "tcga-coread", `Ovary` = "tcga-ov"))

N <- cnv_prot_rna_cor_batch %>%
  group_by(batch) %>%
  tally() %>%
  mutate(batch = fct_recode(batch, Brain = "cbttc", `Breast (ccle)` = "ccle-breast", `Colorectal (ccle)` = "ccle-colorectal", `Colon`= "colon-opportunities", `Kidney` = "discovery-ccrcc", `Uterus` = "discovery-ucec", `Liver` = "hcc-proteogenomics", `Breast (tcga)` = "tcga-brca", `Colorectal (tcga)` = "tcga-coread", `Ovary` = "tcga-ov")) %>%
  mutate(batch = factor(x = as.character(batch), levels = sort(unique(as.character(batch)))))

cnv_prot_rna_cor_batch_plot <- cnv_prot_rna_cor_batch %>%
  select(-pval_cnv_prot, -pval_cnv_rna) %>%
  pivot_longer(-c(gene, batch), names_to = "corr", values_to = "pearson_r") %>%
  mutate(batch = fct_recode(batch, Brain = "cbttc", `Breast (ccle)` = "ccle-breast", `Colorectal (ccle)` = "ccle-colorectal", `Colon`= "colon-opportunities", `Kidney` = "discovery-ccrcc", `Uterus` = "discovery-ucec", `Liver` = "hcc-proteogenomics", `Breast (tcga)` = "tcga-brca", `Colorectal (tcga)` = "tcga-coread", `Ovary` = "tcga-ov")) %>%
  mutate(batch = factor(x = as.character(batch), levels = sort(unique(as.character(batch))))) %>%
  ggplot(mapping = aes(x = pearson_r, color = corr, fill = corr)) +
  geom_histogram(mapping = aes(y = stat(density)), position = "identity", alpha = 0.8, bins = 30, size = 0.1) +
  geom_density(alpha = 0.4, size = 0.2) +
  geom_text(data = N, mapping = aes(x = -0.5, y = 3, label = n), inherit.aes = F, size = 4) +
  facet_wrap(~ batch, nrow = 2, scales = "fixed") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 16, colour = "black"),
    legend.position = "none") +
  scale_fill_viridis(discrete = T, label = c("corr_cnv_prot" = "CNV vs Protein", "corr_cnv_rna" = "CNV vs mRNA")) +
  scale_color_viridis(discrete = T, label = c("corr_cnv_prot" = "CNV vs Protein", "corr_cnv_rna" = "CNV vs mRNA")) +
  labs(x = "Pearson's r", y = "Density", fill = "Correlations", color = "Correlations")

#+ fig.width=8, fig.height=3
cnv_prot_rna_cor_batch_plot

ggsave(filename = "prot_rna_cnv_cor_bygene_by_study2.png", plot = cnv_prot_rna_cor_batch_plot, path = "./output/plots/cnv_protein_rna_correlation/", height = 3, width = 8)
ggsave(filename = "prot_rna_cnv_cor_bygene_by_study2.pdf", plot = cnv_prot_rna_cor_batch_plot, path = "./output/plots/cnv_protein_rna_correlation/", height = 3, width = 8)
