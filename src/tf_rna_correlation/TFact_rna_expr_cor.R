#' ---
#' title: "Correlation of TF activity to RNA expression"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, tidy.opts=list(width.cutoff=40), tidy=TRUE)


#' load R packages
library(tidyverse)
library(ggpubr)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load CPTAC
source("./src/utils/getSamples.R")
cptac_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")
cptac_samples <- getSamples(cptac_samples, c("protein","mRNA")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue)


#' load TF activity estimates
tf <- read_csv(file = "./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  rename(tf=X1) %>%
  pivot_longer(-tf, names_to = "sample", values_to = "tf_activity")


#' load rna abundance
rna <- read_tsv(file = "./output/files/transcriptomics_log2fc.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "log2fc")


#' # scatterplot - correlation of TF activity to RNA expression
tf_rna_expr <- tf %>%
  inner_join(rna, by = c("tf" = "gene", "sample")) %>%
  filter(!is.na(log2fc)) %>%
  ggplot(mapping = aes(x = tf_activity, y = log2fc)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  theme_classic() +
  stat_cor(size = 3) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10),
    aspect.ratio = 1) +
  labs(x = "TF activity estimate", y = "mRNA expression (log2fc)")


#+ fig.width=4, fig.height=4
tf_rna_expr

ggsave(filename = "TFact_rna_expr_scatterplot.png", plot = tf_rna_expr, path = "./output/plots/tf_rna_correlation/", height = 4, width = 4)
unlink("TFact_rna_expr_scatterplot.png")


#' # scatterplot - correlation of TF activity to RNA expression by batch
tf_rna_expr <- tf %>%
  inner_join(rna, by = c("tf" = "gene", "sample")) %>%
  filter(!is.na(log2fc)) %>%
  inner_join(cptac_samples, by = "sample") %>%
  ggplot(mapping = aes(x = tf_activity, y = log2fc)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~ batch, scales = "free") +
  theme_classic() +
  stat_cor(size = 3) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10),
    aspect.ratio = 1) +
  labs(x = "TF activity estimate", y = "mRNA expression (log2fc)")


#+ fig.width=6, fig.height=6
tf_rna_expr

ggsave(filename = "TFact_rna_expr_scatterplot_batch.png", plot = tf_rna_expr, path = "./output/plots/tf_rna_correlation/", height = 6, width = 6)
unlink("TFact_rna_expr_scatterplot_batch.png")


#' # boxplot - distribution of TF activity to RNA expression correlation across samples
tf_rna_expr <- tf %>%
  inner_join(rna, by = c("tf" = "gene", "sample")) %>%
  filter(!is.na(log2fc)) %>%
  group_by(tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$tf_activity, .x$log2fc)))) %>%
  select(-data) %>%
  unnest(cor) %>%
  ggplot(mapping = aes(x = factor(0), y = estimate)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "", y = "Pearson's r (TF activity vs RNA expression)")


#+ fig.width=4, fig.height=4
tf_rna_expr

ggsave(filename = "TFact_rna_expr_cor_boxpl.png", plot = tf_rna_expr, path = "./output/plots/tf_rna_correlation/", height = 4, width = 4)
unlink("TFact_rna_expr_cor_boxpl.png")


#' # boxplot - distribution of TF activity to RNA expression correlation across samples by batch
tf_rna_expr <- tf %>%
  inner_join(rna, by = c("tf" = "gene", "sample")) %>%
  filter(!is.na(log2fc)) %>%
  inner_join(cptac_samples, by = "sample") %>%
  group_by(batch, tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$tf_activity, .x$log2fc)))) %>%
  select(-data) %>%
  unnest(cor) %>%
  ggplot(mapping = aes(x = batch, y = estimate, fill = batch)) +
  geom_boxplot(show.legend = F) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "", y = "Pearson's r (TF activity vs RNA expression)")


#+ fig.width=4, fig.height=4
tf_rna_expr

ggsave(filename = "TFact_rna_expr_cor_boxpl_batch.png", plot = tf_rna_expr, path = "./output/plots/tf_rna_correlation/", height = 4, width = 4)
unlink("TFact_rna_expr_cor_boxpl_batch.png")

