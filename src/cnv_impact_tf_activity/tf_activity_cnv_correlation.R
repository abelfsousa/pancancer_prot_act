#' ---
#' title: "Correlation of TF activity to copy-number variation"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=6, fig.height=6)


#' load R packages
library(tidyverse)
library(ggpubr)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load CPTAC samples
source("./src/utils/getSamples.R")
cptac_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")
cptac_samples <- getSamples(cptac_samples, c("protein","mRNA", "cnv")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue)


#' load cnv data
cnv <- read_tsv(file = "./output/files/cnv.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "cnv")


#' load CPTAC TF activities
tf_activity <- read_csv(file = "./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "tf_activity") %>%
  rename(tf=X1) %>%
  select(sample, tf, tf_activity)


#' TF activity distribution across CNV levels - all samples
tf_cnv <- tf_activity %>%
  inner_join(cnv, by = c("tf" = "gene", "sample")) %>%
  mutate(cnv = as.character(cnv)) %>%
  mutate(cnv = fct_relevel(.f=cnv, "-2", "-1", "0", "1", "2")) %>%
  ggplot(mapping = aes(x = cnv, y = tf_activity)) +
  geom_boxplot(mapping = aes(fill = cnv), notch = T, lwd=0.3, show.legend = F, outlier.size = 0.5) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "CNV", y = "TF activity")

#+ fig.width=3, fig.height=3
tf_cnv

ggsave(filename = "tf_activity_cnv_distribution.png", plot = tf_cnv, path = "./output/plots/cnv_impact_tf_activity/", height = 3, width = 3)
unlink("tf_activity_cnv_distribution.png")


#' TF activity distribution across CNV levels - by batch
tf_cnv <- tf_activity %>%
  inner_join(cnv, by = c("tf" = "gene", "sample")) %>%
  mutate(cnv = as.character(cnv)) %>%
  mutate(cnv = fct_relevel(.f=cnv, "-2", "-1", "0", "1", "2")) %>%
  inner_join(cptac_samples, by = "sample") %>%
  ggplot(mapping = aes(x = cnv, y = tf_activity)) +
  geom_boxplot(mapping = aes(fill = cnv), notch = T, lwd=0.3, show.legend = F, outlier.size = 0.5) +
  facet_wrap(~ batch) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 8),
    axis.title = element_text(colour = "black", size = 8),
    axis.text = element_text(colour = "black", size = 6)) +
  labs(x = "CNV", y = "TF activity")

#+ fig.width=5, fig.height=5
tf_cnv

ggsave(filename = "tf_activity_cnv_distribution_batch.png", plot = tf_cnv, path = "./output/plots/cnv_impact_tf_activity/", height = 5, width = 5)
unlink("tf_activity_cnv_distribution_batch.png")



#' TF activity to CNV correlation across samples - all samples
tf_cnv <- tf_activity %>%
  inner_join(cnv, by = c("tf" = "gene", "sample")) %>%
  group_by(tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$tf_activity, .x$cnv)))) %>%
  select(-data) %>%
  unnest(cor) %>%
  ggplot(mapping = aes(x = factor(0), y = estimate)) +
  geom_boxplot(fill = "grey") +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  labs(y = "TF activity vs CNV (Pearson's r)")

#+ fig.width=2, fig.height=4
tf_cnv

ggsave(filename = "tf_activity_cnv_cor_distribution.png", plot = tf_cnv, path = "./output/plots/cnv_impact_tf_activity/", height = 4, width = 2)
unlink("tf_activity_cnv_cor_distribution.png")


#' TF activity to CNV correlation across samples - by batch
tf_cnv <- tf_activity %>%
  inner_join(cnv, by = c("tf" = "gene", "sample")) %>%
  inner_join(cptac_samples, by = "sample") %>%
  group_by(batch, tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$tf_activity, .x$cnv)))) %>%
  select(-data) %>%
  unnest(cor) %>%
  ggplot(mapping = aes(x = batch, y = estimate)) +
  geom_boxplot(mapping = aes(fill = batch), show.legend = F) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "Batch", y = "TF activity vs CNV (Pearson's r)")

#+ fig.width=6, fig.height=4
tf_cnv

ggsave(filename = "tf_activity_cnv_cor_distribution_batch.png", plot = tf_cnv, path = "./output/plots/cnv_impact_tf_activity/", height = 4, width = 6)
unlink("tf_activity_cnv_cor_distribution_batch.png")
