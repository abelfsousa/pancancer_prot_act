#' ---
#' title: "Correlation of kinase activity to copy-number variation"
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
cptac_samples <- getSamples(cptac_samples, c("cnv", "phosphorylation")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue)


#' load cnv data
cnv <- read_tsv(file = "./output/files/cnv.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "cnv")


#' load kinase-activity inference data\
#' (quantile-normalized protein regressed-out phosphorylation data)
ka <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz")

k_subN <- 3
k_sampN <- 0
ka_dbt <- ka %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  group_by(kinase) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > k_sampN) %>%
  select(-n)


#' KA distribution across cnv levels - all samples
ka_cnv <- ka_dbt %>%
  inner_join(cnv, by = c("kinase" = "gene", "sample")) %>%
  mutate(cnv = as.character(cnv)) %>%
  mutate(cnv = fct_relevel(.f=cnv, "-2", "-1", "0", "1", "2")) %>%
  ggplot(mapping = aes(x = cnv, y = log10P)) +
  geom_boxplot(mapping = aes(fill = cnv), notch = T, lwd=0.3, show.legend = F, outlier.size = 0.5) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "CNV", y = "Kinase activity (log10P)")

#+ fig.width=3, fig.height=3
ka_cnv

ggsave(filename = "KA_cnv_distribution.png", plot = ka_cnv, path = "./output/plots/cnv_impact_kin_activity/", height = 3, width = 3)
unlink("KA_cnv_distribution.png")


#' KA distribution across cnv levels - by batch
ka_cnv <- ka_dbt %>%
  inner_join(cnv, by = c("kinase" = "gene", "sample")) %>%
  mutate(cnv = as.character(cnv)) %>%
  mutate(cnv = fct_relevel(.f=cnv, "-2", "-1", "0", "1", "2")) %>%
  inner_join(cptac_samples, by = "sample") %>%
  ggplot(mapping = aes(x = cnv, y = log10P)) +
  geom_boxplot(mapping = aes(fill = cnv), notch = T, lwd=0.3, show.legend = F, outlier.size = 0.5) +
  facet_wrap(~ batch) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 8),
    axis.title = element_text(colour = "black", size = 8),
    axis.text = element_text(colour = "black", size = 6)) +
  labs(x = "CNV", y = "Kinase activity (log10P)")

#+ fig.width=5, fig.height=5
ka_cnv

ggsave(filename = "KA_cnv_distribution_batch.png", plot = ka_cnv, path = "./output/plots/cnv_impact_kin_activity/", height = 5, width = 5)
unlink("KA_cnv_distribution_batch.png")



#' KA to cnv correlation across samples - all samples
ka_cnv <- ka_dbt %>%
  inner_join(cnv, by = c("kinase" = "gene", "sample")) %>%
  group_by(kinase) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > 100) %>%
  select(-n) %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$log10P, .x$cnv)))) %>%
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
  labs(y = "KA vs CNV (Pearson's r)")

#+ fig.width=2, fig.height=4
ka_cnv

ggsave(filename = "KA_cnv_cor_distribution.png", plot = ka_cnv, path = "./output/plots/cnv_impact_kin_activity/", height = 4, width = 2)
unlink("KA_cnv_cor_distribution.png")


#' KA to cnv correlation across samples - by batch
ka_cnv <- ka_dbt %>%
  inner_join(cnv, by = c("kinase" = "gene", "sample")) %>%
  inner_join(cptac_samples, by = "sample") %>%
  group_by(batch, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$log10P, .x$cnv)))) %>%
  select(-data) %>%
  unnest(cor) %>%
  ggplot(mapping = aes(x = batch, y = estimate)) +
  geom_boxplot(mapping = aes(fill = batch), show.legend = F) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "Batch", y = "KA vs CNV (Pearson's r)")

#+ fig.width=6, fig.height=4
ka_cnv

ggsave(filename = "KA_cnv_cor_distribution_batch.png", plot = ka_cnv, path = "./output/plots/cnv_impact_kin_activity/", height = 4, width = 6)
unlink("KA_cnv_cor_distribution_batch.png")
