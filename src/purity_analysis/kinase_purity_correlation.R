#' ---
#' title: "Correlation between kinase activity and purity scores"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=12, fig.height=6)


#' load R packages
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(ggrepel)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load samples annotation
source("./src/utils/getSamples.R")
samples_annotation <- read_tsv(file = "./output/files/all_samples_annotation.txt")
samples_annotation <- getSamples(samples_annotation, c("phosphorylation")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue)


#' load kinase activities and imputed values
kin_matrix <- data.table::fread("./data/Danish/kinaseActMatImputed.tsv") %>%
  as.data.frame() %>%
  column_to_rownames(var = "V1") %>%
  t()

kinases <- kin_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  as_tibble() %>%
  arrange(sample)


#' load purity scores
purity <- read_tsv("./output/files/samples_purity.txt") %>%
  select(-metric)


#' ## correlate each kinase with purity scores
kin_purity <- kinases %>%
  inner_join(purity, by = "sample")

cor_mat <- cor(kin_purity[, -c(1, length(kin_purity))], kin_purity[, c(length(kin_purity))], use = "pairwise.complete.obs")


#' ## heatmap of the correlations
heatmap <- Heatmap(
  matrix = cor_mat,
  name = "Pearson's r",
  cluster_columns = F,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8))

#+ fig.width=3, fig.height=10
heatmap

pdf(file="./output/plots/purity_analysis/kinases_purity_cor.pdf", height=10, width=3)
print(heatmap)
dev.off()

png(file="./output/plots/purity_analysis/kinases_purity_cor.png", units = "in", res = 300, height=10, width=3)
print(heatmap)
dev.off()

cors <- kin_purity %>%
  pivot_longer(-c(sample, score), names_to = "kinase", values_to = "activity") %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(data, ~ cor.test(.x$activity, .x$score) %>% broom::tidy() %>% select(estimate, p.value))) %>%
  unnest(cor) %>%
  arrange(desc(estimate))

#+ fig.width=6, fig.height=3
cors %>%
  filter(kinase == "PAK4") %>%
  unnest() %>%
  inner_join(samples_annotation, by = "sample") %>%
  ggplot(mapping = aes(x = score, y = activity)) +
  geom_point(mapping = aes(color = batch)) +
  ggpubr::stat_cor() +
  geom_smooth(method = "lm")
