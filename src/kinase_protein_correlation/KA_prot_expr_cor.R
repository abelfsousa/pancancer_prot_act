#' ---
#' title: "Correlation of kinase activity to protein expression"
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


# load CPTAC/CCLE samples
cptac_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
  filter(data == "phosphorylation") %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue)


#' load kinase-activity inference data\
#' (quantile-normalized protein regressed-out phosphorylation data)
ka <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= 3) %>%
  select(-source_type, -n)

# ka <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ.txt.gz") %>%
#   filter(source_type == "DB_text-mining") %>%
#   filter(n >= 3) %>%
#   select(-source_type, -n)

# ka <- data.table::fread(file = "./data/Danish/kinaseActMatImputed.tsv") %>%
#   as.tibble() %>%
#   rename(kinase=V1) %>%
#   pivot_longer(-kinase, names_to = "sample", values_to = "log10P")


#' load protein abundance
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "log2fc")


#' # scatterplot - correlation of kinase activity to protein expression
ka_prot_expr <- ka %>%
  inner_join(protein, by = c("kinase" = "gene", "sample")) %>%
  filter(!is.na(log2fc)) %>%
  ggplot(mapping = aes(x = log10P, y = log2fc)) +
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
  labs(x = "Kinase activity score (log10)", y = "Protein expression (log2fc)")

#+ fig.width=4, fig.height=4
ka_prot_expr

ggsave(filename = "KA_prot_expr_cor_scatter.png", plot = ka_prot_expr, path = "./output/plots/kinase_protein_correlation/", height = 4, width = 4)
unlink("KA_prot_expr_cor_scatter.png")


#' # scatterplot - correlation of kinase activity to protein expression by batch
ka_prot_expr <- ka %>%
  inner_join(protein, by = c("kinase" = "gene", "sample")) %>%
  filter(!is.na(log2fc)) %>%
  inner_join(cptac_samples, by = "sample") %>%
  ggplot(mapping = aes(x = log10P, y = log2fc)) +
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
  labs(x = "Kinase activity score (log10)", y = "Protein expression (log2fc)")

#+ fig.width=6, fig.height=6
ka_prot_expr

ggsave(filename = "KA_prot_expr_cor_scatter_batch.png", plot = ka_prot_expr, path = "./output/plots/kinase_protein_correlation/", height = 6, width = 6)
unlink("KA_prot_expr_cor_scatter_batch.png")


#' # boxplot - distribution of kinase activity to protein expression correlation across samples
ka_prot_expr <- ka %>%
  inner_join(protein, by = c("kinase" = "gene", "sample")) %>%
  filter(!is.na(log2fc)) %>%
  group_by(kinase) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > 50) %>%
  select(-n) %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$log10P, .x$log2fc)))) %>%
  select(-data) %>%
  unnest(cor) %>%
  ggplot(mapping = aes(x = factor(0), y = estimate)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "", y = "Pearson's r")

#+ fig.width=4, fig.height=4
ka_prot_expr

ggsave(filename = "KA_prot_expr_cor_boxpl.png", plot = ka_prot_expr, path = "./output/plots/kinase_protein_correlation/", height = 4, width = 4)
unlink("KA_prot_expr_cor_boxpl.png")



#' # boxplot - distribution of kinase activity to protein expression correlation across samples by batch
ka_prot_expr <- ka %>%
  inner_join(protein, by = c("kinase" = "gene", "sample")) %>%
  filter(!is.na(log2fc)) %>%
  inner_join(cptac_samples, by = "sample") %>%
  group_by(batch, kinase) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > 10) %>%
  select(-n) %>%
  group_by(batch, kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$log10P, .x$log2fc)))) %>%
  select(-data) %>%
  unnest(cor) %>%
  ggplot(mapping = aes(x = batch, y = estimate, fill = batch)) +
  geom_boxplot(show.legend = F) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "", y = "Pearson's r")

#+ fig.width=4, fig.height=4
ka_prot_expr

ggsave(filename = "KA_prot_expr_cor_boxpl_batch.png", plot = ka_prot_expr, path = "./output/plots/kinase_protein_correlation/", height = 4, width = 4)
unlink("KA_prot_expr_cor_boxpl_batch.png")

