#' ---
#' title: "Comparison of phosphorylation data rna/protein/not-regressed out"
#' author: "abelsousa"
#' ---

#+ setup, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)


#' load R packages
library(tidyverse)

source("./src/utils/getSamples.R")


#' load all samples gathered in this study
all_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")


#' get samples with phosphorylation, protein and mRNA data
samples <- getSamples(all_samples, c("phosphorylation", "mRNA", "protein")) %>%
  select(sample)


# load quantile normalized phosphorylation data (log2FC)
# remove some phosphosites with NAs in all samples
phosphoQ <- read_tsv(file = "./output/files/phosphoproteomicsQ.txt.gz") %>%
  pivot_longer(names_to = "sample", values_to = "phos_log2fc", -c(gene, psite, psites)) %>%
  select(sample, gene, psite, psites, phos_log2fc) %>%
  group_by(gene, psite, psites) %>%
  filter( !(sum(is.na(phos_log2fc)) == n()) ) %>%
  ungroup() %>%
  pivot_wider(names_from = "sample", values_from = "phos_log2fc")

#' load rna regressed-out quantile normalized phosphorylation data (log2FC)
phosphoQrna_out <- read_tsv(file = "./output/files/phosphoproteomicsQ_RNAreg.txt.gz")

#' load protein regressed-out quantile normalized phosphorylation data (log2FC)
phosphoQprot_out <- read_tsv(file = "./output/files/phosphoproteomicsQ_Protreg.txt.gz")


measurements <- tibble(
  data = c("rna_reg", "prot_reg", "not_reg"),
  n = c(
    sum(!is.na(phosphoQrna_out[,-c(1:3)])),
    sum(!is.na(phosphoQprot_out[,-c(1:3)])),
    sum(!is.na(phosphoQ[,-c(1:3)])))) %>%
  mutate(data = fct_reorder(data, n)) %>%
  ggplot(mapping = aes(x = data, y = n)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 10)) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  labs(x = "Phospho data", y = "Number of measurements")


#+ fig.width=5, fig.height=2
measurements

ggsave(filename="phosphoQ_reg_comparison_measurements.png", plot = measurements, path = "./output/plots/data_description/", width=5, height=2)
unlink("phosphoQ_reg_comparison_measurements.png")


# load protein data (log2FC)
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(names_to = "sample", values_to = "prot_log2fc", -gene) %>%
  select(sample, protein=gene, prot_log2fc)

# tidy phosphoQ data
phosphoQ <- phosphoQ %>%
  pivot_longer(names_to = "sample", values_to = "phos_log2fc", -c(gene, psite, psites)) %>%
  select(sample, protein=gene, psite, psites, phos_log2fc)

# join phosphorylation and protein data
phosphoQ <- phosphoQ %>%
  inner_join(protein, by = c("sample", "protein")) %>%
  select(sample, prot=protein, psite, psites, prot_fc=prot_log2fc, psite_fc=phos_log2fc)


phosphoQ %>% filter(!is.na(psite_fc)) %>% nrow()
# 7,980,072

phosphoQ %>% filter(!is.na(psite_fc) & is.na(prot_fc)) %>% nrow()
# 661,034


measurements <- tibble(
  data = c("psite_not_NA", "psite_not_NA_&_prot_NA"),
  n = c(
    phosphoQ %>% filter(!is.na(psite_fc)) %>% nrow(),
    phosphoQ %>% filter(!is.na(psite_fc) & is.na(prot_fc)) %>% nrow())) %>%
  mutate(data = fct_reorder(data, n)) %>%
  ggplot(mapping = aes(x = data, y = n)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  labs(x = "Psite type", y = "Number of measurements")


#+ fig.width=5, fig.height=1
measurements

ggsave(filename="psites_number_NAs.png", plot = measurements, path = "./output/plots/data_description/", width=5, height=1)
unlink("psites_number_NAs.png")

