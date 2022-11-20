#' ---
#' title: "RNA regressed-out phosphorylation data to protein/RNA correlation"
#' author: "abelsousa"
#' ---

#+ setup, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)



#' load R packages
library(RColorBrewer)
library(tidyverse)
library(viridis)

source("./src/utils/cor_functions.R")



#' load all samples with phosphorylation data gathered in this study
phosho_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
  filter(data == "phosphorylation") %>%
  select(-data)


#' load rna regressed-out quantile normalized phosphorylation data (log2FC)
phospho <- read_tsv(file = "./output/files/phosphoproteomicsQ_RNAreg.txt.gz") %>%
  pivot_longer(names_to = "sample", values_to = "phos_log2fc", -c(gene, psite, psites)) %>%
  select(sample, protein=gene, psite, psites, phos_log2fc)


#' load protein data (log2FC)
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(names_to = "sample", values_to = "prot_log2fc", -gene) %>%
  select(sample, protein=gene, prot_log2fc)


#' load mRNA data (log2FC)
rna <- read_tsv(file = "./output/files/transcriptomics_log2fc.txt.gz") %>%
  pivot_longer(names_to = "sample", values_to = "rna_log2fc", -gene) %>%
  select(sample, protein=gene, rna_log2fc)



#' join phosphorylation and protein data
pho_prot <- phospho %>%
  inner_join(phosho_samples[, c("sample", "batch")], by = "sample") %>%
  inner_join(protein, by = c("sample", "protein")) %>%
  select(batch, sample, prot=protein, psite, psites, prot_fc=prot_log2fc, psite_fc=phos_log2fc) %>%
  mutate(prot_data = "protein_level")


#' join phosphorylation and rna data
pho_rna <- phospho %>%
  inner_join(phosho_samples[, c("sample", "batch")], by = "sample") %>%
  inner_join(rna, by = c("sample", "protein")) %>%
  select(batch, sample, prot=protein, psite, psites, prot_fc=rna_log2fc, psite_fc=phos_log2fc) %>%
  mutate(prot_data = "rna_level")


#' bind all data together
phospho <- bind_rows(pho_prot, pho_rna)

rm(protein, rna, pho_prot, pho_rna)


#' set nest and unnest functions to the previous tidyr version
#' current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' # phosphosite to protein/rna correlation by batch\
#' ### select phosphosites with expression in at least 10% of samples of each batch
cor <- phospho %>%
  group_by(prot_data, batch, psite, psites, prot) %>%
  filter(sum(!is.na(prot_fc) & !is.na(psite_fc)) > 0.1*n()) %>%
  nest() %>%
  ungroup() %>%
  mutate(r = map(.x=data, .f = ~ cor_test(.x$prot_fc, .x$psite_fc, method = "pearson"))) %>%
  select(-data) %>%
  unnest() %>%
  select(prot_data, batch, psite, psites, pearson=pearson_estimate, r_p_value=pearson_p.value) %>%
  group_by(prot_data, batch) %>%
  mutate(r_p_adjust = p.adjust(r_p_value, method = "BH")) %>%
  ungroup() %>%
  select(-r_p_value)


plot_data <- bind_rows(
  cor %>% mutate(genes = "all"),
  cor %>% filter(r_p_adjust < 0.05) %>% mutate(genes = "FDR<0.05")) %>%
  mutate(batch = fct_reorder(batch, pearson, .fun = function(x){median(x, na.rm = T)}))


plot <- plot_data %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = batch, y = pearson, fill = genes), outlier.shape = NA, lwd=0.2) +
  geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.4) +
  facet_wrap(~ prot_data,
             labeller=labeller(prot_data = c("protein_level" = "corr(psite, protein)", "rna_level" = "corr(psite, mRNA)"))) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 12),
    legend.title = element_text(color = "black", size = 14),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 14)) +
  labs(x = "Batch", y = "Pearson's r", fill = "Genes")

#+ fig.width=8, fig.height=5
plot

ggsave(filename="PhosQ_rnaregout_ProtRNA_corr_batch.png", plot = plot, path = "./output/plots/phospho_ProtRNA_correlation/", width=8, height=5)
unlink("PhosQ_rnaregout_ProtRNA_corr_batch.png")



#' # phosphosite to protein/rna correlation (all studies together)\
#' ### select phosphosites with expression in at least 1% of samples (n > ~9 samples)
cor <- phospho %>%
  group_by(prot_data, psite, psites, prot) %>%
  filter(sum(!is.na(prot_fc) & !is.na(psite_fc)) > 0.01*n()) %>%
  nest() %>%
  ungroup() %>%
  mutate(r = map(.x=data, .f = ~ cor_test(.x$prot_fc, .x$psite_fc, method = "pearson"))) %>%
  select(-data) %>%
  unnest() %>%
  select(prot_data, psite, psites, pearson=pearson_estimate, r_p_value=pearson_p.value) %>%
  group_by(prot_data) %>%
  mutate(r_p_adjust = p.adjust(r_p_value, method = "BH")) %>%
  ungroup() %>%
  select(-r_p_value)


plot_data <- bind_rows(
  cor %>% mutate(genes = "all"),
  cor %>% filter(r_p_adjust < 0.05) %>% mutate(genes = "FDR<0.05"))


plot <- plot_data %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = factor(0), y = pearson, fill = genes), outlier.shape = NA, lwd=0.2) +
  geom_hline(yintercept = 0, linetype="dashed", color = "black", size=0.4) +
  geom_text(
    data = plot_data %>%
      group_by(prot_data, genes) %>%
      tally() %>%
      ungroup() %>%
      mutate(x = c(seq(0.5, 1, by = 0.5), 0.75)) %>%
      mutate(y = -0.5),
    mapping = aes(x = x, y = y, color = genes, label=n),
    show.legend = F, size=3, nudge_x = 0.2) +
  facet_wrap(~ prot_data,
             labeller=labeller(prot_data = c("protein_level" = "corr(psite, protein)", "rna_level" = "corr(psite, mRNA)"))) +
  theme_classic() +
  theme(
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(color = "black", size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.text = element_text(color = "black", size = 10),
    legend.title = element_text(color = "black", size = 12),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 10)) +
  labs(x = "", y = "Pearson's r", fill = "Genes")

#+ fig.width=5, fig.height=4
plot

ggsave(filename="PhosQ_rnaregout_ProtRNA_corr.png", plot = plot, path = "./output/plots/phospho_ProtRNA_correlation/", width=5, height=4)
unlink("PhosQ_rnaregout_ProtRNA_corr.png")
