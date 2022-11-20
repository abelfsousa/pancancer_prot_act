#' ---
#' title: "Gene-wise protein to mRNA correlation across samples"
#' author: "Abel Sousa"
#' date: "November 18th, 2019"
#' ---

#+ setup, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)


#' load R packages
library(RColorBrewer)
library(tidyverse)
library(cowplot)
library(viridis)


all_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")


#' get samples with protein and mRNA
protein_rna_samples <- all_samples %>%
  filter(data %in% c("protein", "mRNA")) %>%
  rename(info = data) %>%
  group_by(info, batch) %>%
  summarise(sample = list(sample)) %>%
  ungroup() %>%
  spread(key = "info", value = "sample") %>%
  group_by(batch) %>%
  mutate(overlap = list(reduce(compact(list(mRNA[[1]], protein[[1]])), intersect))) %>%
  ungroup() %>%
  select(batch, overlap) %>%
  unnest(cols = overlap) %>%
  rename(sample = overlap)
  
  

#' load protein data
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  gather(key = "sample", value = "prot_log2fc", -gene)


#' load rna data
rna_log2fc <- read_tsv(file = "./output/files/transcriptomics_log2fc.txt.gz") %>%
  gather(key = "sample", value = "rna_log2fc", -gene)

rna_log2fpkm <- read_tsv(file = "./output/files/transcriptomics_log2fpkm.txt.gz") %>%
  gather(key = "sample", value = "rna_log2fpkm", -gene)

rna <- rna_log2fc %>%
  inner_join(rna_log2fpkm, by = c("gene", "sample"))

rm(rna_log2fc, rna_log2fpkm)
gc()


#' add protein and rna data to samples
protein_rna_samples <- protein_rna_samples %>%
  inner_join(protein, by = "sample") %>%
  inner_join(rna, by = c("sample", "gene"))

rm(protein, rna)
gc()



#' # protein to mRNA (log2fc) correlation

#' ## by batch

#' pearson correlation
prot_rnaLog2FC_pearson <- protein_rna_samples %>%
  group_by(batch, gene) %>%
  filter(sum(!(is.na(prot_log2fc) | is.na(rna_log2fc))) > 0.2*n()) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor_test = map(.x=data, .f = ~ broom::tidy(cor.test(.x$prot_log2fc, .x$rna_log2fc, method = c("pearson"))))) %>%
  select(-data) %>%
  unnest(cols = cor_test) %>%
  select(batch, gene, r=estimate, p_value=p.value) %>%
  group_by(batch) %>%
  mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  ungroup()

#' spearman correlation
prot_rnaLog2FC_spearman <- protein_rna_samples %>%
  group_by(batch, gene) %>%
  filter(sum(!(is.na(prot_log2fc) | is.na(rna_log2fc))) > 0.2*n()) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor_test = map(.x=data, .f = ~ broom::tidy(cor.test(.x$prot_log2fc, .x$rna_log2fc, method = c("spearman"))))) %>%
  select(-data) %>%
  unnest(cols = cor_test) %>%
  select(batch, gene, p=estimate, p_value=p.value) %>%
  group_by(batch) %>%
  mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  ungroup()

plot_data <- prot_rnaLog2FC_pearson %>% 
  mutate(genes = "all", cor_type = "Pearson") %>%
  rename(cor = r) %>%
  bind_rows(prot_rnaLog2FC_spearman %>% mutate(genes = "all", cor_type = "Spearman") %>% rename(cor = p)) %>%
  bind_rows(prot_rnaLog2FC_pearson %>% filter(p_adjust<0.05) %>% rename(cor = r) %>% mutate(genes = "FDR<0.05", cor_type = "Pearson")) %>%
  bind_rows(prot_rnaLog2FC_spearman %>% filter(p_adjust<0.05) %>% rename(cor = p) %>% mutate(genes = "FDR<0.05", cor_type = "Spearman")) %>%
  mutate(batch = fct_reorder(batch, cor, .fun = function(x) {median(x, na.rm = T)}))

plot <- plot_data %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = batch, y = cor, fill = genes), outlier.shape = NA, lwd=0.2) +
  geom_text(
    data = plot_data %>%
      group_by(cor_type, batch, genes) %>%
      tally() %>%
      ungroup() %>%
      mutate(x = c(seq(0.5, 11, by = 0.5), seq(0.5, 11, by = 0.5))) %>%
      mutate(y = -0.6),
    mapping = aes(x = x, y = y, color = genes, label=n),
    show.legend = F, size=3, nudge_x = 0.2) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 12),
    legend.title = element_text(color = "black", size = 14),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 14)) +
  facet_wrap(~ cor_type) +
  labs(x = "Batch", y = "Correlation", fill = "Genes")

#+ fig.width=8, fig.height=5
plot

ggsave(filename="ProtRNAlog2fc_correlation_batch.png", plot = plot, path = "./output/plots/protein_rna_correlation/", width=8, height=5)
unlink("ProtRNAlog2fc_correlation_batch.png")


#' ## all samples

samples <- length(unique(protein_rna_samples$sample))

#' pearson correlation
prot_rnaLog2FC_pearson <- protein_rna_samples %>%
  filter(!(is.na(prot_log2fc) | is.na(rna_log2fc))) %>%
  group_by(gene) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > samples*0.2)) %>%
  mutate(cor_test = map(.x=data, .f = ~ broom::tidy(cor.test(.x$prot_log2fc, .x$rna_log2fc, method = c("pearson"))))) %>%
  select(-data) %>%
  unnest(cols = cor_test) %>%
  select(gene, r=estimate, p_value=p.value) %>%
  mutate(p_adjust = p.adjust(p_value, method = "BH"))

#' spearman correlation
prot_rnaLog2FC_spearman <- protein_rna_samples %>%
  filter(!(is.na(prot_log2fc) | is.na(rna_log2fc))) %>%
  group_by(gene) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > samples*0.2)) %>%
  mutate(cor_test = map(.x=data, .f = ~ broom::tidy(cor.test(.x$prot_log2fc, .x$rna_log2fc, method = c("spearman"))))) %>%
  select(-data) %>%
  unnest(cols = cor_test) %>%
  select(gene, p=estimate, p_value=p.value) %>%
  mutate(p_adjust = p.adjust(p_value, method = "BH"))

plot_data <- prot_rnaLog2FC_pearson %>% 
  mutate(genes = "all", cor_type = "Pearson") %>%
  rename(cor = r) %>%
  bind_rows(prot_rnaLog2FC_spearman %>% mutate(genes = "all", cor_type = "Spearman") %>% rename(cor = p)) %>%
  bind_rows(prot_rnaLog2FC_pearson %>% filter(p_adjust<0.05) %>% rename(cor = r) %>% mutate(genes = "FDR<0.05", cor_type = "Pearson")) %>%
  bind_rows(prot_rnaLog2FC_spearman %>% filter(p_adjust<0.05) %>% rename(cor = p) %>% mutate(genes = "FDR<0.05", cor_type = "Spearman"))

plot <- plot_data %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = genes, y = cor, fill = genes), outlier.shape = NA, lwd=0.2) +
  geom_text(
    data = plot_data %>%
      group_by(cor_type, genes) %>%
      tally() %>%
      ungroup() %>%
      mutate(x = c(0.8, 1.8, 0.8, 1.8)) %>%
      mutate(y = -0.2),
    mapping = aes(x = x, y = y, color = genes, label=n),
    show.legend = F, size=4, nudge_x = 0.2) +
  theme_classic() +
  theme(
    axis.text.y = element_text(color = "black", size = 12),
    axis.title.y = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 12),
    legend.title = element_text(color = "black", size = 14),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 14)) +
  facet_wrap(~ cor_type) +
  scale_fill_discrete(guide = F) +
  labs(x = "FDR", y = "Correlation", fill = "Genes")

#+ fig.width=3, fig.height=5
plot

ggsave(filename="ProtRNAlog2fc_correlation.png", plot = plot, path = "./output/plots/protein_rna_correlation/", width=3, height=5)
unlink("ProtRNAlog2fc_correlation.png")






#' # protein to mRNA (log2fpkm) correlation

#' pearson correlation
prot_rnaLog2FPKM_pearson <- protein_rna_samples %>%
  group_by(batch, gene) %>%
  filter(sum(!(is.na(prot_log2fc) | is.na(rna_log2fpkm))) > 0.2*n()) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor_test = map(.x=data, .f = ~ broom::tidy(cor.test(.x$prot_log2fc, .x$rna_log2fpkm, method = c("pearson"))))) %>%
  select(-data) %>%
  unnest(cols = cor_test) %>%
  select(batch, gene, r=estimate, p_value=p.value) %>%
  group_by(batch) %>%
  mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  ungroup()

#' spearman correlation
prot_rnaLog2FPKM_spearman <- protein_rna_samples %>%
  group_by(batch, gene) %>%
  filter(sum(!(is.na(prot_log2fc) | is.na(rna_log2fpkm))) > 0.2*n()) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor_test = map(.x=data, .f = ~ broom::tidy(cor.test(.x$prot_log2fc, .x$rna_log2fpkm, method = c("spearman"))))) %>%
  select(-data) %>%
  unnest(cols = cor_test) %>%
  select(batch, gene, p=estimate, p_value=p.value) %>%
  group_by(batch) %>%
  mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  ungroup()

plot_data <- prot_rnaLog2FPKM_pearson %>% 
  mutate(genes = "all", cor_type = "Pearson") %>%
  rename(cor = r) %>%
  bind_rows(prot_rnaLog2FPKM_spearman %>% mutate(genes = "all", cor_type = "Spearman") %>% rename(cor = p)) %>%
  bind_rows(prot_rnaLog2FPKM_pearson %>% filter(p_adjust<0.05) %>% rename(cor = r) %>% mutate(genes = "FDR<0.05", cor_type = "Pearson")) %>%
  bind_rows(prot_rnaLog2FPKM_spearman %>% filter(p_adjust<0.05) %>% rename(cor = p) %>% mutate(genes = "FDR<0.05", cor_type = "Spearman")) %>%
  mutate(batch = fct_reorder(batch, cor, .fun = function(x) {median(x, na.rm = T)}))

plot <- plot_data %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = batch, y = cor, fill = genes), outlier.shape = NA, lwd=0.2) +
  geom_text(
    data = plot_data %>%
      group_by(cor_type, batch, genes) %>%
      tally() %>%
      ungroup() %>%
      mutate(x = c(seq(0.5, 11, by = 0.5), seq(0.5, 11, by = 0.5))) %>%
      mutate(y = -0.6),
    mapping = aes(x = x, y = y, color = genes, label=n),
    show.legend = F, size=3, nudge_x = 0.2) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 12),
    legend.title = element_text(color = "black", size = 14),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 14)) +
  facet_wrap(~ cor_type) +
  labs(x = "Batch", y = "Correlation", fill = "Genes")

#+ fig.width=8, fig.height=5
plot

ggsave(filename="ProtRNAlog2fpkm_correlation_batch.png", plot = plot, path = "./output/plots/protein_rna_correlation/", width=8, height=5)
unlink("ProtRNAlog2fpkm_correlation_batch.png")




#' # load multi-omics data from MCP paper

#' corrected data (batch effects regressed-out)
mcp_corrected <- read_tsv(file = "./data/mcp_paper/protein_attenuation_cnv_rna_protein.txt")
metadata <- read_tsv(file = "./data/mcp_paper/metadata_cptac_cellLines.txt") %>%
  select(sample, batch) %>%
  mutate(batch  = str_replace(batch, "LAW", "ccle")) %>%
  mutate(batch  = str_replace(batch, "LPK", "ccle")) %>%
  mutate(batch  = str_replace(batch, "RMLT", "ccle")) 

mcp_corrected <- mcp_corrected %>%
  inner_join(metadata[, c("sample", "batch")], by = "sample") %>%
  group_by(batch, gene) %>%
  filter(sum(!is.na(prot_log2FC)) > 0.2*n()) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor_test = map(.x=data, .f = ~ broom::tidy(cor.test(.x$prot_log2FC, .x$rna_log2CPM, method = c("pearson"))))) %>%
  select(-data) %>%
  unnest(cols = cor_test) %>%
  select(batch, gene, r=estimate, p_value=p.value) %>%
  group_by(batch) %>%
  mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  ungroup()

mcp_corrected_plot <- mcp_corrected %>%
  mutate(batch = fct_reorder(batch, r, mean)) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = batch, y = r), outlier.shape = NA, lwd=0.2, fill = "grey") +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12)) +
  labs(x = "", y = "r", fill = "Genes")

#+ fig.width=4, fig.height=2
mcp_corrected_plot

ggsave(filename="mcp_corrected_plot.png", plot = mcp_corrected_plot, path = "./output/plots/protein_rna_correlation/", width=4, height=2)
unlink("mcp_corrected_plot.png")



#' uncorrected data
mcp_prot <- read_tsv(file = "./data/mcp_paper/proteomicsQ_cptac_cellLines.txt") %>%
  gather(key = "sample", value = "prot_log2FC", -gene)
mcp_rna <- read_tsv(file = "./data/mcp_paper/rna_tcga_cellLines.txt") %>%
  gather(key = "sample", value = "rna_log2CPM", -gene)

mcp_uncorrected <- mcp_prot %>%
  inner_join(mcp_rna, by = c("gene", "sample")) %>%
  inner_join(metadata[, c("sample", "batch")], by = "sample") %>%
  group_by(batch, gene) %>%
  filter(sum(!is.na(prot_log2FC)) > 0.2*n()) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor_test = map(.x=data, .f = ~ broom::tidy(cor.test(.x$prot_log2FC, .x$rna_log2CPM, method = c("pearson"))))) %>%
  select(-data) %>%
  unnest(cols = cor_test) %>%
  select(batch, gene, r=estimate, p_value=p.value) %>%
  group_by(batch) %>%
  mutate(p_adjust = p.adjust(p_value, method = "BH")) %>%
  ungroup()

mcp_uncorrected_plot <- mcp_uncorrected %>%
  mutate(batch = fct_reorder(batch, r, mean)) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = batch, y = r), outlier.shape = NA, lwd=0.2, fill = "grey") +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12)) +
  labs(x = "", y = "r", fill = "Genes")

#+ fig.width=4, fig.height=2
mcp_uncorrected_plot

ggsave(filename="mcp_uncorrected_plot.png", plot = mcp_uncorrected_plot, path = "./output/plots/protein_rna_correlation/", width=4, height=2)
unlink("mcp_uncorrected_plot.png")
