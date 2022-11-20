#' ---
#' title: "Phosphosite to protein/mRNA correlation across samples"
#' author: "abelsousa"
#' ---

#+ setup, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)



#' load R packages
library(RColorBrewer)
library(tidyverse)
library(cowplot)
library(viridis)

source("./src/utils/cor_functions.R")
source("./src/utils/getSamples.R")


#' load all samples gathered in this study
all_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")


#' get samples with phosphorylation and protein data
samplesPP <- getSamples(all_samples, c("phosphorylation", "protein", "clinical")) %>%
  select(sample, batch) %>%
  mutate(data = "Phospho/Protein")
nrow(samplesPP)

#' get samples with phosphorylation and mRNA data
samplesPR <- getSamples(all_samples, c("phosphorylation", "mRNA", "clinical")) %>%
  select(sample, batch) %>%
  mutate(data = "Phospho/mRNA")
nrow(samplesPR)


#' plot number of samples by batch
plot <- samplesPP %>%
  bind_rows(samplesPR) %>%
  mutate(data = fct_rev(fct_infreq(data))) %>%
  mutate(batch = fct_rev(fct_infreq(batch))) %>%
  ggplot() +
  geom_bar(mapping = aes(x = batch, y = ..count.., fill = data), position = "dodge") +
  theme_classic() +
  coord_flip() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 12),
    legend.title = element_blank()) +
  scale_fill_manual(values = c("#d8b365", "#5ab4ac")) +
  labs(x = "Study", y = "Number of samples")

#+ fig.width=6, fig.height=5
plot

ggsave(filename="PhoProtRNA_samples.png", plot = plot, path = "./output/plots/phospho_ProtRNA_correlation/", width=6, height=5)
unlink("PhoProtRNA_samples.png")



#' load protein data (log2FC)
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(names_to = "sample", values_to = "prot_log2fc", -gene) %>%
  select(sample, protein=gene, prot_log2fc)


#' load mRNA data (log2FC)
rna <- read_tsv(file = "./output/files/transcriptomics_log2fc.txt.gz") %>%
  pivot_longer(names_to = "sample", values_to = "rna_log2fc", -gene) %>%
  select(sample, protein=gene, rna_log2fc)


#' load quantile normalized phosphorylation data (log2FC)\
#' remove some phosphosites with NAs in all samples
phospho <- read_tsv(file = "./output/files/phosphoproteomicsQ.txt.gz") %>%
  pivot_longer(names_to = "sample", values_to = "phos_log2fc", -c(gene, psite, psites)) %>%
  select(sample, protein=gene, psite, psites, phos_log2fc) %>%
  group_by(protein, psite, psites) %>%
  filter( !(sum(is.na(phos_log2fc)) == n()) ) %>%
  ungroup()



#' bind protein and rna data together
prot_rna_data <- bind_rows(
  protein %>%
    rename(prot_fc=prot_log2fc) %>%
    mutate(prot_data = "protein_level"),
  rna %>%
    rename(prot_fc=rna_log2fc) %>%
    mutate(prot_data = "rna_level"))


#' # number of proteins by protein and rna data (filter at the protein level)
proteins <- bind_rows(
  prot_rna_data %>%
    group_by(prot_data) %>%
    summarise(n = length(unique(protein))) %>%
    mutate(filter = "None"),
  prot_rna_data %>%
    group_by(prot_data, protein) %>%
    filter(sum(!is.na(prot_fc)) > 0.1*n()) %>%
    ungroup() %>%
    group_by(prot_data) %>%
    summarise(n = length(unique(protein))) %>%
    mutate(filter = "expression in 10% samples"),
  prot_rna_data %>%
    group_by(prot_data, protein) %>%
    filter(sum(!is.na(prot_fc)) > 0.2*n()) %>%
    ungroup() %>%
    group_by(prot_data) %>%
    summarise(n = length(unique(protein))) %>%
    mutate(filter = "expression in 20% samples"),
  prot_rna_data %>%
    group_by(prot_data, protein) %>%
    filter(sum(!is.na(prot_fc)) > 0.3*n()) %>%
    ungroup() %>%
    group_by(prot_data) %>%
    summarise(n = length(unique(protein))) %>%
    mutate(filter = "expression in 30% samples"),
  prot_rna_data %>%
    group_by(prot_data, protein) %>%
    filter(sum(!is.na(prot_fc)) > 0.4*n()) %>%
    ungroup() %>%
    group_by(prot_data) %>%
    summarise(n = length(unique(protein))) %>%
    mutate(filter = "expression in 40% samples"),
  prot_rna_data %>%
    group_by(prot_data, protein) %>%
    filter(sum(!is.na(prot_fc)) > 0.5*n()) %>%
    ungroup() %>%
    group_by(prot_data) %>%
    summarise(n = length(unique(protein))) %>%
    mutate(filter = "expression in 50% samples"))


#' plot number of proteins by data
plot <- proteins %>%
  mutate(prot_data = fct_reorder(.f=prot_data, .x=n, .fun=function(x){median(x)})) %>%
  mutate(filter = fct_reorder(.f=filter, .x=n, .fun=function(x){median(x)})) %>%
  ggplot() +
  geom_bar(mapping = aes(x = prot_data, y = n, fill = filter), position = "dodge", stat = "identity") +
  theme_classic() +
  coord_flip() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 12)) +
  scale_fill_viridis(discrete = T) +
  scale_x_discrete(labels = c("protein", "mRNA")) +
  labs(x = "Protein level data", y = "Number of proteins")

#+ fig.width=8, fig.height=4
plot

ggsave(filename="protein_RNA_data_proteins.png", plot = plot, path = "./output/plots/phospho_ProtRNA_correlation/", width=8, height=4)
unlink("protein_RNA_data_proteins.png")

rm(prot_rna_data)


#' # number of psites/proteins in the phosphorylation data (filter at the psite level)
psites <- bind_rows(
  phospho %>%
    summarise(psite = length(unique(psite)), protein = length(unique(protein))) %>%
    mutate(filter = "none"),
  phospho %>%
    group_by(psite) %>%
    filter(sum(!is.na(phos_log2fc)) > 0.1*n()) %>%
    ungroup() %>%
    summarise(psite = length(unique(psite)), protein = length(unique(protein))) %>%
    mutate(filter = "10% samples"),
  phospho %>%
    group_by(psite) %>%
    filter(sum(!is.na(phos_log2fc)) > 0.2*n()) %>%
    ungroup() %>%
    summarise(psite = length(unique(psite)), protein = length(unique(protein))) %>%
    mutate(filter = "20% samples"),
  phospho %>%
    group_by(psite) %>%
    filter(sum(!is.na(phos_log2fc)) > 0.3*n()) %>%
    ungroup() %>%
    summarise(psite = length(unique(psite)), protein = length(unique(protein))) %>%
    mutate(filter = "30% samples"),
  phospho %>%
    group_by(psite) %>%
    filter(sum(!is.na(phos_log2fc)) > 0.4*n()) %>%
    ungroup() %>%
    summarise(psite = length(unique(psite)), protein = length(unique(protein))) %>%
    mutate(filter = "40% samples"),
  phospho %>%
    group_by(psite) %>%
    filter(sum(!is.na(phos_log2fc)) > 0.5*n()) %>%
    ungroup() %>%
    summarise(psite = length(unique(psite)), protein = length(unique(protein))) %>%
    mutate(filter = "50% samples")) %>%
  select(filter, everything()) %>%
  pivot_longer(-filter, names_to="data", values_to="n")



#' plot number of proteins and phosphosites
plot <- psites %>%
  mutate(data = fct_reorder(.f=data, .x=n, .fun=function(x){median(x)})) %>%
  mutate(filter = fct_reorder(.f=filter, .x=n, .fun=function(x){median(x)})) %>%
  ggplot() +
  geom_bar(mapping = aes(x = filter, y = n, fill = data), position = "dodge", stat = "identity") +
  facet_wrap(~ data, scales = "free_x") +
  theme_classic() +
  coord_flip() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 12),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 14)) +
  scale_fill_viridis(discrete = T, guide=F) +
  labs(x = "Filter", y = "Number of features")

#+ fig.width=8, fig.height=4
plot

ggsave(filename="phospho_number_psites_prots.png", plot = plot, path = "./output/plots/phospho_ProtRNA_correlation/", width=8, height=4)
unlink("phospho_number_psites_prots.png")



#' add phosphorylation and protein data to respective samples
pho_prot <- samplesPP %>%
  select(-data) %>%
  inner_join(protein, by = "sample") %>%
  inner_join(phospho, by = c("sample", "protein")) %>%
  select(batch, sample, prot=protein, psite, psites, prot_fc=prot_log2fc, psite_fc=phos_log2fc) %>%
  mutate(prot_data = "protein_level")


#' add phosphorylation and mRNA data to respective samples
pho_rna <- samplesPR %>%
  select(-data) %>%
  inner_join(rna, by = "sample") %>%
  inner_join(phospho, by = c("sample", "protein")) %>%
  select(batch, sample, prot=protein, psite, psites, prot_fc=rna_log2fc, psite_fc=phos_log2fc) %>%
  mutate(prot_data = "rna_level")


#' bind all data together
phospho_data <- bind_rows(pho_prot, pho_rna)

rm(phospho, protein, rna, pho_prot, pho_rna)



#' # number of proteins by phospho/protein and phospho/rna data\
#' # (filter at the phosphosite level)
proteins <- bind_rows(
  phospho_data %>%
    group_by(prot_data) %>%
    summarise(n = length(unique(prot))) %>%
    mutate(filter = "None"),
  phospho_data %>%
    group_by(prot_data, psite) %>%
    filter(sum(!is.na(prot_fc) & !is.na(psite_fc)) > 0.1*n()) %>%
    ungroup() %>%
    group_by(prot_data) %>%
    summarise(n = length(unique(prot))) %>%
    mutate(filter = "expression in 10% samples"),
  phospho_data %>%
    group_by(prot_data, psite) %>%
    filter(sum(!is.na(prot_fc) & !is.na(psite_fc)) > 0.2*n()) %>%
    ungroup() %>%
    group_by(prot_data) %>%
    summarise(n = length(unique(prot))) %>%
    mutate(filter = "expression in 20% samples"),
  phospho_data %>%
    group_by(prot_data, psite) %>%
    filter(sum(!is.na(prot_fc) & !is.na(psite_fc)) > 0.3*n()) %>%
    ungroup() %>%
    group_by(prot_data) %>%
    summarise(n = length(unique(prot))) %>%
    mutate(filter = "expression in 30% samples"),
  phospho_data %>%
    group_by(prot_data, psite) %>%
    filter(sum(!is.na(prot_fc) & !is.na(psite_fc)) > 0.4*n()) %>%
    ungroup() %>%
    group_by(prot_data) %>%
    summarise(n = length(unique(prot))) %>%
    mutate(filter = "expression in 40% samples"),
  phospho_data %>%
    group_by(prot_data, psite) %>%
    filter(sum(!is.na(prot_fc) & !is.na(psite_fc)) > 0.5*n()) %>%
    ungroup() %>%
    group_by(prot_data) %>%
    summarise(n = length(unique(prot))) %>%
    mutate(filter = "expression in 50% samples"))


#' plot number of proteins by data
plot <- proteins %>%
  mutate(prot_data = fct_reorder(.f=prot_data, .x=n, .fun=function(x){median(x)})) %>%
  mutate(filter = fct_reorder(.f=filter, .x=n, .fun=function(x){median(x)})) %>%
  ggplot() +
  geom_bar(mapping = aes(x = prot_data, y = n, fill = filter), position = "dodge", stat = "identity") +
  theme_classic() +
  coord_flip() +
  theme(
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 12)) +
  scale_fill_viridis(discrete = T) +
  scale_x_discrete(labels = c("protein", "mRNA")) +
  labs(x = "Protein level data", y = "Number of proteins")

#+ fig.width=8, fig.height=4
plot

ggsave(filename="PhoProtRNA_proteins.png", plot = plot, path = "./output/plots/phospho_ProtRNA_correlation/", width=8, height=4)
unlink("PhoProtRNA_proteins.png")






#' # phosphosite to protein/rna correlation


#' ## correlation across samples by study
cor <- phospho_data %>%
  group_by(prot_data, batch, psite, psites) %>%
  filter(sum(!(is.na(prot_fc) | is.na(psite_fc))) > 0.1*n()) %>%
  nest() %>%
  ungroup() %>%
  mutate(r = map(.x=data, .f = ~ cor_test(.x$prot_fc, .x$psite_fc, method = "pearson"))) %>%
  select(-data) %>%
  unnest(cols = r) %>%
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
  geom_hline(yintercept = 0.5, linetype="dashed", color = "black", size=0.4) +
  geom_text(
    data = plot_data %>%
      group_by(prot_data, batch, genes) %>%
      tally() %>%
      ungroup() %>%
      mutate(x = c(seq(0.5, 10, by = 0.5), seq(0.5, 10, by = 0.5))) %>%
      mutate(y = -0.6),
    mapping = aes(x = x, y = y, color = genes, label=n),
    show.legend = F, size=3, nudge_x = 0.2) +
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

ggsave(filename="PhosphoProtRNAlog2fc_correlation_batch.png", plot = plot, path = "./output/plots/phospho_ProtRNA_correlation/", width=8, height=5)
unlink("PhosphoProtRNAlog2fc_correlation_batch.png")



#' ## correlation across samples (all studies together)
cor <- phospho_data %>%
  group_by(prot_data, psite, psites) %>%
  filter(sum(!(is.na(prot_fc) | is.na(psite_fc))) > 0.1*n()) %>%
  nest() %>%
  ungroup() %>%
  mutate(r = map(.x=data, .f = ~ cor_test(.x$prot_fc, .x$psite_fc, method = "pearson"))) %>%
  select(-data) %>%
  unnest(cols = r) %>%
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
  geom_hline(yintercept = 0.5, linetype="dashed", color = "black", size=0.4) +
  geom_text(
    data = plot_data %>%
      group_by(prot_data, genes) %>%
      tally() %>%
      ungroup() %>%
      mutate(x = c(seq(0.5, 1, by = 0.5), seq(0.5, 1, by = 0.5))) %>%
      mutate(y = -0.4),
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

ggsave(filename="PhosphoProtRNAlog2fc_correlation.png", plot = plot, path = "./output/plots/phospho_ProtRNA_correlation/", width=5, height=4)
unlink("PhosphoProtRNAlog2fc_correlation.png")


