#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)

#' load R packages
library(tidyverse)

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load CPTAC samples
source("./src/utils/getSamples.R")
cptac_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")
cptac_samples <- getSamples(cptac_samples, c("protein")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = str_c(batch, tissue, sep = "-")) %>%
  select(-tissue)


#' load protein abundance
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "log2fc") %>%
  filter(!is.na(log2fc))


#' load CPTAC mutations
cptac_mutations <- read_tsv(file = "./output/files/mutations_protpos.txt.gz")


#' select only nonsense mutations\
#' select only the gene, sample, protein position and protein sequence\
#' remove duplicates (same gene in the same sample can have multiple nonsense mutations that map to the same protein position)
nonsense <- cptac_mutations %>%
  filter(variant_class == "Nonsense_Mutation") %>%
  select(sample, gene=gene_symbol, prot_pos_i=prot_pos, seq) %>%
  distinct()


#' select only frameshift mutations\
#' select only the gene, sample, protein position and protein sequence\
#' remove duplicates (same gene in the same sample can have multiple frameshift mutations that map to the same protein position)
frameshift <- cptac_mutations %>%
  filter(variant_class == "Frame_Shift_Del" | variant_class == "Frame_Shift_Ins") %>%
  select(sample, gene=gene_symbol, prot_pos_i=prot_pos, seq) %>%
  distinct()


#' # -- assess the impact of Nonsense mutations on protein abundance
nonsense_impact <- nonsense %>%
  inner_join(protein, by = c("sample", "gene")) %>%
  group_by(sample, gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  mutate(protei_len = nchar(seq), pos_p = (prot_pos_i/protei_len)*100) %>%
  pivot_longer(cols = c("pos_p", "prot_pos_i"), names_to = "position_type", values_to = "position") %>%
  ggplot(mapping = aes(x = position, y = log2fc)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggpubr::stat_cor(method = "pearson") +
  facet_wrap(~ position_type, ncol = 2, scales = "free", labeller = labeller(position_type = c("pos_p"="%", "prot_pos_i"="aa")))

#+ fig.width=8, fig.height=4
nonsense_impact

ggsave(filename = "protein_abundance_nonsense_position_cor.png", plot = nonsense_impact, path = "./output/plots/mutation_impact_protein_abundance/", height = 4, width = 8)
unlink("protein_abundance_nonsense_position_cor.png")


#' # -- assess the impact of Frameshift mutations on TF activities
frameshift_impact <- frameshift %>%
  inner_join(protein, by = c("sample", "gene")) %>%
  group_by(sample, gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  mutate(protei_len = nchar(seq), pos_p = (prot_pos_i/protei_len)*100) %>%
  pivot_longer(cols = c("pos_p", "prot_pos_i"), names_to = "position_type", values_to = "position") %>%
  ggplot(mapping = aes(x = position, y = log2fc)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggpubr::stat_cor(method = "pearson") +
  facet_wrap(~ position_type, ncol = 2, scales = "free", labeller = labeller(position_type = c("pos_p"="%", "prot_pos_i"="aa")))

#+ fig.width=8, fig.height=4
frameshift_impact

ggsave(filename = "protein_abundance_frameshift_position_cor.png", plot = frameshift_impact, path = "./output/plots/mutation_impact_protein_abundance/", height = 4, width = 8)
unlink("protein_abundance_frameshift_position_cor.png")
