#' ---
#' title: "Protein abundance distribution - nonsense vs missense mutations"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=6, fig.height=6)


#' load R packages
library(tidyverse)
library(ggpubr)

set.seed(123)

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load mutation data
#mut <- read_tsv(file = "./output/files/mutations_protpos.txt.gz") %>%
#  filter(variant_class == "Nonsense_Mutation") %>%
#  select(sample, gene = gene_symbol, variant_type, dna_pos=POS, ref=REF, alt = ALT, prot_pos, aa_wt, aa_mut)
nonsense_mut <- read_tsv(file = "./output/files/mutations_protpos.txt.gz") %>%
  filter(variant_class == "Nonsense_Mutation") %>%
  select(sample, gene = gene_symbol) %>%
  distinct()

missense_mut <- read_tsv(file = "./output/files/mutations_protpos.txt.gz") %>%
  filter(variant_class == "Missense_Mutation") %>%
  select(sample, gene = gene_symbol) %>%
  distinct()

exclude <- intersect(missense_mut, nonsense_mut)

nonsense_mut <- nonsense_mut %>%
  anti_join(exclude, by = c("sample", "gene"))

missense_mut <- missense_mut %>%
  anti_join(exclude, by = c("sample", "gene"))


#' load protein abundance
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "log2fc")


nonsense_samples <- protein %>%
  semi_join(nonsense_mut, by = c("sample", "gene")) %>%
  mutate(class = "nonsense_mutation") %>%
  filter(!is.na(log2fc))

missense_samples <- protein %>%
  semi_join(missense_mut, by = c("sample", "gene")) %>%
  mutate(class = "missense_mutation") %>%
  filter(!is.na(log2fc))
  #sample_n(nrow(nonsense_samples))


nonsense_missense <- nonsense_samples %>%
  bind_rows(missense_samples) %>%
  ggplot(mapping = aes(x = log2fc, fill = class, color = class)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "Protein abundance (log2fc)", y = "Density")

#+ fig.width=8, fig.height=2
nonsense_missense

ggsave(filename = "prot_abundance_nonsense_missense_density1.png", plot = nonsense_missense, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)
unlink("prot_abundance_nonsense_missense_density1.png")


nonsense_missense <- nonsense_samples %>%
  bind_rows(missense_samples) %>%
  mutate(class = fct_recode(class, Nonsense = "nonsense_mutation", Missense = "missense_mutation")) %>%
  ggplot(mapping = aes(x = class, y = log2fc, fill = class)) +
  geom_boxplot(notch = T, alpha = 0.8, color = "grey60") +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()) +
  coord_flip() +
  scale_fill_discrete(guide=F) +
  stat_compare_means(method = "t.test", label.x = 2.4, label.y = 10, size = 3) +
  labs(x = "", y = "Protein abundance (log2fc)")

#+ fig.width=8, fig.height=2
nonsense_missense

ggsave(filename = "prot_abundance_nonsense_missense_box1.png", plot = nonsense_missense, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)
unlink("prot_abundance_nonsense_missense_box1.png")


nonsense_missense <- nonsense_samples %>%
  bind_rows(missense_samples) %>%
  ggplot(mapping = aes(x = log2fc, fill = class, color = class)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "Protein abundance (log2fc)", y = "Density") +
  scale_x_continuous(limits = c(-2,2))

#+ fig.width=8, fig.height=2
nonsense_missense

ggsave(filename = "prot_abundance_nonsense_missense_density2.png", plot = nonsense_missense, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)
unlink("prot_abundance_nonsense_missense_density2.png")


nonsense_missense <- nonsense_samples %>%
  bind_rows(missense_samples) %>%
  mutate(class = fct_recode(class, Nonsense = "nonsense_mutation", Missense = "missense_mutation")) %>%
  ggplot(mapping = aes(x = class, y = log2fc, fill = class)) +
  geom_boxplot(notch = T, alpha = 0.8, color = "grey60") +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()) +
  coord_flip() +
  scale_fill_discrete(guide=F) +
  scale_y_continuous(limits = c(-2,2)) +
  stat_compare_means(method = "t.test", label.x = 2.4, label.y = 1, size = 3) +
  labs(x = "", y = "Protein abundance (log2fc)")

#+ fig.width=8, fig.height=2
nonsense_missense

ggsave(filename = "prot_abundance_nonsense_missense_box2.png", plot = nonsense_missense, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)
unlink("prot_abundance_nonsense_missense_box2.png")
