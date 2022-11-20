#' ---
#' title: "Missense mutation SIFT score to protein abundance correlation"
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


#' load SIFT data
sift <- data.table::fread(file = "./data/mutfunc/human/sift.tab.gz") %>%
  as_tibble()

# load uniprot to gene name mapping
sift_map <- read_tsv(file = "./data/mutfunc/human/sift_uniprot_acc_genename.tab") %>%
  rename(acc = From, gene = To)

sift <- sift %>%
  inner_join(sift_map, by = c("acc")) %>%
  select(acc, gene, everything())

exclude <- sift %>%
  select(gene, acc) %>%
  distinct() %>%
  group_by(gene) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  filter(n > 1)

sift <- sift %>%
  #anti_join(exclude, by = "gene")
  filter(!(gene %in% exclude$gene)) %>%
  mutate(score = -log10(score+1e-10)) %>%
  select(gene, pos, ref, alt, score)


#' load mutation data
mut <- read_tsv(file = "./output/files/mutations_protpos.txt.gz")

mut <- mut %>%
  filter(variant_class == "Missense_Mutation" & variant_type == "SNP") %>%
  select(-c(variant_type, variant_class, CHROM, POS, REF, ALT, HGVSp, seq)) %>%
  distinct() %>%
  rename(gene = gene_symbol, pos = prot_pos, ref = aa_wt, alt = aa_mut) %>%
  group_by(sample, gene) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)


#' load protein abundance
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "log2fc")


#' load rna abundance
rna <- read_tsv(file = "./output/files/transcriptomics_log2fc.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "log2fc")


#' ## correlate SIFT score (-log10) to protein abundance (log2fc)
sift2 <- sift %>%
  filter(gene %in% protein$gene)

prot_sift <- protein %>%
  inner_join(mut, by = c("sample", "gene")) %>%
  inner_join(sift2, by = c("gene", "pos", "ref", "alt")) %>%
  filter(!is.na(log2fc)) %>%
  ggplot(mapping = aes(x = score, y = log2fc)) +
  geom_point(size = 0.5) +
  theme_classic() +
  stat_cor(size = 3) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10),
    aspect.ratio = 1) +
  labs(x = "SIFT score (-log10)", y = "Protein abundance (log2fc)")


#+ fig.width=4, fig.height=4
prot_sift

ggsave(filename = "protein_abundance_missense_sift_cor.png", plot = prot_sift, path = "./output/plots/mutation_impact_protein_abundance/", height = 4, width = 4)
unlink("protein_abundance_missense_sift_cor.png")



#' ## correlate SIFT score (-log10) to protein/rna abundance ratio (log2fc)
protein <- protein %>%
  rename(p_log2fc = log2fc) %>%
  inner_join(rna, by = c("gene", "sample")) %>%
  rename(r_log2fc = log2fc) %>%
  mutate(pr_log2fc = p_log2fc - r_log2fc)

sift2 <- sift %>%
  filter(gene %in% protein$gene)

prot_sift <- protein %>%
  inner_join(mut, by = c("sample", "gene")) %>%
  inner_join(sift2, by = c("gene", "pos", "ref", "alt")) %>%
  filter(!is.na(pr_log2fc)) %>%
  ggplot(mapping = aes(x = score, y = pr_log2fc)) +
  geom_point(size = 0.5) +
  theme_classic() +
  stat_cor(size = 3) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10),
    aspect.ratio = 1) +
  labs(x = "SIFT score (-log10)", y = "Protein/rna abundance ratio (log2fc)")


#+ fig.width=4, fig.height=4
prot_sift

ggsave(filename = "protrna_abundance_ratio_missense_sift_cor.png", plot = prot_sift, path = "./output/plots/mutation_impact_protein_abundance/", height = 4, width = 4)
unlink("protrna_abundance_ratio_missense_sift_cor.png")
