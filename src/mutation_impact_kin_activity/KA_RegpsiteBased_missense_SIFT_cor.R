#' ---
#' title: "Missense mutation SIFT score to kinase activity correlation (regulatory phosphosite-based kinase activity)"
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
mut <- read_tsv(file = "./output/files/mutations_protpos.txt.gz") %>%
  filter(variant_class == "Missense_Mutation" & variant_type == "SNP") %>%
  select(sample, gene_symbol, prot_pos, aa_wt, aa_mut) %>%
  distinct() %>%
  rename(gene = gene_symbol, pos = prot_pos, ref = aa_wt, alt = aa_mut) %>%
  group_by(sample, gene) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)


#' load kinase-activity inference data\
#' (quantile-normalized protein regressed-out phosphorylation data)
ka <- read_tsv(file = "./output/files/CPTAC_KA_kin_psites_phosphoQ_Protregout_allSamp_withZtransf.txt.gz")

k <- 3
ka <- ka %>%
  rename(kinase = gene) %>%
  filter(n >= k) %>%
  select(-n)

#' ## correlate SIFT score (-log10) to kinase activity
sift2 <- sift %>%
  filter(gene %in% ka$kinase)

ka_sift <- ka %>%
  inner_join(mut, by = c("sample", "kinase" = "gene")) %>%
  inner_join(sift2, by = c("kinase" = "gene", "pos", "ref", "alt")) %>%
  ggplot(mapping = aes(x = score, y = log10P)) +
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
  labs(x = "SIFT score (-log10)", y = "Kinase activity (log10)")


#+ fig.width=4, fig.height=4
ka_sift

ggsave(filename = "kinase_activity_psite_based_missense_sift_cor.png", plot = ka_sift, path = "./output/plots/mutation_impact_kin_activity/", height = 4, width = 4)
unlink("kinase_activity_psite_based_missense_sift_cor.png")
