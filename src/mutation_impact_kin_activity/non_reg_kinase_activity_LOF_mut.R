#' ---
#' title: "Kinase activity in samples with LOF mutations - non-reg-out phospho data"
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
mut <- read_tsv(file = "./output/files/mutations_protpos.txt.gz")

mut1 <- mut %>%
  select(sample, variant_class, gene = gene_symbol) %>%
  distinct()

mut2 <- mut %>%
  select(sample, variant_class, gene = gene_symbol, prot_pos, seq) %>%
  distinct() %>%
  group_by(sample, gene, variant_class) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)


#' load kinase-activity inference data\
#' (quantile-normalized non-regressed-out phosphorylation data)
ka <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ.txt.gz")

k_subN <- 3
k_sampN <- 0
ka_dbt <- ka %>%
  filter(source_type == "DB_text-mining") %>%
  #filter(source_type == "in_vivo") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  group_by(kinase) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > k_sampN) %>%
  select(-n) %>%
  #pivot_wider(names_from = "sample", values_from = "log10P") %>%
  #pivot_longer(cols = -c(kinase), names_to = "sample", values_to = "log10P") %>%
  select(kinase, sample, log10P)



ka_mut <- ka_dbt %>%
  inner_join(mut1, by = c("kinase" = "gene", "sample")) %>%
  mutate(variant_class = fct_reorder(.f=variant_class, .x=log10P, .fun=function(x) median(x, na.rm = T))) %>%
  ggplot(mapping = aes(x = variant_class, y = log10P)) +
  geom_boxplot(notch = T, outlier.size = 0.1, lwd=0.3) +
  theme_classic() +
  coord_flip() +
  stat_compare_means(method = "wilcox.test", comparisons = list(c(3, 7))) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "Variant class", y = "log10P")

#+ fig.width=8, fig.height=2
ka_mut

ggsave(filename = "non_reg_kinase_activity_lof_mut.png", plot = ka_mut, path = "./output/plots/mutation_impact_kin_activity/", height = 2, width = 8)
unlink("non_reg_kinase_activity_lof_mut.png")


ka_mut <- ka_dbt %>%
  inner_join(mut1, by = c("kinase" = "gene", "sample")) %>%
  mutate(is_missense = if_else(variant_class == "Missense_Mutation", "yes", "no")) 

ka_mut_no <- ka_mut %>%
  filter(is_missense == "no") %>%
  filter(!is.na(log10P))

ka_mut_yes <- ka_mut %>%
  filter(is_missense == "yes") %>%
  filter(!is.na(log10P)) %>%
  sample_n(nrow(ka_mut_no))

ka_mut <- ka_mut_no %>%
  bind_rows(ka_mut_yes) %>%
  ggplot(mapping = aes(x = log10P, fill = is_missense, color = is_missense)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "log10P", y = "Density")

#+ fig.width=8, fig.height=2
ka_mut

ggsave(filename = "non_reg_kinase_activity_isMissense.png", plot = ka_mut, path = "./output/plots/mutation_impact_kin_activity/", height = 2, width = 8)
unlink("non_reg_kinase_activity_isMissense.png")




ka_mutLOF <- ka_dbt %>%
  inner_join(mut1, by = c("kinase" = "gene", "sample")) %>%
  filter(str_detect(variant_class, "Frame|Nonsense")) %>%
  mutate(class = "lof_mutation") %>%
  filter(!is.na(log10P))

ka_bk <- ka_dbt %>%
  anti_join(ka_mutLOF, by = c("kinase", "sample", "log10P")) %>%
  mutate(class = "background") %>%
  filter(!is.na(log10P)) %>%
  sample_n(nrow(ka_mutLOF))

ka_mut <- bind_rows(
  ka_mutLOF, ka_bk) %>%
  ggplot(mapping = aes(x = log10P, fill = class, color = class)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  #stat_compare_means(method = "wilcox.test", comparisons = list(c(1, 7), c(2, 7), c(5, 7))) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "log10P", y = "Density")

#+ fig.width=8, fig.height=2
ka_mut

ggsave(filename = "non_reg_kinase_activity_lofmut_background.png", plot = ka_mut, path = "./output/plots/mutation_impact_kin_activity/", height = 2, width = 8)
unlink("non_reg_kinase_activity_lofmut_background.png")



ka_mut <- ka_dbt %>%
  inner_join(mut2, by = c("kinase" = "gene", "sample")) %>%
  filter(variant_class == "Nonsense_Mutation", !is.na(log10P)) %>%
  mutate(prot_length = nchar(seq)) %>%
  #filter(prot_length < 1000) %>%
  mutate(prot_pos_p = prot_pos/prot_length) %>%
  select(-seq, -prot_length) %>%
  pivot_longer(c(prot_pos, prot_pos_p), names_to = "position", values_to = "value") %>%
  ggplot(mapping = aes(x = value, y = log10P)) +
  facet_wrap( ~ position, scales = "free_x") +
  geom_point(size = 0.5) +
  theme_classic() +
  stat_cor(size = 3) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10),
    aspect.ratio = 1) +
  labs(x = "Position", y = "log10P")

#+ fig.width=4, fig.height=4
ka_mut

ggsave(filename = "non_reg_kinase_activity_nonsense_pos_cor.png", plot = ka_mut, path = "./output/plots/mutation_impact_kin_activity/", height = 4, width = 4)
unlink("non_reg_kinase_activity_nonsense_pos_cor.png")


