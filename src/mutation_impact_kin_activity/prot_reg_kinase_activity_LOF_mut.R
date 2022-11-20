#' ---
#' title: "Kinase activity in samples with LOF mutations - protein reg-out phospho data"
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


#' load tumour purity scores
#purity <- read_tsv(file = "./output/files/samples_purity.txt")
#sel_samples <- purity %>%
#  filter(score > 0.9)


#' load kinase-activity inference data\
#' (quantile-normalized protein regressed-out phosphorylation data)
ka <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz")

k_subN <- 3
ka_dbt <- ka %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  select(kinase, sample, log10P)


#' compare kinase activation distribution between LOF mutations
mut1 <- mut %>%
  select(sample, variant_class, gene = gene_symbol) %>%
  distinct() %>%
  #filter(sample %in% sel_samples$sample)
  # remove cases of samples with multiple mutations on the same gene
  group_by(sample, gene) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)

ka_mut <- ka_dbt %>%
  inner_join(mut1, by = c("kinase" = "gene", "sample")) %>%
  mutate(variant_class = fct_reorder(.f=variant_class, .x=log10P, .fun=function(x) median(x, na.rm = T))) %>%
  ggplot(mapping = aes(x = variant_class, y = log10P)) +
  geom_boxplot(mapping = aes(fill = variant_class), notch = F, outlier.size = 0.1, lwd=0.3, outlier.shape = NA, show.legend = F) +
  geom_jitter(width = 0.1, alpha = 0.2, size = 0.1) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.title = element_text(colour = "black", size = 8),
    axis.text = element_text(colour = "black", size = 6)) +
  scale_y_continuous(limits = c(-6,6)) +
  labs(x = "Variant class", y = "Kinase activity (log10P)")

#+ fig.width=5, fig.height=2
ka_mut

ggsave(filename = "prot_reg_kinase_activity_lof_mut.png", plot = ka_mut, path = "./output/plots/mutation_impact_kin_activity/", height = 2, width = 5)
unlink("prot_reg_kinase_activity_lof_mut.png")


ka_mut <- ka_dbt %>%
  inner_join(mut1, by = c("kinase" = "gene", "sample")) %>%
  mutate(is_missense = if_else(variant_class == "Missense_Mutation", "yes", "no")) %>%
  ggplot(mapping = aes(x = log10P, fill = is_missense, color = is_missense)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  scale_x_continuous(limits = c(-6,6)) +
  scale_y_continuous(limits = c(0,0.8)) +
  labs(x = "Kinase activity (log10P)", y = "Density", fill = "Missense mutation", color = "Missense mutation")

#+ fig.width=8, fig.height=2
ka_mut

ggsave(filename = "prot_reg_kinase_activity_isMissense.png", plot = ka_mut, path = "./output/plots/mutation_impact_kin_activity/", height = 2, width = 8)
unlink("prot_reg_kinase_activity_isMissense.png")


ka_lof <- ka_dbt %>%
  inner_join(mut1, by = c("kinase" = "gene", "sample")) %>%
  filter(str_detect(variant_class, "Frame|Nonsense")) %>%
  mutate(class = "lof_mutation")

ka_bk <- ka_dbt %>%
  anti_join(ka_lof, by = c("kinase", "sample", "log10P")) %>%
  mutate(class = "background")

ka_mut <- bind_rows(ka_lof, ka_bk) %>%
  ggplot(mapping = aes(x = log10P, fill = class, color = class)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  scale_fill_discrete(label = c("Background", "LoF mutations:\nIn Frame\nFrame shit\nNonsense")) +
  scale_color_discrete(label = c("Background", "LoF mutations:\nIn Frame\nFrame shit\nNonsense")) +
  labs(x = "Kinase activity (log10P)", y = "Density")

#+ fig.width=8, fig.height=2
ka_mut

ggsave(filename = "prot_reg_kinase_activity_lofmut_background.png", plot = ka_mut, path = "./output/plots/mutation_impact_kin_activity/", height = 2, width = 8)
unlink("prot_reg_kinase_activity_lofmut_background.png")


ka_lof <- ka_dbt %>%
  inner_join(mut1, by = c("kinase" = "gene", "sample")) %>%
  mutate(class = "lof_mutation")

ka_bk <- ka_dbt %>%
  anti_join(ka_lof, by = c("kinase", "sample", "log10P")) %>%
  mutate(class = "background")

ka_mut <- bind_rows(ka_lof, ka_bk) %>%
  ggplot(mapping = aes(x = log10P, fill = class, color = class)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  scale_fill_discrete(label = c("Background", "LoF mutations")) +
  scale_color_discrete(label = c("Background", "LoF mutations")) +
  labs(x = "Kinase activity (log10P)", y = "Density")

#+ fig.width=8, fig.height=2
ka_mut

ggsave(filename = "prot_reg_kinase_activity_lofmut_background2.png", plot = ka_mut, path = "./output/plots/mutation_impact_kin_activity/", height = 2, width = 8)
unlink("prot_reg_kinase_activity_lofmut_background2.png")


#' correlate kinase activity with the position of nonsense mutations
mut2 <- mut %>%
  filter(variant_class == "Nonsense_Mutation") %>%
  select(sample, variant_class, gene = gene_symbol, prot_pos, seq) %>%
  distinct() %>%
  # remove cases of samples with multiple mutations on the same gene
  group_by(sample, gene) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)
  #filter(sample %in% sel_samples$sample)


ka_mut <- ka_dbt %>%
  inner_join(mut2, by = c("kinase" = "gene", "sample")) %>%
  mutate(prot_length = nchar(seq)) %>%
  mutate(prot_pos_p = prot_pos/prot_length) %>%
  select(-seq, -prot_length) %>%
  pivot_longer(c(prot_pos, prot_pos_p), names_to = "position", values_to = "value") %>%
  ggplot(mapping = aes(x = value, y = log10P)) +
  facet_wrap( ~ position, scales = "free_x", labeller = labeller(position = c("prot_pos" = "Protein position (aa)", "prot_pos_p" = "Protein position (%)"))) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  theme_classic() +
  stat_cor(size = 3) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10),
    panel.spacing = unit(1, "lines")) +
  labs(x = "Nonsense mutation position", y = "Kinase activity (log10P)")

#+ fig.width=4, fig.height=4
ka_mut

ggsave(filename = "prot_reg_kinase_activity_nonsense_pos_cor.png", plot = ka_mut, path = "./output/plots/mutation_impact_kin_activity/", height = 4, width = 4)
unlink("prot_reg_kinase_activity_nonsense_pos_cor.png")


#' correlate kinase activity with the position of nonsense mutations
#' by gene
mut2 <- mut %>%
  filter(variant_class == "Nonsense_Mutation") %>%
  select(sample, gene = gene_symbol, prot_pos) %>%
  inner_join(ka_dbt, by = c("gene" = "kinase", "sample")) %>%
  group_by(gene) %>%
  mutate(n = n()) %>%
  filter(n == 7) %>%
  ungroup()



#' # kinase activity distribution between LoF mutations inside/outside kinase domains

#' load kinase domains
kinDomains <- read_tsv(file = "./output/files/kinase_domain_positions.txt") %>%
  group_by(kinase) %>%
  nest(.key = "domains") %>%
  ungroup()

#' ## across all mutations
mutKin <- mut %>%
  select(sample, variant_class, gene = gene_symbol, pos = prot_pos) %>%
  distinct() %>%
  filter(!is.na(pos)) %>%
  # remove cases of samples with multiple mutations on the same gene
  #group_by(sample, gene) %>%
  #mutate(n = n()) %>%
  #ungroup() %>%
  #filter(n == 1) %>%
  #select(-n) %>%
  inner_join(kinDomains, by = c("gene" = "kinase")) %>%
  mutate(inside = map2_lgl(.x = pos, .y = domains, .f = ~ .y %>% mutate(pos = .x, inside = (pos >= start & pos <= end)) %>% pull(inside) %>% any())) %>%
  group_by(sample, gene) %>%
  summarise(inside = any(inside)) %>%
  ungroup() %>%
  mutate_if(is_logical, as.character) %>%
  rename(kinase = gene)


#' inside vs outside lof mutation
ka_mut <- ka_dbt %>%
  inner_join(mutKin, by = c("kinase", "sample")) %>%
  ggplot(mapping = aes(x = inside, y = log10P, fill = inside)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, notch = T) +
  geom_jitter(alpha = 0.3, width = 0.1, color = "black") +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", comparisons = list(c(1, 2)), label.y = 2) +
  scale_fill_discrete(guide = F) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank()) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(x = "Mutation in kinase domain", y = "Kinase activity (log10P)")

#+ fig.width=2, fig.height=4
ka_mut


#' inside lof mutation vs background
ka_mut <- ka_dbt %>%
  left_join(mutKin, by = c("kinase", "sample")) %>%
  mutate(inside = replace_na(inside, "background"))
  #mutate(inside = str_replace(inside, "FALSE", "background"))

ka_mut <- ka_mut %>%
  ggplot(mapping = aes(x = inside, y = log10P, fill = inside)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, notch = T) +
  geom_jitter(data = ka_mut %>% mutate(log10P = ifelse(inside == "background", NA, log10P)), alpha = 0.3, width = 0.1, show.legend = F) +
  #geom_violin() +
  #geom_boxplot(width=0.1, fill = "grey", outlier.size = 1) +
  theme_classic() +
  stat_compare_means(method = "wilcox.test", comparisons = list(c(1, 3), c(2, 3)), label.y = c(2, 1.5), tip.length = 0.005) +
  scale_fill_discrete(guide = F) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank()) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(x = "Mutation in kinase domain", y = "Kinase activity (log10P)")

#+ fig.width=2, fig.height=4
ka_mut

ggsave(filename = "prot_reg_kinase_activity_bk_ins_out_kinase_domain.png", plot = ka_mut, path = "./output/plots/mutation_impact_kin_activity/", height = 4, width = 2)
unlink("prot_reg_kinase_activity_bk_ins_out_kinase_domain.png")



#' ## stratified by mutation type
mutKin <- mut %>%
  select(sample, variant_class, gene = gene_symbol, pos = prot_pos) %>%
  distinct() %>%
  filter(!is.na(pos)) %>%
  # remove cases of samples with multiple mutations on the same gene
  #group_by(sample, gene) %>%
  #mutate(n = n()) %>%
  #ungroup() %>%
  #filter(n == 1) %>%
  #select(-n) %>%
  inner_join(kinDomains, by = c("gene" = "kinase")) %>%
  mutate(inside = map2_lgl(.x = pos, .y = domains, .f = ~ .y %>% mutate(pos = .x, inside = (pos >= start & pos <= end)) %>% pull(inside) %>% any())) %>%
  select(-domains) %>%
  group_by(variant_class) %>%
  filter(any(inside)) %>%
  ungroup() %>%
  group_by(variant_class, sample, gene) %>%
  summarise(inside = any(inside)) %>%
  ungroup() %>%
  mutate_if(is_logical, as.character) %>%
  rename(kinase = gene) %>%
  filter(!variant_class %in% c("Splice_Site", "In_Frame_Ins"))


labels <- mutKin %>%
  inner_join(ka_dbt, by = c("kinase", "sample")) %>%
  group_by(variant_class, inside) %>%
  summarise(n = n(), median = median(log10P)) %>%
  ungroup() %>%
  mutate(label = str_c(n, round(median, 3), sep = ": "))


#' inside vs outside lof mutation
ka_mut <- mutKin %>%
  inner_join(ka_dbt, by = c("kinase", "sample")) %>%
  ggplot(mapping = aes(x = inside, y = log10P, fill = inside)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(alpha = 0.4, width = 0.1, color = "black") +
  #geom_text(data = labels, mapping = aes(x = inside, y = 2, label = label)) +
  theme_classic() +
  facet_wrap(~ variant_class, scales = "free", nrow = 1) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c(1, 2))) +
  scale_fill_discrete(guide = F) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank()) +
  #scale_y_continuous(limits = c(-0.2, 0.2)) +
  labs(x = "Mutation inside kinase domain", y = "log10P")

#+ fig.width=8, fig.height=4
ka_mut


#' inside lof mutation vs background
ka_mut <- mutKin %>%
  group_by(variant_class) %>%
  nest() %>%
  ungroup() %>%
  mutate(data2 = map(.x = data, .f = ~ .x %>% right_join(ka_dbt, by = c("kinase", "sample")) %>% mutate(inside = replace_na(inside, "background"), inside = str_replace(inside, "FALSE", "background")))) %>%
  select(-data) %>%
  unnest() %>%
  ggplot(mapping = aes(x = log10P, fill = inside, color = inside)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  facet_wrap(~ variant_class) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank()) +
  labs(x = "log10P", y = "Density")

#+ fig.width=8, fig.height=2
ka_mut

