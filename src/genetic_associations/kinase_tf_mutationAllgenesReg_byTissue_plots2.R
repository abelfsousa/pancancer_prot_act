#

# load R packages
library(tidyverse)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load mutation to kinase/TF activity associations
kin_tissue <- read_tsv(file = str_c("./output/files/ka_mutStatus_allGenes_byTissue", 2, ".txt.gz"), progress = T, col_types = cols())
tf_tissue <- read_tsv(file = str_c("./output/files/tf_mutStatus_allGenes_byTissue", 4, ".txt.gz"), progress = T, col_types = cols())

kin_all <- read_tsv(file = "./output/files/ka_mutStatus_allGenes.txt.gz", progress = T, col_types = cols())
tf_all <- read_tsv(file = "./output/files/tf_mutStatus_allGenes.txt.gz", progress = T, col_types = cols())


# compare all associations 
kin_associations <- kin_all %>%
  rename(estimate_all = estimate, p_value_all = p.value, p_adjust_all = p.adjust) %>%
  inner_join(kin_tissue %>% rename(estimate_tissue = estimate, p_value_tissue = p.value, p_adjust_tissue = p.adjust), by = c("kinase", "gene")) %>%
  mutate(across(.cols = starts_with("p_"), .fns = ~ -log10(.x)))

kin_associations_tidy <- kin_associations %>%
  select(kinase, gene, tissue, batch, everything()) %>%
  pivot_longer(cols = -c(kinase, gene, tissue, batch), names_to = "stat", values_to = "value")

kin_associations2 <- tibble(stat1 = c("estimate_all", "p_value_all", "p_adjust_all"), stat2 = c("estimate_tissue", "p_value_tissue", "p_adjust_tissue")) %>%
  inner_join(kin_associations_tidy, by = c("stat1" = "stat")) %>%
  rename(value1 = value) %>%
  inner_join(kin_associations_tidy, by = c("stat2" = "stat", "kinase", "gene", "tissue", "batch")) %>%
  rename(value2 = value) %>%
  unite(col = "comparison", "stat1", "stat2", sep = " vs ") %>%
  rename(protein = kinase) %>%
  mutate(protein_type = "kinase") %>%
  select(protein_type, everything())

tf_associations <- tf_all %>%
  rename(estimate_all = estimate, p_value_all = p.value, p_adjust_all = p.adjust) %>%
  inner_join(tf_tissue %>% rename(estimate_tissue = estimate, p_value_tissue = p.value, p_adjust_tissue = p.adjust), by = c("tf", "gene")) %>%
  mutate(across(.cols = starts_with("p_"), .fns = ~ -log10(.x)))

tf_associations_tidy <- tf_associations %>%
  select(tf, gene, tissue, batch, everything()) %>%
  pivot_longer(cols = -c(tf, gene, tissue, batch), names_to = "stat", values_to = "value")

tf_associations2 <- tibble(stat1 = c("estimate_all", "p_value_all", "p_adjust_all"), stat2 = c("estimate_tissue", "p_value_tissue", "p_adjust_tissue")) %>%
  inner_join(tf_associations_tidy, by = c("stat1" = "stat")) %>%
  rename(value1 = value) %>%
  inner_join(tf_associations_tidy, by = c("stat2" = "stat", "tf", "gene", "tissue", "batch")) %>%
  rename(value2 = value) %>%
  unite(col = "comparison", "stat1", "stat2", sep = " vs ") %>%
  rename(protein = tf) %>%
  mutate(protein_type = "tf") %>%
  select(protein_type, everything())

kin_plot <- kin_associations2 %>%
  mutate(comparison = fct_relevel(comparison, "estimate_all vs estimate_tissue", "p_value_all vs p_value_tissue")) %>%
  mutate(comparison = fct_recode(comparison, `Effect size` = "estimate_all vs estimate_tissue", `P-value (-log10)` = "p_value_all vs p_value_tissue", `FDR (-log10)` = "p_adjust_all vs p_adjust_tissue")) %>%
  ggplot(mapping = aes(x = value1, y = value2)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(comparison), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 14)) +
  labs(x = "All tissues", y = "By tissue", title = "Kinases")
ggsave("corr_stats_all_tissues_by_tissue_all_assoc_kinases.png", plot = kin_plot, path = "./output/plots/genetic_associations/", width = 7, height = 3)
#ggsave("corr_stats_all_tissues_by_tissue_all_assoc_kinases.pdf", plot = kin_plot, path = "./output/plots/genetic_associations/", width = 7, height = 3)

kin_plot <- kin_associations2 %>%
  mutate(comparison = fct_relevel(comparison, "estimate_all vs estimate_tissue", "p_value_all vs p_value_tissue")) %>%
  mutate(comparison = fct_recode(comparison, `Effect size` = "estimate_all vs estimate_tissue", `P-value (-log10)` = "p_value_all vs p_value_tissue", `FDR (-log10)` = "p_adjust_all vs p_adjust_tissue")) %>%
  filter(comparison == "Effect size") %>%
  ggplot(mapping = aes(x = value1, y = value2)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 14)) +
  labs(x = "Beta (all tissues)", y = "Beta (by tissue)", title = "Kinases")
ggsave("corr_stats_all_tissues_by_tissue_all_assoc_kinases_beta.png", plot = kin_plot, path = "./output/plots/genetic_associations/", width = 3, height = 3)

kin_plot <- kin_associations2 %>%
  mutate(dataset = str_c(tissue, batch, sep = "_")) %>%
  mutate(comparison = fct_relevel(comparison, "estimate_all vs estimate_tissue", "p_value_all vs p_value_tissue")) %>%
  mutate(comparison = fct_recode(comparison, `Effect size` = "estimate_all vs estimate_tissue", `P-value (-log10)` = "p_value_all vs p_value_tissue", `FDR (-log10)` = "p_adjust_all vs p_adjust_tissue")) %>%
  ggplot(mapping = aes(x = value1, y = value2)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_grid(rows = vars(dataset), cols = vars(comparison), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 14)) +
  labs(x = "All tissues", y = "By tissue", title = "Kinases")

kin_plot <- kin_associations2 %>%
  mutate(dataset = str_c(tissue, batch, sep = "_")) %>%
  mutate(comparison = fct_relevel(comparison, "estimate_all vs estimate_tissue", "p_value_all vs p_value_tissue")) %>%
  mutate(comparison = fct_recode(comparison, `Effect size` = "estimate_all vs estimate_tissue", `P-value (-log10)` = "p_value_all vs p_value_tissue", `FDR (-log10)` = "p_adjust_all vs p_adjust_tissue")) %>%
  filter(comparison == "Effect size") %>%
  ggplot(mapping = aes(x = value1, y = value2)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(dataset), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 14)) +
  labs(x = "All tissues", y = "By tissue", title = "Kinases")
ggsave("corr_stats_all_tissues_by_tissue_all_assoc_kinases_BY_TISSUE_beta.png", plot = kin_plot, path = "./output/plots/genetic_associations/", width = 7, height = 7)

kin_plot <- kin_associations2 %>%
  mutate(dataset = str_c(tissue, batch, sep = "_")) %>%
  mutate(comparison = fct_relevel(comparison, "estimate_all vs estimate_tissue", "p_value_all vs p_value_tissue")) %>%
  mutate(comparison = fct_recode(comparison, `Effect size` = "estimate_all vs estimate_tissue", `P-value (-log10)` = "p_value_all vs p_value_tissue", `FDR (-log10)` = "p_adjust_all vs p_adjust_tissue")) %>%
  filter(comparison == "P-value (-log10)") %>%
  ggplot(mapping = aes(x = value1, y = value2)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(dataset), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 14)) +
  labs(x = "All tissues", y = "By tissue", title = "Kinases")
ggsave("corr_stats_all_tissues_by_tissue_all_assoc_kinases_BY_TISSUE_pval.png", plot = kin_plot, path = "./output/plots/genetic_associations/", width = 7, height = 7)

kin_plot <- kin_associations2 %>%
  mutate(dataset = str_c(tissue, batch, sep = "_")) %>%
  mutate(comparison = fct_relevel(comparison, "estimate_all vs estimate_tissue", "p_value_all vs p_value_tissue")) %>%
  mutate(comparison = fct_recode(comparison, `Effect size` = "estimate_all vs estimate_tissue", `P-value (-log10)` = "p_value_all vs p_value_tissue", `FDR (-log10)` = "p_adjust_all vs p_adjust_tissue")) %>%
  filter(comparison == "FDR (-log10)") %>%
  ggplot(mapping = aes(x = value1, y = value2)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(dataset), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 14)) +
  labs(x = "All tissues", y = "By tissue", title = "Kinases")
ggsave("corr_stats_all_tissues_by_tissue_all_assoc_kinases_BY_TISSUE_fdr.png", plot = kin_plot, path = "./output/plots/genetic_associations/", width = 7, height = 7)

kin_plot <- kin_associations %>%
  ggplot(mapping = aes(x = estimate_tissue, y = estimate_all)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  #ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  theme_classic() +
  labs(x = "All tissues", y = "By tissue")

tf_plot <- tf_associations2 %>%
  mutate(comparison = fct_relevel(comparison, "estimate_all vs estimate_tissue", "p_value_all vs p_value_tissue")) %>%
  mutate(comparison = fct_recode(comparison, `Effect size` = "estimate_all vs estimate_tissue", `P-value (-log10)` = "p_value_all vs p_value_tissue", `FDR (-log10)` = "p_adjust_all vs p_adjust_tissue")) %>%
  ggplot(mapping = aes(x = value1, y = value2)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(comparison), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 14)) +
  labs(x = "All tissues", y = "By tissue", title = "TFs")
ggsave("corr_stats_all_tissues_by_tissue_all_assoc_tfs.png", plot = tf_plot, path = "./output/plots/genetic_associations/", width = 7, height = 3)
#ggsave("corr_stats_all_tissues_by_tissue_all_assoc_tfs.pdf", plot = tf_plot, path = "./output/plots/genetic_associations/", width = 7, height = 3)

tf_plot <- tf_associations2 %>%
  mutate(comparison = fct_relevel(comparison, "estimate_all vs estimate_tissue", "p_value_all vs p_value_tissue")) %>%
  mutate(comparison = fct_recode(comparison, `Effect size` = "estimate_all vs estimate_tissue", `P-value (-log10)` = "p_value_all vs p_value_tissue", `FDR (-log10)` = "p_adjust_all vs p_adjust_tissue")) %>%
  filter(comparison == "Effect size") %>%
  ggplot(mapping = aes(x = value1, y = value2)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 14)) +
  labs(x = "Beta (all tissues)", y = "Beta (by tissue)", title = "TFs")
ggsave("corr_stats_all_tissues_by_tissue_all_assoc_tfs_beta.png", plot = tf_plot, path = "./output/plots/genetic_associations/", width = 3, height = 3)

tf_plot <- tf_associations2 %>%
  mutate(dataset = str_c(tissue, batch, sep = "_")) %>%
  mutate(comparison = fct_relevel(comparison, "estimate_all vs estimate_tissue", "p_value_all vs p_value_tissue")) %>%
  mutate(comparison = fct_recode(comparison, `Effect size` = "estimate_all vs estimate_tissue", `P-value (-log10)` = "p_value_all vs p_value_tissue", `FDR (-log10)` = "p_adjust_all vs p_adjust_tissue")) %>%
  ggplot(mapping = aes(x = value1, y = value2)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_grid(rows = vars(dataset), cols = vars(comparison), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 14)) +
  labs(x = "All tissues", y = "By tissue", title = "TFs")

tf_plot <- tf_associations2 %>%
  mutate(dataset = str_c(tissue, batch, sep = "_")) %>%
  mutate(comparison = fct_relevel(comparison, "estimate_all vs estimate_tissue", "p_value_all vs p_value_tissue")) %>%
  mutate(comparison = fct_recode(comparison, `Effect size` = "estimate_all vs estimate_tissue", `P-value (-log10)` = "p_value_all vs p_value_tissue", `FDR (-log10)` = "p_adjust_all vs p_adjust_tissue")) %>%
  filter(comparison == "Effect size") %>%
  ggplot(mapping = aes(x = value1, y = value2)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(dataset), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 14)) +
  labs(x = "All tissues", y = "By tissue", title = "TFs")
ggsave("corr_stats_all_tissues_by_tissue_all_assoc_tfs_BY_TISSUE_beta.png", plot = tf_plot, path = "./output/plots/genetic_associations/", width = 7, height = 7)

tf_plot <- tf_associations2 %>%
  mutate(dataset = str_c(tissue, batch, sep = "_")) %>%
  mutate(comparison = fct_relevel(comparison, "estimate_all vs estimate_tissue", "p_value_all vs p_value_tissue")) %>%
  mutate(comparison = fct_recode(comparison, `Effect size` = "estimate_all vs estimate_tissue", `P-value (-log10)` = "p_value_all vs p_value_tissue", `FDR (-log10)` = "p_adjust_all vs p_adjust_tissue")) %>%
  filter(comparison == "P-value (-log10)") %>%
  ggplot(mapping = aes(x = value1, y = value2)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(dataset), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 14)) +
  labs(x = "All tissues", y = "By tissue", title = "TFs")
ggsave("corr_stats_all_tissues_by_tissue_all_assoc_tfs_BY_TISSUE_pval.png", plot = tf_plot, path = "./output/plots/genetic_associations/", width = 7, height = 7)

tf_plot <- tf_associations2 %>%
  mutate(dataset = str_c(tissue, batch, sep = "_")) %>%
  mutate(comparison = fct_relevel(comparison, "estimate_all vs estimate_tissue", "p_value_all vs p_value_tissue")) %>%
  mutate(comparison = fct_recode(comparison, `Effect size` = "estimate_all vs estimate_tissue", `P-value (-log10)` = "p_value_all vs p_value_tissue", `FDR (-log10)` = "p_adjust_all vs p_adjust_tissue")) %>%
  filter(comparison == "FDR (-log10)") %>%
  ggplot(mapping = aes(x = value1, y = value2)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(dataset), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 14)) +
  labs(x = "All tissues", y = "By tissue", title = "TFs")
ggsave("corr_stats_all_tissues_by_tissue_all_assoc_tfs_BY_TISSUE_fdr.png", plot = tf_plot, path = "./output/plots/genetic_associations/", width = 7, height = 7)

tf_plot <- tf_associations %>%
  ggplot(mapping = aes(x = estimate_tissue, y = estimate_all)) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  #ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  theme_classic() +
  labs(x = "All tissues", y = "By tissue")


# compare the significant associations
kin_associations <- kin_all %>%
  filter(p.adjust < 0.05) %>%
  rename(estimate_all = estimate, p_value_all = p.value, p_adjust_all = p.adjust) %>%
  inner_join(kin_tissue %>% filter(p.adjust < 0.05) %>% rename(estimate_tissue = estimate, p_value_tissue = p.value, p_adjust_tissue = p.adjust), by = c("kinase", "gene")) %>%
  mutate(across(.cols = starts_with("p_"), .fns = ~ -log10(.x)))

kin_associations_tidy <- kin_associations %>%
  select(kinase, gene, tissue, batch, everything()) %>%
  pivot_longer(cols = -c(kinase, gene, tissue, batch), names_to = "stat", values_to = "value")

kin_associations2 <- tibble(stat1 = c("estimate_all", "p_value_all", "p_adjust_all"), stat2 = c("estimate_tissue", "p_value_tissue", "p_adjust_tissue")) %>%
  inner_join(kin_associations_tidy, by = c("stat1" = "stat")) %>%
  rename(value1 = value) %>%
  inner_join(kin_associations_tidy, by = c("stat2" = "stat", "kinase", "gene", "tissue", "batch")) %>%
  rename(value2 = value) %>%
  unite(col = "comparison", "stat1", "stat2", sep = " vs ") %>%
  rename(protein = kinase) %>%
  mutate(protein_type = "kinase") %>%
  select(protein_type, everything())

tf_associations <- tf_all %>%
  filter(p.adjust < 0.05) %>%
  rename(estimate_all = estimate, p_value_all = p.value, p_adjust_all = p.adjust) %>%
  inner_join(tf_tissue %>% filter(p.adjust < 0.05) %>% rename(estimate_tissue = estimate, p_value_tissue = p.value, p_adjust_tissue = p.adjust), by = c("tf", "gene")) %>%
  mutate(across(.cols = starts_with("p_"), .fns = ~ -log10(.x)))

tf_associations_tidy <- tf_associations %>%
  select(tf, gene, tissue, batch, everything()) %>%
  pivot_longer(cols = -c(tf, gene, tissue, batch), names_to = "stat", values_to = "value")

tf_associations2 <- tibble(stat1 = c("estimate_all", "p_value_all", "p_adjust_all"), stat2 = c("estimate_tissue", "p_value_tissue", "p_adjust_tissue")) %>%
  inner_join(tf_associations_tidy, by = c("stat1" = "stat")) %>%
  rename(value1 = value) %>%
  inner_join(tf_associations_tidy, by = c("stat2" = "stat", "tf", "gene", "tissue", "batch")) %>%
  rename(value2 = value) %>%
  unite(col = "comparison", "stat1", "stat2", sep = " vs ") %>%
  rename(protein = tf) %>%
  mutate(protein_type = "tf") %>%
  select(protein_type, everything())

associations <- bind_rows(kin_associations2, tf_associations2)

associations %>%
  ggplot(mapping = aes(x = value1, y = value2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_grid(rows = vars(comparison), cols = vars(protein_type), scales = "free", space = "fixed") +
  theme_classic()

kin_plot <- kin_associations2 %>%
  mutate(comparison = fct_relevel(comparison, "estimate_all vs estimate_tissue", "p_value_all vs p_value_tissue")) %>%
  mutate(comparison = fct_recode(comparison, `Effect size` = "estimate_all vs estimate_tissue", `P-value (-log10)` = "p_value_all vs p_value_tissue", `FDR (-log10)` = "p_adjust_all vs p_adjust_tissue")) %>%
  ggplot(mapping = aes(x = value1, y = value2)) +
  geom_point(alpha = 0.5, size = 3) +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(comparison), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank()) +
  labs(x = "All tissues", y = "By tissue")
ggsave("corr_stats_all_tissues_by_tissue_signf_assoc_kinases.png", plot = kin_plot, path = "./output/plots/genetic_associations/", width = 7, height = 3)
ggsave("corr_stats_all_tissues_by_tissue_signf_assoc_kinases.pdf", plot = kin_plot, path = "./output/plots/genetic_associations/", width = 7, height = 3)

tf_plot <- tf_associations2 %>%
  mutate(comparison = fct_relevel(comparison, "estimate_all vs estimate_tissue", "p_value_all vs p_value_tissue")) %>%
  mutate(comparison = fct_recode(comparison, `Effect size` = "estimate_all vs estimate_tissue", `P-value (-log10)` = "p_value_all vs p_value_tissue", `FDR (-log10)` = "p_adjust_all vs p_adjust_tissue")) %>%
  ggplot(mapping = aes(x = value1, y = value2)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = F) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(comparison), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank()) +
  labs(x = "All tissues", y = "By tissue")
ggsave("corr_stats_all_tissues_by_tissue_signf_assoc_tfs.png", plot = tf_plot, path = "./output/plots/genetic_associations/", width = 7, height = 3)
ggsave("corr_stats_all_tissues_by_tissue_signf_assoc_tfs.pdf", plot = tf_plot, path = "./output/plots/genetic_associations/", width = 7, height = 3)

