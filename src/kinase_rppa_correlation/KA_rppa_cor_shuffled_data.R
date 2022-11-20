#' ---
#' title: "Correlation of kinase activity to RPPA phosphorylation data"
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

set.seed(123)

#' load kinase-activity inference data
k_subN <- 3
k_psit <- 3

ka_protreg_dbt <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-n, -source_type)

ka_nonreg_dbt <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-n, -source_type)

ka_kin_psite_protreg <- read_tsv(file = "./output/files/CPTAC_KA_kin_psites_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(n >= k_psit) %>%
  select(-n) %>%
  rename(kinase=gene)

shuffle_tibble <- function(tib){
  tib_shuf <- tib %>%
    mutate(index = sample(seq_len(nrow(tib)), nrow(tib), replace = FALSE)) %>%
    mutate(sample = sample[index]) %>%
    select(-index)
}

ka_protreg_dbt_shuff <- ka_protreg_dbt %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(data2 = map(.x=data, .f=shuffle_tibble)) %>%
  select(-data) %>%
  unnest()

ka_nonreg_dbt_shuff <- ka_nonreg_dbt %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(data2 = map(.x=data, .f=shuffle_tibble)) %>%
  select(-data) %>%
  unnest()

ka_kin_psite_protreg_shuff <- ka_kin_psite_protreg %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(data2 = map(.x=data, .f=shuffle_tibble)) %>%
  select(-data) %>%
  unnest()


#' load RPPA phosphorylation data
load("./data/protein/tcga/rppa/pan_can_tcga_rppa_mat_updated.Rdata")
pan_can_tcga_rppa_mat <- pan_can_tcga_rppa_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  as_tibble() %>%
  select_at(.vars = c(1, which(str_sub(colnames(.),14,15) == "01"))) %>%
  pivot_longer(-feature, names_to = "sample", values_to = "expr") %>%
  mutate(sample = str_replace_all(str_sub(sample, 1, 12), "-", ".")) %>%
  filter(sample %in% unique(ka_protreg_dbt$sample))

load("./data/protein/tcga/rppa/matched_prot_phosphoprot_list_updated.Rdata")
matched_prot_phosphoprot_list <- matched_prot_phosphoprot_list %>%
  as_tibble()

load("./data/protein/tcga/rppa/kin_ps_relation.Rdata")
kin_ps_relation <- kin_ps_relation %>%
  as_tibble() %>%
  select(kinase=Protein, psite=Ab, activity=Activity)


rppa_pho <- pan_can_tcga_rppa_mat %>%
  inner_join(matched_prot_phosphoprot_list[, c("fil_rppa_phospho_ab", "fil_gene_prot_list")], by = c("feature" = "fil_rppa_phospho_ab")) %>%
  select(psite = feature, protein = fil_gene_prot_list, sample, psite_expr = expr) %>%
  mutate_if(is.factor, as.character)

rppa_prot <- pan_can_tcga_rppa_mat %>%
  inner_join(distinct(matched_prot_phosphoprot_list[, c("fil_rppa_prot_ab", "fil_gene_prot_list")]), by = c("feature" = "fil_rppa_prot_ab")) %>%
  select(protein = fil_gene_prot_list, sample, protein_expr = expr) %>%
  mutate_if(is.factor, as.character)

#shuffle_tibble <- function(tib){
#  tib_shuf <- tib %>%
#    mutate(index = sample(seq_len(nrow(tib)), nrow(tib), replace = FALSE)) %>%
#    mutate(protein = protein[index], psite = psite[index]) %>%
#    select(-index)
#}

rppa_pho_shuff <- rppa_pho %>%
  group_by(psite, protein) %>%
  nest() %>%
  ungroup() %>%
  mutate(data2 = map(.x=data, .f=shuffle_tibble)) %>%
  select(-data) %>%
  unnest()


KA_rppa_pho_cor <- ka_protreg_dbt %>%
  inner_join(rppa_pho, by = c("kinase" = "protein", "sample")) %>%
  select(sample, kinase, psite, log10P, psite_expr) %>%
  group_by(psite, kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$log10P, .x$psite_expr, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(psite, kinase, cor=estimate) %>%
  #inner_join(kin_ps_relation, by = c("kinase", "psite")) %>%
  mutate(class1 = "True kinase", class2="Substrate-based KA")

KA_rppa_pho_cor_shuff <- ka_protreg_dbt_shuff %>%
  inner_join(rppa_pho_shuff, by = c("kinase" = "protein", "sample")) %>%
  select(sample, kinase, psite, log10P, psite_expr) %>%
  group_by(psite, kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$log10P, .x$psite_expr, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(psite, kinase, cor=estimate) %>%
  mutate(class1 = "Shuffled kinase", class2 = "Substrate-based KA")



KA_psites_rppa_pho_cor <- ka_kin_psite_protreg %>%
  inner_join(rppa_pho, by = c("kinase" = "protein", "sample")) %>%
  select(sample, kinase, psite, log10P, psite_expr) %>%
  group_by(psite, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) >= 5)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$log10P, .x$psite_expr, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(psite, kinase, cor=estimate) %>%
  mutate(class1 = "True kinase", class2="RegPsites-based KA")

KA_psites_rppa_pho_cor_shuff <- ka_kin_psite_protreg_shuff %>%
  inner_join(rppa_pho_shuff, by = c("kinase" = "protein", "sample")) %>%
  select(sample, kinase, psite, log10P, psite_expr) %>%
  group_by(psite, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) >= 5)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$log10P, .x$psite_expr, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(psite, kinase, cor=estimate) %>%
  mutate(class1 = "Shuffled kinase", class2 = "RegPsites-based KA")



KA_rppa_pho_cor_plot <- KA_rppa_pho_cor %>%
  bind_rows(KA_rppa_pho_cor_shuff) %>%
  bind_rows(KA_psites_rppa_pho_cor) %>%
  bind_rows(KA_psites_rppa_pho_cor_shuff) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = class1, y = cor), notch = F, outlier.shape = NA) +
  geom_jitter(mapping = aes(x = class1, y = cor, color = class1), alpha = 0.7, width = 0.1) +
  facet_wrap( ~ class2) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 12)) +
  stat_compare_means(mapping = aes(x = class1, y = cor), label.y = -0.1, label.x = 2.5, size=3) +
  coord_flip() +
  scale_color_discrete(guide=F) +
  labs(x = "", y = "Pearson's r (kinase activity vs RPPA phosphosite)")

#+ fig.width=7, fig.height=2
KA_rppa_pho_cor_plot

ggsave(filename = "KA_prot_reg_rppa_phos_cor_distribution.png", plot = KA_rppa_pho_cor_plot, path = "./output/plots/kinase_rppa_correlation/", height = 2, width = 7)
unlink("KA_prot_reg_rppa_phos_cor_distribution.png")



KA_rppa_pho_cor <- ka_nonreg_dbt %>%
  inner_join(rppa_pho, by = c("kinase" = "protein", "sample")) %>%
  select(sample, kinase, psite, log10P, psite_expr) %>%
  group_by(psite, kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$log10P, .x$psite_expr, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(psite, kinase, cor=estimate) %>%
  #inner_join(kin_ps_relation, by = c("kinase", "psite")) %>%
  mutate(class = "True kinase")

KA_rppa_pho_cor_shuff <- ka_nonreg_dbt_shuff %>%
  inner_join(rppa_pho_shuff, by = c("kinase" = "protein", "sample")) %>%
  select(sample, kinase, psite, log10P, psite_expr) %>%
  group_by(psite, kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$log10P, .x$psite_expr, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(psite, kinase, cor=estimate) %>%
  mutate(class = "Shuffled kinase")


KA_rppa_pho_cor_plot <- KA_rppa_pho_cor %>%
  bind_rows(KA_rppa_pho_cor_shuff) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = class, y = cor), notch = TRUE, outlier.shape = NA) +
  geom_jitter(mapping = aes(x = class, y = cor, color = class), alpha = 0.7, width = 0.1) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  stat_compare_means(mapping = aes(x = class, y = cor), label.y = -0.1, label.x = 2.5) +
  coord_flip() +
  theme_classic() +
  scale_color_discrete(guide=F) +
  labs(x = "", y = "Pearson's r (kinase activity vs RPPA phosphosite)")

#+ fig.width=6, fig.height=2
KA_rppa_pho_cor_plot

ggsave(filename = "KA_non_reg_rppa_phos_cor_distribution.png", plot = KA_rppa_pho_cor_plot, path = "./output/plots/kinase_rppa_correlation/", height = 2, width = 6)
unlink("KA_non_reg_rppa_phos_cor_distribution.png")
