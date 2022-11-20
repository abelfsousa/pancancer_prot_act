#' ---
#' title: "Correlation of CPTAC to TCGA RPPA protein/phosphorylation data"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=6, fig.height=6)


#' load R packages
library(tidyverse)
library(ggpubr)
library(ggbeeswarm)

source("./src/utils/matchTool.R")

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load CPTAC protein/phosphorylation data
phospho <- read_tsv(file = "./output/files/phosphoproteomics.txt.gz") %>%
  pivot_longer(-c(gene, psite, psites), names_to = "sample", values_to = "cptac_expr") %>%
  rename(protein_cptac = gene, psite_cptac = psite) %>%
  filter(!is.na(cptac_expr)) %>%
  #filter(psites == 1) %>%
  select(-psites)
  
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "cptac_expr") %>%
  rename(protein_cptac = gene) %>%
  filter(!is.na(cptac_expr))


#' load RPPA rotein/phosphorylation data
load("./data/protein/tcga/rppa/pan_can_tcga_rppa_mat_updated.Rdata")
pan_can_tcga_rppa_mat <- pan_can_tcga_rppa_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  as_tibble() %>%
  select_at(.vars = c(1, which(str_sub(colnames(.),14,15) == "01"))) %>%
  pivot_longer(-feature, names_to = "sample", values_to = "expr") %>%
  mutate(sample = str_replace_all(str_sub(sample, 1, 12), "-", ".")) %>%
  filter(sample %in% unique(protein$sample))

load("./data/protein/tcga/rppa/matched_prot_phosphoprot_list_updated.Rdata")
matched_prot_phosphoprot_list <- matched_prot_phosphoprot_list %>%
  as_tibble()

rppa_pho <- pan_can_tcga_rppa_mat %>%
  inner_join(matched_prot_phosphoprot_list[, c("fil_rppa_phospho_ab", "fil_gene_prot_list")], by = c("feature" = "fil_rppa_phospho_ab")) %>%
  select(psite = feature, protein = fil_gene_prot_list, sample, rppa_expr = expr) %>%
  mutate(psite = map_chr(.x = psite, .f = ~ str_split_fixed(.x, "_p", 2)[,2])) %>%
  mutate(psite = paste(protein, psite, sep = "_")) %>%
  mutate(psite = map_chr(.x=psite, .f = ~ if(str_count(.x, "\\_") == 2) {strReverse(str_replace(strReverse(.x), "_", ""))} else{.x})) %>% 
  select(protein, psite, everything()) %>%
  rename(protein_rppa = protein, psite_rppa = psite)

rppa_prot <- pan_can_tcga_rppa_mat %>%
  inner_join(distinct(matched_prot_phosphoprot_list[, c("fil_rppa_prot_ab", "fil_gene_prot_list")]), by = c("feature" = "fil_rppa_prot_ab")) %>%
  select(protein = fil_gene_prot_list, sample, rppa_expr = expr) %>%
  rename(protein_rppa = protein)


#' correlate CPTAC to RPPA phosphosite quantifications on the same phosphosite
cptac_rppa_psite_corr_same <- phospho %>%
  inner_join(rppa_pho, by = c("protein_cptac" = "protein_rppa", "psite_cptac" = "psite_rppa", "sample")) %>%
  group_by(protein_cptac, psite_cptac) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) >= 5)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$cptac_expr, .x$rppa_expr, method = "pearson")))) %>%
  select(-data) %>%
  unnest() %>%
  mutate(protein_rppa = protein_cptac, psite_rppa = psite_cptac) %>%
  select(protein_cptac, psite_cptac, protein_rppa, psite_rppa, corr = estimate, p_value = p.value) %>%
  mutate(class = "Same")


#' correlate CPTAC to RPPA phosphosite quantifications between different phosphosites
psites <- cptac_rppa_psite_corr_same$psite_cptac

cptac_rppa_psite_corr_all <- phospho %>%
  filter(psite_cptac %in% psites) %>%
  inner_join(rppa_pho %>% filter(psite_rppa %in% psites), by = c("sample")) %>%
  select(sample, protein_cptac, psite_cptac, protein_rppa, psite_rppa, cptac_expr, rppa_expr) %>%
  group_by(protein_cptac, psite_cptac, protein_rppa, psite_rppa) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$cptac_expr, .x$rppa_expr, method = "pearson")))) %>%
  select(-data) %>%
  unnest() %>%
  select(protein_cptac, psite_cptac, protein_rppa, psite_rppa, corr = estimate, p_value = p.value) %>%
  filter(psite_rppa != psite_cptac) %>%
  mutate(class = "Others")


#' plot the distributions
# cptac_rppa_psite_corr_plot <- cptac_rppa_psite_corr_same %>%
#   bind_rows(cptac_rppa_psite_corr_all) %>%
#   ggplot() +
#   geom_boxplot(mapping = aes(x = class, y = corr, fill = class), alpha = 0.8, outlier.shape = NA, notch = F) +
#   geom_jitter(mapping = aes(x = class, y = corr, alpha = class), width = 0.1, color = "black") +
#   stat_compare_means(mapping = aes(x = class, y = corr), method = "t.test", label.y = -0.7, label.x = 2.4, size = 5) +
#   scale_alpha_manual(values = c(0.3, 0.5), guide = FALSE) +
#   scale_fill_manual(values = c("#bdbdbd", "#67a9cf"), guide = F) +
#   coord_flip() +
#   theme_classic() +
#   theme(
#     #plot.margin = unit(c(0,1,0,0), "cm"),
#     axis.title.x = element_text(colour = "black", size = 17),
#     axis.title.y = element_text(colour = "black", size = 17, hjust = 1),
#     axis.text = element_text(colour = "black", size = 16)) +
#   labs(x = "Phosphosite pair", y = "Pearson's r (CPTAC vs RPPA phosphosite)")
# 
# #+ fig.width=6, fig.height=2
# cptac_rppa_psite_corr_plot
# 
# ggsave(filename = "cptac_rppa_psite_cor.png", plot = cptac_rppa_psite_corr_plot, path = "./output/plots/cptac_rppa_phospho_cor/", height = 2, width = 6)
# ggsave(filename = "cptac_rppa_psite_cor.pdf", plot = cptac_rppa_psite_corr_plot, path = "./output/plots/cptac_rppa_phospho_cor/", height = 2, width = 6)
# unlink("cptac_rppa_psite_cor.png")
# unlink("cptac_rppa_psite_cor.pdf")


#' plot the distributions
cptac_rppa_psite_corr_plot <- cptac_rppa_psite_corr_same %>%
  bind_rows(cptac_rppa_psite_corr_all) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = class, y = corr), outlier.shape = NA, color = "black", fill = "#f0f0f0") +
  geom_beeswarm(mapping = aes(x = class, y = corr, color = class, alpha = class)) +
  stat_compare_means(mapping = aes(x = class, y = corr), method = "wilcox.test", label.y = 1, label.x = 0.9, size = 5) +
  scale_alpha_manual(values = c(0.4, 0.8), guide = FALSE) +
  scale_color_manual(values = c("#525252", "#e31a1c"), guide = F) +
  #scale_color_manual(values = c("#a6611a", "#018571"), guide = F) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 18),
    axis.text = element_text(colour = "black", size = 16)) +
  labs(x = "Phosphosite pair", y = "Pearson's r (CPTAC vs RPPA phosphosite)")

#+ fig.width=3, fig.height=6
cptac_rppa_psite_corr_plot

ggsave(filename = "cptac_rppa_psite_cor.png", plot = cptac_rppa_psite_corr_plot, path = "./output/plots/cptac_rppa_phospho_cor/", height = 6, width = 3)
ggsave(filename = "cptac_rppa_psite_cor.pdf", plot = cptac_rppa_psite_corr_plot, path = "./output/plots/cptac_rppa_phospho_cor/", height = 6, width = 3)


#' plot the distributions
cptac_rppa_psite_corr_plot <- cptac_rppa_psite_corr_same %>%
  bind_rows(cptac_rppa_psite_corr_all) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = class, y = corr, fill = class), outlier.shape = NA, alpha = 0.8) +
  geom_beeswarm(mapping = aes(x = class, y = corr, alpha = class), color = "black") +
  stat_compare_means(mapping = aes(x = class, y = corr, label = paste0("p = ", ..p.format..)), method = "wilcox.test", label.y = 1, label.x = 0.9, size = 5) +
  scale_alpha_manual(values = c(0.2, 0.8), guide = FALSE) +
  #scale_fill_manual(values = c("#525252", "#e31a1c")) +
  scale_fill_manual(values = c("#bdbdbd", "#e31a1c")) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 20),
    axis.text = element_text(colour = "black", size = 18),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 18)) +
  labs(x = "Phosphosite pair", y = "Pearson's r (CPTAC vs RPPA phosphosite)")
  #guides(fill = guide_legend(override.aes = list(alpha=1)))

#+ fig.width=2, fig.height=6
cptac_rppa_psite_corr_plot

ggsave(filename = "cptac_rppa_psite_cor.png", plot = cptac_rppa_psite_corr_plot, path = "./output/plots/cptac_rppa_phospho_cor/", height = 6, width = 3)
ggsave(filename = "cptac_rppa_psite_cor.pdf", plot = cptac_rppa_psite_corr_plot, path = "./output/plots/cptac_rppa_phospho_cor/", height = 6, width = 3)


#' correlate CPTAC to RPPA protein quantifications on the same protein
cptac_rppa_protein_corr_same <- protein %>%
  inner_join(rppa_prot, by = c("protein_cptac" = "protein_rppa", "sample")) %>%
  group_by(protein_cptac) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$cptac_expr, .x$rppa_expr, method = "pearson")))) %>%
  select(-data) %>%
  unnest() %>%
  mutate(protein_rppa = protein_cptac) %>%
  select(protein_cptac, protein_rppa, corr = estimate, p_value = p.value) %>%
  mutate(class = "Same")


#' correlate CPTAC to RPPA protein quantifications between different proteins
proteins <- cptac_rppa_protein_corr_same$protein_cptac

cptac_rppa_protein_corr_all <- protein %>%
  filter(protein_cptac %in% proteins) %>%
  inner_join(rppa_prot %>% filter(protein_rppa %in% proteins), by = c("sample")) %>%
  select(sample, protein_cptac, protein_rppa, cptac_expr, rppa_expr) %>%
  group_by(protein_cptac, protein_rppa) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$cptac_expr, .x$rppa_expr, method = "pearson")))) %>%
  select(-data) %>%
  unnest() %>%
  select(protein_cptac, protein_rppa, corr = estimate, p_value = p.value) %>%
  filter(protein_cptac != protein_rppa) %>%
  mutate(class = "Others")


#' plot the distributions
cptac_rppa_protein_corr_plot <- cptac_rppa_protein_corr_same %>%
  bind_rows(cptac_rppa_protein_corr_all) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = class, y = corr, fill = class), alpha = 0.8, outlier.shape = NA, notch = F) +
  geom_jitter(mapping = aes(x = class, y = corr, alpha = class), width = 0.1, color = "black") +
  stat_compare_means(mapping = aes(x = class, y = corr), method = "t.test", label.y = 0.9, size = 3.3) +
  scale_alpha_manual(values = c(0.05, 0.5), guide = FALSE) +
  scale_fill_manual(values = c("#bdbdbd", "#67a9cf"), guide = F) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 15),
    axis.text = element_text(colour = "black", size = 13)) +
  labs(x = "Protein pair", y = "Pearson's r (CPTAC vs RPPA protein quantifications)")

#+ fig.width=2, fig.height=6
cptac_rppa_protein_corr_plot

ggsave(filename = "cptac_rppa_prot_cor.png", plot = cptac_rppa_protein_corr_plot, path = "./output/plots/cptac_rppa_phospho_cor/", height = 6, width = 2)
ggsave(filename = "cptac_rppa_prot_cor.pdf", plot = cptac_rppa_protein_corr_plot, path = "./output/plots/cptac_rppa_phospho_cor/", height = 6, width = 2)
unlink("cptac_rppa_prot_cor.png")
unlink("cptac_rppa_prot_cor.pdf")

