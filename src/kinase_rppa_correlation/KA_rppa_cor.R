#' ---
#' title: "Correlation of kinase activity to RPPA phosphorylation data"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=6, fig.height=6)


#' load R packages
library(tidyverse)
library(viridis)
library(ggpubr)
library(Hmisc)
library(ggbeeswarm)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy

set.seed(123)


#' load kinase-activity inference data
k_subN <- 3
k_psit <- 3

ka_protreg <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-n, -source_type)

ka_kinPsite_protreg <- read_tsv(file = "./output/files/CPTAC_KA_kin_psites_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(n >= k_psit) %>%
  select(-n) %>%
  rename(kinase=gene)

ka_nonreg <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-n, -source_type)

ka_kinPsite_nonreg <- read_tsv(file = "./output/files/CPTAC_KA_kin_psites_phosphoQ.txt.gz") %>%
  filter(n >= k_psit) %>%
  select(-n) %>%
  rename(kinase=gene)


#' load RPPA phosphorylation data
load("./data/protein/tcga/rppa/pan_can_tcga_rppa_mat_updated.Rdata")
pan_can_tcga_rppa_mat <- pan_can_tcga_rppa_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  as_tibble() %>%
  select_at(.vars = c(1, which(str_sub(colnames(.),14,15) == "01"))) %>%
  pivot_longer(-feature, names_to = "sample", values_to = "expr") %>%
  mutate(sample = str_replace_all(str_sub(sample, 1, 12), "-", ".")) %>%
  filter(sample %in% unique(ka_protreg$sample))

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


#' correlate RPPA phosphorylation with RPPA protein abundance
rppa_pho_prot <- rppa_pho %>%
  inner_join(rppa_prot, by = c("protein", "sample")) %>%
  group_by(protein, psite) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map_dbl(data, ~ cor.test(.x$psite_expr, .x$protein_expr)$estimate)) %>%
  select(-data) %>%
  ggplot(mapping = aes(x = cor)) +
  geom_density()

#+ fig.width=3, fig.height=2
rppa_pho_prot

ggsave(filename = "rppa_protein_phospho_correlation.png", plot = rppa_pho_prot, path = "./output/plots/kinase_rppa_correlation/", height = 2, width = 3)
ggsave(filename = "rppa_protein_phospho_correlation.pdf", plot = rppa_pho_prot, path = "./output/plots/kinase_rppa_correlation/", height = 2, width = 3)
unlink("rppa_protein_phospho_correlation.png")
unlink("rppa_protein_phospho_correlation.pdf")


#' correlate kinase activities with TCGA RPPA phosphosites on the same kinase
KA_protreg_rppa_pho_cor <- ka_protreg %>%
  inner_join(rppa_pho, by = c("kinase" = "protein", "sample")) %>%
  select(sample, kinase, psite, activity=log10P, psite_expr) %>%
  group_by(psite, kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$psite_expr, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(psite, kinase, pearson_r=estimate, p_value=p.value) %>%
  #inner_join(kin_ps_relation, by = c("kinase", "psite")) %>%
  mutate(class1 = "Kinase-RPPA", class2="Substrate-based", class3 = "Protein reg-out") %>%
  rename(kinaseA=kinase) %>%
  mutate(kinaseB=kinaseA)

KA_protreg_KinPsite_rppa_pho_cor <- ka_kinPsite_protreg %>%
  inner_join(rppa_pho, by = c("kinase" = "protein", "sample")) %>%
  select(sample, kinase, psite, activity=log10P, psite_expr) %>%
  group_by(psite, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) >= 5)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$psite_expr, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(psite, kinase, pearson_r=estimate, p_value=p.value) %>%
  mutate(class1 = "Kinase-RPPA", class2="Psite-based", class3 = "Protein reg-out") %>%
  rename(kinaseA=kinase) %>%
  mutate(kinaseB=kinaseA)

KA_nonreg_rppa_pho_cor <- ka_nonreg %>%
  inner_join(rppa_pho, by = c("kinase" = "protein", "sample")) %>%
  select(sample, kinase, psite, activity=log10P, psite_expr) %>%
  group_by(psite, kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$psite_expr, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(psite, kinase, pearson_r=estimate, p_value=p.value) %>%
  #inner_join(kin_ps_relation, by = c("kinase", "psite")) %>%
  mutate(class1 = "Kinase-RPPA", class2="Substrate-based", class3 = "Protein non reg-out") %>%
  rename(kinaseA=kinase) %>%
  mutate(kinaseB=kinaseA)

KA_nonreg_KinPsite_rppa_pho_cor <- ka_kinPsite_nonreg %>%
  inner_join(rppa_pho, by = c("kinase" = "protein", "sample")) %>%
  select(sample, kinase, psite, activity=log10P, psite_expr) %>%
  group_by(psite, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) >= 5)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$psite_expr, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(psite, kinase, pearson_r=estimate, p_value=p.value) %>%
  #inner_join(kin_ps_relation, by = c("kinase", "psite")) %>%
  mutate(class1 = "Kinase-RPPA", class2="Psite-based", class3 = "Protein non reg-out") %>%
  rename(kinaseA=kinase) %>%
  mutate(kinaseB=kinaseA)


#' correlate kinase activities with all CPTAC kinase phosphosites (all pairs)

# set up a function to convert a squared matrix into a data frame (tibble)
matrix_to_tb <- function(x, y, k, p){
  
  mat <- x
  mat_data <- y
  KIN <- k
  PHO <- p
  
  mat <- mat[p,k]
  df <- mat %>%
    as.data.frame() %>%
    rownames_to_column(var="kinaseA") %>%
    as_tibble() %>%
    pivot_longer(-c("kinaseA"), names_to="kinaseB", values_to=mat_data) %>%
    separate(col = "kinaseA", into = c("kinaseA", "psite"), sep = "-")
  
  df
}

# set up a function to join the rcorr output into a single tibble
joinMatrices <- function(x, kin, pho){
  
  rcorr_matrices <- x
  rcorr_names <- names(rcorr_matrices)
  
  KIN <- kin
  PHO <- pho
  
  dfs <- map2(.x = rcorr_matrices, .y = rcorr_names, .f = matrix_to_tb, k = KIN, p = PHO)
  
  tb <- reduce(.x = dfs, .f = ~ inner_join(.x, .y, by = c("kinaseA", "psite", "kinaseB")))
  
  tb
}


# substrate-based protein reg-out
proteins <- intersect(unique(ka_protreg$kinase), unique(rppa_pho$protein))
samples <- intersect(unique(ka_protreg$sample), unique(rppa_pho$sample))

ka_mat <- ka_protreg %>%
  #semi_join(rppa_pho, by = c("sample")) %>%
  semi_join(rppa_pho, by = c("kinase" = "protein", "sample")) %>%
  arrange(sample) %>%
  pivot_wider(names_from = "sample", values_from = "log10P") %>%
  as.data.frame() %>%
  column_to_rownames(var = "kinase") %>%
  as.matrix() %>%
  t()

phospho_mat <- rppa_pho %>%
  semi_join(ka_protreg, by = c("protein" = "kinase", "sample")) %>%
  arrange(sample) %>%
  pivot_wider(names_from = "sample", values_from = "psite_expr") %>%
  unite(protein, psite, col = psite, sep = "-") %>%
  as.data.frame() %>%
  column_to_rownames(var = "psite") %>%
  as.matrix() %>%
  t()

kinases <- colnames(ka_mat)
psites <- colnames(phospho_mat)

KA_protreg_rppa_pho_cor_all <- rcorr(phospho_mat, ka_mat) %>%
  joinMatrices(kinases, psites) %>%
  select(-n, pearson_r=r, p_value=P) %>%
  mutate(class1 = "All pairs", class2="Substrate-based", class3 = "Protein reg-out") %>%
  filter(kinaseA != kinaseB)

# substrate-based protein non reg-out
proteins <- intersect(unique(ka_nonreg$kinase), unique(rppa_pho$protein))
samples <- intersect(unique(ka_nonreg$sample), unique(rppa_pho$sample))

ka_mat <- ka_nonreg %>%
  #semi_join(rppa_pho, by = c("sample")) %>%
  semi_join(rppa_pho, by = c("kinase" = "protein", "sample")) %>%
  arrange(sample) %>%
  pivot_wider(names_from = "sample", values_from = "log10P") %>%
  as.data.frame() %>%
  column_to_rownames(var = "kinase") %>%
  as.matrix() %>%
  t()

phospho_mat <- rppa_pho %>%
  semi_join(ka_nonreg, by = c("protein" = "kinase", "sample")) %>%
  arrange(sample) %>%
  pivot_wider(names_from = "sample", values_from = "psite_expr") %>%
  unite(protein, psite, col = psite, sep = "-") %>%
  as.data.frame() %>%
  column_to_rownames(var = "psite") %>%
  as.matrix() %>%
  t()

kinases <- colnames(ka_mat)
psites <- colnames(phospho_mat)

KA_nonreg_rppa_pho_cor_all <- rcorr(phospho_mat, ka_mat) %>%
  joinMatrices(kinases, psites) %>%
  select(-n, pearson_r=r, p_value=P) %>%
  mutate(class1 = "All pairs", class2="Substrate-based", class3 = "Protein non reg-out") %>%
  filter(kinaseA != kinaseB)

# regulatory phosphosite-based protein reg-out
proteins <- intersect(unique(ka_kinPsite_protreg$kinase), unique(rppa_pho$protein))
samples <- intersect(unique(ka_kinPsite_protreg$sample), unique(rppa_pho$sample))

ka_mat <- ka_kinPsite_protreg %>%
  #semi_join(rppa_pho, by = c("sample")) %>%
  semi_join(rppa_pho, by = c("kinase" = "protein", "sample")) %>%
  arrange(sample) %>%
  pivot_wider(names_from = "sample", values_from = "log10P") %>%
  as.data.frame() %>%
  column_to_rownames(var = "kinase") %>%
  as.matrix() %>%
  t()

phospho_mat <- rppa_pho %>%
  semi_join(ka_kinPsite_protreg, by = c("protein" = "kinase", "sample")) %>%
  arrange(sample) %>%
  pivot_wider(names_from = "sample", values_from = "psite_expr") %>%
  unite(protein, psite, col = psite, sep = "-") %>%
  as.data.frame() %>%
  column_to_rownames(var = "psite") %>%
  as.matrix() %>%
  t()

kinases <- colnames(ka_mat)
psites <- colnames(phospho_mat)

KA_protreg_KinPsite_rppa_pho_cor_all <- rcorr(phospho_mat, ka_mat) %>%
  joinMatrices(kinases, psites) %>%
  filter(n >= 5) %>%
  select(-n, pearson_r=r, p_value=P) %>%
  mutate(class1 = "All pairs", class2="Psite-based", class3 = "Protein reg-out") %>%
  filter(kinaseA != kinaseB)

# regulatory phosphosite-based protein non reg-out
proteins <- intersect(unique(ka_kinPsite_nonreg$kinase), unique(rppa_pho$protein))
samples <- intersect(unique(ka_kinPsite_nonreg$sample), unique(rppa_pho$sample))

ka_mat <- ka_kinPsite_nonreg %>%
  #semi_join(rppa_pho, by = c("sample")) %>%
  semi_join(rppa_pho, by = c("kinase" = "protein", "sample")) %>%
  arrange(sample) %>%
  pivot_wider(names_from = "sample", values_from = "log10P") %>%
  as.data.frame() %>%
  column_to_rownames(var = "kinase") %>%
  as.matrix() %>%
  t()

phospho_mat <- rppa_pho %>%
  semi_join(ka_kinPsite_nonreg, by = c("protein" = "kinase", "sample")) %>%
  arrange(sample) %>%
  pivot_wider(names_from = "sample", values_from = "psite_expr") %>%
  unite(protein, psite, col = psite, sep = "-") %>%
  as.data.frame() %>%
  column_to_rownames(var = "psite") %>%
  as.matrix() %>%
  t()

kinases <- colnames(ka_mat)
psites <- colnames(phospho_mat)

KA_nonreg_KinPsite_rppa_pho_cor_all <- rcorr(phospho_mat, ka_mat) %>%
  joinMatrices(kinases, psites) %>%
  filter(n >= 5) %>%
  select(-n, pearson_r=r, p_value=P) %>%
  mutate(class1 = "All pairs", class2="Psite-based", class3 = "Protein non reg-out") %>%
  filter(kinaseA != kinaseB)


#' join all correlations
all_cor <- bind_rows(KA_protreg_rppa_pho_cor,
                     KA_protreg_KinPsite_rppa_pho_cor,
                     KA_nonreg_rppa_pho_cor,
                     KA_nonreg_KinPsite_rppa_pho_cor,
                     KA_protreg_rppa_pho_cor_all,
                     KA_protreg_KinPsite_rppa_pho_cor_all,
                     KA_nonreg_rppa_pho_cor_all,
                     KA_nonreg_KinPsite_rppa_pho_cor_all) %>%
  select(kinaseB, kinaseA, psite, everything())
write_tsv(x = all_cor, file = "./output/files/corr_kin_activity_rppa.txt")

#' plot the distributions
KA_rppa_pho_cor_plot <- all_cor %>%
  filter(class3 == "Protein reg-out") %>%
  mutate(class1 = fct_recode(class1, Others = "All pairs", `Same kinase` = "Kinase-RPPA")) %>%
  mutate(class2 = fct_recode(class2, `Phosphosite-based activities` = "Psite-based", `Substrate-based activities` = "Substrate-based")) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = class1, y = pearson_r), notch = F, outlier.shape = NA, color = "grey60") +
  geom_jitter(mapping = aes(x = class1, y = pearson_r, color = class1), alpha = 0.3, width = 0.1) +
  facet_wrap( ~ class2) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()) +
  stat_compare_means(method = "t.test", mapping = aes(x = class1, y = pearson_r), label.x = 2.4, label.y = -0.6, size = 3) +
  coord_flip() +
  scale_color_manual(values = c("#bdbdbd", "#67a9cf"), guide = F) +
  labs(x = "", y = "Pearson's r (kinase activity vs RPPA phosphosite)")

#+ fig.width=7, fig.height=2
KA_rppa_pho_cor_plot

ggsave(filename = "kinase_protreg_rppa_psite_correlation_box.png", plot = KA_rppa_pho_cor_plot, path = "./output/plots/kinase_rppa_correlation/", height = 2, width = 7)
ggsave(filename = "kinase_protreg_rppa_psite_correlation_box.pdf", plot = KA_rppa_pho_cor_plot, path = "./output/plots/kinase_rppa_correlation/", height = 2, width = 7)
unlink("kinase_protreg_rppa_psite_correlation_box.png")
unlink("kinase_protreg_rppa_psite_correlation_box.pdf")


KA_rppa_pho_cor_plot <- all_cor %>%
  filter(class3 == "Protein reg-out") %>%
  mutate(class1 = fct_recode(class1, Others = "All pairs", `Same kinase` = "Kinase-RPPA")) %>%
  mutate(class2 = fct_recode(class2, `Phosphosite-based activities` = "Psite-based", `Substrate-based activities` = "Substrate-based")) %>%
  ggplot() +
  geom_histogram(mapping = aes(x = pearson_r, y = stat(density), fill = class1), position = "identity", alpha = 0.5, bins = 12) +
  facet_wrap( ~ class2) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10)) +
  scale_fill_manual(values = c("#bdbdbd", "#67a9cf")) +
  labs(x = "Pearson's r (kinase activity vs RPPA phosphosite)", y = "Density", fill = "Class")

#+ fig.width=7, fig.height=2
KA_rppa_pho_cor_plot

ggsave(filename = "kinase_protreg_rppa_psite_correlation_hist.png", plot = KA_rppa_pho_cor_plot, path = "./output/plots/kinase_rppa_correlation/", height = 2, width = 7)
ggsave(filename = "kinase_protreg_rppa_psite_correlation_hist.pdf", plot = KA_rppa_pho_cor_plot, path = "./output/plots/kinase_rppa_correlation/", height = 2, width = 7)
unlink("kinase_protreg_rppa_psite_correlation_hist.png")
unlink("kinase_protreg_rppa_psite_correlation_hist.pdf")


KA_rppa_pho_cor_plot <- all_cor %>%
  filter(class3 == "Protein non reg-out") %>%
  mutate(class1 = fct_recode(class1, Others = "All pairs", `Same kinase` = "Kinase-RPPA")) %>%
  mutate(class2 = fct_recode(class2, `Phosphosite-based activities` = "Psite-based", `Substrate-based activities` = "Substrate-based")) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = class1, y = pearson_r), notch = F, outlier.shape = NA, color = "grey60") +
  geom_jitter(mapping = aes(x = class1, y = pearson_r, color = class1), alpha = 0.3, width = 0.1) +
  facet_wrap( ~ class2) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()) +
  stat_compare_means(method = "t.test", mapping = aes(x = class1, y = pearson_r), label.x = 2.4, label.y = -0.6, size = 3) +
  coord_flip() +
  scale_color_manual(values = c("#bdbdbd", "#67a9cf"), guide = F) +
  labs(x = "", y = "Pearson's r (kinase activity vs RPPA phosphosite)")

#+ fig.width=7, fig.height=2
KA_rppa_pho_cor_plot

ggsave(filename = "kinase_nonreg_rppa_psite_correlation_box.png", plot = KA_rppa_pho_cor_plot, path = "./output/plots/kinase_rppa_correlation/", height = 2, width = 7)
ggsave(filename = "kinase_nonreg_rppa_psite_correlation_box.pdf", plot = KA_rppa_pho_cor_plot, path = "./output/plots/kinase_rppa_correlation/", height = 2, width = 7)
unlink("kinase_nonreg_rppa_psite_correlation_box.png")
unlink("kinase_nonreg_rppa_psite_correlation_box.pdf")


KA_rppa_pho_cor_plot <- all_cor %>%
  filter(class3 == "Protein non reg-out") %>%
  mutate(class1 = fct_recode(class1, Others = "All pairs", `Auto-psite` = "Kinase-RPPA")) %>%
  mutate(class2 = fct_recode(class2, `Phosphosite` = "Psite-based", `Substrate` = "Substrate-based")) %>%
  mutate(class2 = fct_relevel(class2, "Substrate")) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = class2, y = pearson_r, fill = class1), outlier.shape = NA, color = "black") +
  geom_beeswarm(mapping = aes(x = class2, y = pearson_r, alpha = class1), dodge.width = 0.75, show.legend = F) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(colour = "black", size = 20),
    legend.text = element_text(colour = "black", size = 18),
    axis.title = element_text(colour = "black", size = 20),
    axis.text = element_text(colour = "black", size = 18)) +
  stat_compare_means(mapping = aes(x = class2, y = pearson_r, group = class1, label = paste0("p = ", ..p.format..)), method = "wilcox.test", label.y = 0.7, size = 6) +
  scale_alpha_manual(values = c(0.3, 1)) +
  scale_fill_manual(values = c("#bdbdbd", "#67a9cf")) +
  guides(fill=guide_legend(ncol=2)) +
  labs(x = "Activity type", y = "Pearson's r (kinase activity vs RPPA psite)", fill = "Kinase-psite")

#+ fig.width=4, fig.height=6
KA_rppa_pho_cor_plot

ggsave(filename = "kinase_nonreg_rppa_psite_correlation_box2.png", plot = KA_rppa_pho_cor_plot, path = "./output/plots/kinase_rppa_correlation/", height = 6, width = 5)
ggsave(filename = "kinase_nonreg_rppa_psite_correlation_box2.pdf", plot = KA_rppa_pho_cor_plot, path = "./output/plots/kinase_rppa_correlation/", height = 6, width = 5)


KA_rppa_pho_cor_plot <- all_cor %>%
  filter(class3 == "Protein non reg-out") %>%
  mutate(class1 = fct_recode(class1, Others = "All pairs", `Same kinase` = "Kinase-RPPA")) %>%
  mutate(class2 = fct_recode(class2, `Phosphosite-based activities` = "Psite-based", `Substrate-based activities` = "Substrate-based")) %>%
  ggplot() +
  geom_histogram(mapping = aes(x = pearson_r, y = stat(density), fill = class1), position = "identity", alpha = 0.5, bins = 15) +
  facet_wrap( ~ class2) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10)) +
  scale_fill_manual(values = c("#bdbdbd", "#67a9cf")) +
  labs(x = "Pearson's r (kinase activity vs RPPA phosphosite)", y = "Density", fill = "Class")

#+ fig.width=7, fig.height=2
KA_rppa_pho_cor_plot

ggsave(filename = "kinase_nonreg_rppa_psite_correlation_hist.png", plot = KA_rppa_pho_cor_plot, path = "./output/plots/kinase_rppa_correlation/", height = 2, width = 7)
ggsave(filename = "kinase_nonreg_rppa_psite_correlation_hist.pdf", plot = KA_rppa_pho_cor_plot, path = "./output/plots/kinase_rppa_correlation/", height = 2, width = 7)
unlink("kinase_nonreg_rppa_psite_correlation_hist.png")
unlink("kinase_nonreg_rppa_psite_correlation_hist.pdf")

