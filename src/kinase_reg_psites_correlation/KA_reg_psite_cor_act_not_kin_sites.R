#' ---
#' title: "Correlation of kinase activity to CPTAC regulatory phosphosites"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=6, fig.height=6)


#' load R packages
library(tidyverse)
library(ggpubr)
library(viridis)
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

ka_protreg <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf_notKinSites.txt.gz") %>%
#ka_protreg <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-n, -source_type)

ka_kinPsite_protreg <- read_tsv(file = "./output/files/CPTAC_KA_kin_psites_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(n >= k_psit) %>%
  select(-n) %>%
  rename(kinase=gene)


#' load CPTAC phosphorylation data
phospho <- read_tsv(file = "./output/files/phosphoproteomicsQ_Protreg_allSamp_withZtransf.txt.gz") %>%
  pivot_longer(-c(gene, psite, psites), names_to = "sample", values_to = "log2fc") %>%
  filter(!is.na(log2fc)) %>%
  filter(psites == 1) %>%
  select(-gene, -psites) %>%
  separate(col = "psite", into = c("gene", "psite"), sep = "_") %>%
  mutate(residue = str_extract(psite, "^[STY]{1}")) %>%
  mutate(position = as.numeric(str_extract(psite, "[0-9]+"))) %>%
  select(-psite) %>%
  select(gene, position, residue, everything())


#' load CPTAC kinase phosphosites with respective PhosphoSitePlus regulatory status and functional score
regPsite <- read_tsv(file = "./output/files/cptac_kin_psites_funscoR.txt.gz")


#' load PhosphoSitePlus regulatory sites annotation

# function to evaluate the phosphosites regulatory status
status <- function(strg){
  if (is.na(strg)){
    "unknown"
  }
  else if(str_detect(string = strg, pattern = "inhibited") & str_detect(string = strg, pattern = "induced")){
    "both"
  } else if(str_detect(string = strg, pattern = "inhibited")){
    "inhibition"
  } else if(str_detect(string = strg, pattern = "induced")){
    "activation"
  } else{
    "unknown"
  }
}

psp_reg <- read_tsv(file = "./data/PhosphoSitePlus/Regulatory_sites.gz", skip = 3) %>%
  filter(ORGANISM == "human") %>%
  select(gene=GENE, psite=MOD_RSD, reg_fun = ON_FUNCTION) %>%
  na.exclude() %>%
  separate(col = psite, into = c("psite", "modification")) %>%
  filter(modification == "p") %>%
  select(-modification) %>%
  mutate(residue = str_extract(psite, "^[A-Z]{1}"), position = str_extract(psite, "[0-9]+")) %>%
  select(-psite) %>%
  select(gene, residue, position, reg_fun) %>%
  distinct() %>%
  mutate(position = as.numeric(position)) %>%
  mutate(reg_status = map_chr(.x = reg_fun, .f = status)) %>%
  select(-reg_fun)


#' correlate kinase activities with CPTAC kinase phosphosites on the same kinase

# set a correlation function
corr <- function(df, x, y, method){
  a = df[[x]]
  b = df[[y]]
  cor <- broom::tidy(cor.test(a, b, method = method))
  
  res <- cor %>% select(estimate, p.value)
  
  res
}


KA_RegPsite_cor <- ka_protreg %>%
  inner_join(regPsite, by = c("kinase" = "gene")) %>%
  select(sample, kinase, activity=log10P, position, residue, funscoR=probabilities, PSP_reg_status) %>%
  inner_join(phospho, by = c("sample", "kinase" = "gene", "position", "residue")) %>%
  select(sample, kinase, activity, position, residue, log2fc, everything()) %>%
  group_by(kinase, position, residue, PSP_reg_status, funscoR) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(cor = map(.x = data, .f = corr, x="activity", y = "log2fc", method = "pearson")) %>%
  select(-data) %>%
  unnest() %>%
  rename(pearson_r=estimate) %>%
  mutate(class1 = "Kinase-Phosphosite", class2="Substrate-based") %>%
  rename(kinaseA=kinase) %>%
  mutate(kinaseB=kinaseA)

KA_KinPsite_RegPsite_cor <- ka_kinPsite_protreg %>%
  inner_join(regPsite, by = c("kinase" = "gene")) %>%
  select(sample, kinase, activity=log10P, position, residue, funscoR=probabilities, PSP_reg_status) %>%
  inner_join(phospho, by = c("sample", "kinase" = "gene", "position", "residue")) %>%
  select(sample, kinase, activity, position, residue, log2fc, everything()) %>%
  group_by(kinase, position, residue, PSP_reg_status, funscoR) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(cor = map(.x = data, .f = corr, x="activity", y = "log2fc", method = "pearson")) %>%
  select(-data) %>%
  unnest() %>%
  rename(pearson_r=estimate) %>%
  mutate(class1 = "Kinase-Phosphosite", class2="Psite-based") %>%
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
    separate(col = "kinaseA", into = c("kinaseA", "position", "residue"), sep = "_") %>%
    mutate(position = as.numeric(position))
  
  df
}

# set up a function to join the rcorr output into a single tibble
joinMatrices <- function(x, kin, pho){
  
  rcorr_matrices <- x
  rcorr_names <- names(rcorr_matrices)
  
  KIN <- kin
  PHO <- pho
  
  dfs <- map2(.x = rcorr_matrices, .y = rcorr_names, .f = matrix_to_tb, k = KIN, p = PHO)
  
  tb <- reduce(.x = dfs, .f = ~ inner_join(.x, .y, by = c("kinaseA", "position", "residue", "kinaseB")))
  
  tb
}


# substrate-based
ka_mat <- ka_protreg %>%
  semi_join(regPsite, by = c("kinase" = "gene")) %>%
  arrange(sample) %>%
  pivot_wider(names_from = "sample", values_from = "log10P") %>%
  as.data.frame() %>%
  column_to_rownames(var = "kinase") %>%
  as.matrix() %>%
  t()

kinases <- colnames(ka_mat)

phospho_mat <- phospho %>%
  semi_join(regPsite, by = c("gene", "position", "residue")) %>%
  arrange(sample) %>%
  pivot_wider(names_from = "sample", values_from = "log2fc") %>%
  unite(gene, position, residue, col = psite) %>%
  as.data.frame() %>%
  column_to_rownames(var = "psite") %>%
  as.matrix() %>%
  t()

psites <- colnames(phospho_mat)

KA_RegPsite_cor_all <- rcorr(phospho_mat, ka_mat) %>%
  joinMatrices(kinases, psites) %>%
  filter(n > 10) %>%
  select(-n, pearson_r=r, p.value=P) %>%
  inner_join(regPsite, by = c("kinaseA"="gene", "position", "residue")) %>%
  rename(funscoR=probabilities) %>%
  mutate(class1 = "All pairs", class2="Substrate-based") %>%
  filter(kinaseA != kinaseB)


# regulatory phosphosite-based
ka_mat <- ka_kinPsite_protreg %>%
  semi_join(regPsite, by = c("kinase" = "gene")) %>%
  arrange(sample) %>%
  pivot_wider(names_from = "sample", values_from = "log10P") %>%
  as.data.frame() %>%
  column_to_rownames(var = "kinase") %>%
  as.matrix() %>%
  t()

kinases <- colnames(ka_mat)

phospho_mat <- phospho %>%
  semi_join(regPsite, by = c("gene", "position", "residue")) %>%
  filter(sample %in% rownames(ka_mat)) %>%
  arrange(sample) %>%
  pivot_wider(names_from = "sample", values_from = "log2fc") %>%
  unite(gene, position, residue, col = psite) %>%
  as.data.frame() %>%
  column_to_rownames(var = "psite") %>%
  as.matrix() %>%
  t()

psites <- colnames(phospho_mat)

KA_KinPsite_RegPsite_cor_all <- rcorr(phospho_mat, ka_mat) %>%
  joinMatrices(kinases, psites) %>%
  filter(n > 10) %>%
  select(-n, pearson_r=r, p.value=P) %>%
  inner_join(regPsite, by = c("kinaseA"="gene", "position", "residue")) %>%
  rename(funscoR=probabilities) %>%
  mutate(class1 = "All pairs", class2="Psite-based") %>%
  filter(kinaseA != kinaseB)


#' join all correlations
all_cor <- bind_rows(KA_RegPsite_cor, KA_KinPsite_RegPsite_cor, KA_RegPsite_cor_all, KA_KinPsite_RegPsite_cor_all) %>%
  select(kinaseB, kinaseA, position, residue, everything())


#' plot the distributions
KA_RegPsite_cor_distribution <- all_cor %>%
  filter(class1 == "Kinase-Phosphosite", class2 == "Substrate-based") %>%
  mutate(PSP_reg_status = fct_recode(PSP_reg_status, PsP = "known", Others = "unknown")) %>%
  #mutate(PSP_reg_status = pmap_chr(.l = ., .f = ~ if(..5 == "unknown"){"unknown"}else if(..5 == "known" & ..6 < 0.8){"unknown"}else{"known"})) %>%
  ggplot(mapping = aes(x = PSP_reg_status, y = pearson_r)) +
  geom_boxplot(notch = T, outlier.shape = NA, color = "grey60") +
  geom_point(mapping = aes(color = PSP_reg_status), position = position_jitterdodge(jitter.width = 0.3), alpha = 0.1) +
  coord_flip() +
  scale_color_manual(values = c("#018571", "#a6611a"), guide = F) +
  stat_compare_means(method = "t.test", label.x = 2.4, label.y = 0.5, size = 3.5) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 10),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()) +
  labs(x = "", y = "Pearson's r (kinase activity vs kinase phosphosite)")

#+ fig.width=7, fig.height=2
KA_RegPsite_cor_distribution

ggsave(filename = "kinase_phosphosite_correlation_psp_status_box.png", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 2, width = 7)
ggsave(filename = "kinase_phosphosite_correlation_psp_status_box.pdf", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 2, width = 7)
unlink("kinase_phosphosite_correlation_psp_status_box.png")
unlink("kinase_phosphosite_correlation_psp_status_box.pdf")

KA_RegPsite_cor_distribution <- all_cor %>%
  filter(class1 == "Kinase-Phosphosite", class2 == "Substrate-based") %>%
  left_join(psp_reg, by = c("kinaseA" = "gene", "residue", "position")) %>%
  mutate(reg_status = replace_na(reg_status, "unknown")) %>%
  filter(pmap_lgl(.l = ., .f = ~ if(..5 == "unknown"){TRUE}else{if(..5 == "known" & ..11 == "activation"){TRUE}else{FALSE}})) %>%
  mutate(PSP_reg_status = fct_recode(PSP_reg_status, "Activating" = "known", "Unknown" = "unknown")) %>%
  mutate(PSP_reg_status = fct_relevel(PSP_reg_status, "Unknown")) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = PSP_reg_status, y = pearson_r), notch = T, outlier.shape = NA, color = "black", fill = "#f0f0f0") +
  geom_beeswarm(mapping = aes(x = PSP_reg_status, y = pearson_r, color = PSP_reg_status, alpha = PSP_reg_status), alpha = 0.4) +
  scale_color_manual(values = c("#a6611a", "#018571"), guide = F) +
  stat_compare_means(mapping = aes(x = PSP_reg_status, y = pearson_r), method = "wilcox.test", label.y = 1, size = 5, label.sep = ", ") +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 18),
    axis.text = element_text(colour = "black", size = 16)) +
  labs(x = "PsP annotation", y = "Pearson's r\n(kinase activity vs kinase phosphosite)")

#+ fig.width=4, fig.height=6
KA_RegPsite_cor_distribution

ggsave(filename = "kinase_phosphosite_correlation_psp_status_box2.png", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 6, width = 4)
ggsave(filename = "kinase_phosphosite_correlation_psp_status_box2.pdf", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 6, width = 4)


p_value <- all_cor %>%
  filter(class1 == "Kinase-Phosphosite", class2 == "Substrate-based") %>%
  mutate(PSP_reg_status = as.factor(PSP_reg_status)) %>%
  as.data.frame() %>%
  t.test(formula = pearson_r ~ PSP_reg_status, data = .) %>%
  pluck("p.value") %>%
  signif(2)

KA_RegPsite_cor_distribution <- all_cor %>%
  filter(class1 == "Kinase-Phosphosite", class2 == "Substrate-based") %>%
  mutate(PSP_reg_status = fct_recode(PSP_reg_status, Known = "known", Unknown = "unknown")) %>%
  mutate(PSP_reg_status = fct_relevel(PSP_reg_status, "Unknown", "Known")) %>%
  #mutate(PSP_reg_status = pmap_chr(.l = ., .f = ~ if(..5 == "unknown"){"unknown"}else if(..5 == "known" & ..6 < 0.8){"unknown"}else{"known"})) %>%
  ggplot(mapping = aes(x = pearson_r, y = stat(density), fill = PSP_reg_status)) +
  geom_histogram(position = "identity", alpha = 0.8, bins = 30) +
  annotate("text", x = 0.75, y = 1.8, label = str_c("T-test, p = ", p_value), size = 3.5) +
  scale_fill_manual(values = c("Known"="#018571", "Unknown"="#a6611a")) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "Pearson's r (kinase activity vs kinase phosphosite)", y = "Density", fill = "PsP annotation")

#+ fig.width=7, fig.height=2
KA_RegPsite_cor_distribution

ggsave(filename = "kinase_phosphosite_correlation_psp_status_hist.png", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 2, width = 7)
ggsave(filename = "kinase_phosphosite_correlation_psp_status_hist.pdf", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 2, width = 7)
unlink("kinase_phosphosite_correlation_psp_status_hist.png")
unlink("kinase_phosphosite_correlation_psp_status_hist.pdf")


KA_RegPsite_cor_distribution <- all_cor %>%
  filter(class2 == "Substrate-based") %>%
  mutate(class1 = fct_recode(class1, `Auto-phosphosite` = "Kinase-Phosphosite", Others = "All pairs")) %>%
  mutate(PSP_reg_status = fct_recode(PSP_reg_status, Yes = "known", No = "unknown")) %>%
  ggplot(mapping = aes(x = class1, y = pearson_r, fill = PSP_reg_status)) +
  geom_boxplot(notch = T, color = "grey60", alpha = 0.8) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.05) +
  coord_flip() +
  scale_fill_manual(values = c("#018571", "#a6611a")) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Pearson's r (kinase activity vs kinase phosphosite)", fill = "PSP annotation")

#+ fig.width=9, fig.height=3
KA_RegPsite_cor_distribution

ggsave(filename = "kinase_phosphosite_correlation_psp_status_auto_all_pairs.png", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 3, width = 9)
ggsave(filename = "kinase_phosphosite_correlation_psp_status_auto_all_pairs.pdf", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 3, width = 9)
unlink("kinase_phosphosite_correlation_psp_status_auto_all_pairs.png")
unlink("kinase_phosphosite_correlation_psp_status_auto_all_pairs.pdf")


KA_RegPsite_cor_distribution <- all_cor %>%
  filter(class2 == "Substrate-based" & PSP_reg_status == "known") %>%
  mutate(class1 = fct_recode(class1, `Auto-phosphosite` = "Kinase-Phosphosite", Others = "All pairs")) %>%
  ggplot(mapping = aes(x = class1, y = pearson_r, fill = class1)) +
  geom_boxplot(notch = T, color = "grey60", alpha = 0.8) +
  #geom_jitter(width = 0.1, alpha = 0.5, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("#bdbdbd", "#67a9cf"), guide=F) +
  stat_compare_means(method = "t.test", label.x = 2.4, label.y = 0.6, size = 3.5) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 10),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()) +
  labs(x = "", y = "Pearson's r (kinase activity vs kinase regulatory phosphosites)")

#+ fig.width=7, fig.height=2
KA_RegPsite_cor_distribution

ggsave(filename = "kinase_phosphosite_correlation_known_psp_status_auto_all_pairs_box.png", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 2, width = 7)
ggsave(filename = "kinase_phosphosite_correlation_known_psp_status_auto_all_pairs_box.pdf", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 2, width = 7)
unlink("kinase_phosphosite_correlation_known_psp_status_auto_all_pairs_box.png")
unlink("kinase_phosphosite_correlation_known_psp_status_auto_all_pairs_box.pdf")

KA_RegPsite_cor_distribution <- all_cor %>%
  filter(class2 == "Substrate-based" & PSP_reg_status == "known") %>%
  mutate(class1 = fct_recode(class1, `Auto-phosphosite` = "Kinase-Phosphosite", Others = "All pairs")) %>%
  ggplot(mapping = aes(x = pearson_r, y = stat(density), fill = class1)) +
  geom_histogram(position = "identity", alpha = 0.8) +
  scale_fill_manual(values = c("#bdbdbd", "#67a9cf")) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 10)) +
  #scale_x_continuous(limits = c(-1, 1)) +
  labs(x = "Pearson's r (kinase activity vs kinase regulatory phosphosites)", y = "Density", fill = "Class")

#+ fig.width=7, fig.height=2
KA_RegPsite_cor_distribution

ggsave(filename = "kinase_phosphosite_correlation_known_psp_status_auto_all_pairs_hist.png", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 2, width = 7)
ggsave(filename = "kinase_phosphosite_correlation_known_psp_status_auto_all_pairs_hist.pdf", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 2, width = 7)
unlink("kinase_phosphosite_correlation_known_psp_status_auto_all_pairs_hist.png")
unlink("kinase_phosphosite_correlation_known_psp_status_auto_all_pairs_hist.pdf")


KA_RegPsite_cor_distribution <- all_cor %>%
  filter(class2 == "Substrate-based" & PSP_reg_status == "known") %>%
  inner_join(psp_reg, by = c("kinaseA" = "gene", "residue", "position")) %>%
  filter(!(reg_status == "unknown" | reg_status == "both")) %>%
  mutate(class1 = fct_recode(class1, `Auto-phosphosite` = "Kinase-Phosphosite", Others = "All pairs")) %>%
  ggplot(mapping = aes(x = class1, y = pearson_r, fill = reg_status)) +
  geom_boxplot(notch = T, outlier.shape = NA) +
  #geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.4, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("#99d594", "#fc9272")) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Pearson's r (kinase activity vs kinase regulatory phosphosite)", fill = "PSP annotation")

#+ fig.width=9, fig.height=3
KA_RegPsite_cor_distribution

ggsave(filename = "kinase_phosphosite_correlation_known_psp_status_direction_auto_all_pairs_box.png", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 3, width = 9)
ggsave(filename = "kinase_phosphosite_correlation_known_psp_status_direction_auto_all_pairs_box.pdf", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 3, width = 9)
unlink("kinase_phosphosite_correlation_known_psp_status_direction_auto_all_pairs_box.png")
unlink("kinase_phosphosite_correlation_known_psp_status_direction_auto_all_pairs_box.pdf")


KA_RegPsite_cor_distribution <- all_cor %>%
  filter(class2 == "Substrate-based" & PSP_reg_status == "known") %>%
  inner_join(psp_reg, by = c("kinaseA" = "gene", "residue", "position")) %>%
  filter(reg_status == "activation") %>%
  mutate(class1 = fct_recode(class1, `Auto-phosphosite` = "Kinase-Phosphosite", Others = "All pairs")) %>%
  ggplot(mapping = aes(x = class1, y = pearson_r, fill = class1)) +
  geom_boxplot(notch = T, color = "grey60", alpha = 0.8) +
  #geom_jitter(width = 0.1, alpha = 0.5, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("#bdbdbd", "#67a9cf"), guide=F) +
  stat_compare_means(method = "t.test", label.x = 2.4, label.y = 0.6, size = 3.5) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 10),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()) +
  labs(x = "", y = "Pearson's r (kinase activity vs kinase regulatory phosphosites)")

#+ fig.width=7, fig.height=2
KA_RegPsite_cor_distribution

ggsave(filename = "kinase_phosphosite_correlation_known_psp_status_activating_auto_all_pairs_box.png", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 2, width = 7)
ggsave(filename = "kinase_phosphosite_correlation_known_psp_status_activating_auto_all_pairs_box.pdf", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 2, width = 7)
unlink("kinase_phosphosite_correlation_known_psp_status_activating_auto_all_pairs_box.png")
unlink("kinase_phosphosite_correlation_known_psp_status_activating_auto_all_pairs_box.pdf")

KA_RegPsite_cor_distribution <- all_cor %>%
  filter(class2 == "Substrate-based" & PSP_reg_status == "known") %>%
  inner_join(psp_reg, by = c("kinaseA" = "gene", "residue", "position")) %>%
  filter(reg_status == "activation") %>%
  mutate(class1 = fct_recode(class1, `Auto-psite` = "Kinase-Phosphosite", Others = "All pairs")) %>%
  mutate(class1 = fct_relevel(class1, "Others", "Auto-psite")) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = class1, y = pearson_r, fill = class1), notch = T, outlier.shape = NA, color = "black", alpha = 0.8) +
  geom_jitter(mapping = aes(x = class1, y = pearson_r, alpha = class1), width = 0.1, color = "black") +
  scale_fill_manual(values = c("#a6611a", "#018571")) +
  scale_alpha_manual(values = c(0.005,0.2), guide = F) +
  stat_compare_means(mapping = aes(x = class1, y = pearson_r), method = "t.test", label.y = 1.1, size = 5, label.sep = ",\n") +
  theme_classic() +
  theme(
    plot.margin = unit(c(0,0,1.5,0), "cm"),
    axis.title.x = element_text(colour = "black", size = 17, hjust = 1),
    axis.title.y = element_text(colour = "black", size = 17, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 16),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.2,0),
    legend.justification = c(0.4,1.7),
    legend.margin = margin(c(-0.2,0,0,0), unit="cm"),
    legend.text = element_text(colour = "black", size = 16)) +
  labs(x = "Kin-psite pair", y = "Pearson's r (KA vs kinase regulatory phosphosites)") +
  guides(fill = guide_legend(override.aes = list(alpha=1), nrow = 2))

#+ fig.width=2, fig.height=6
KA_RegPsite_cor_distribution

ggsave(filename = "kinase_phosphosite_correlation_known_psp_status_activating_auto_all_pairs_box2.png", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 6, width = 2)
ggsave(filename = "kinase_phosphosite_correlation_known_psp_status_activating_auto_all_pairs_box2.pdf", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 6, width = 2)
unlink("kinase_phosphosite_correlation_known_psp_status_activating_auto_all_pairs_box2.png")
unlink("kinase_phosphosite_correlation_known_psp_status_activating_auto_all_pairs_box2.pdf")


KA_RegPsite_cor_distribution <- all_cor %>%
  filter(class2 == "Substrate-based" & PSP_reg_status == "known") %>%
  inner_join(psp_reg, by = c("kinaseA" = "gene", "residue", "position")) %>%
  filter(reg_status == "activation") %>%
  mutate(class1 = fct_recode(class1, `Auto-phosphosite` = "Kinase-Phosphosite", Others = "All pairs")) %>%
  ggplot(mapping = aes(x = pearson_r, y = stat(density), fill = class1)) +
  geom_histogram(position = "identity", alpha = 0.8) +
  scale_fill_manual(values = c("#bdbdbd", "#67a9cf")) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 10)) +
  #scale_x_continuous(limits = c(-1, 1)) +
  labs(x = "Pearson's r (kinase activity vs kinase regulatory phosphosites)", y = "Density", fill = "Class")

#+ fig.width=7, fig.height=2
KA_RegPsite_cor_distribution

ggsave(filename = "kinase_phosphosite_correlation_known_psp_status_activating_auto_all_pairs_hist.png", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 2, width = 7)
ggsave(filename = "kinase_phosphosite_correlation_known_psp_status_activating_auto_all_pairs_hist.pdf", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 2, width = 7)
unlink("kinase_phosphosite_correlation_known_psp_status_activating_auto_all_pairs_hist.png")
unlink("kinase_phosphosite_correlation_known_psp_status_activating_auto_all_pairs_hist.pdf")

