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


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy

set.seed(123)


#' load kinase-activity inference data
k_subN <- 3
k_psit <- 3

ka_protreg <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf_notKinSites.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-n, -source_type)

ka_kinPsite_protreg <- read_tsv(file = "./output/files/CPTAC_KA_kin_psites_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(n >= k_psit) %>%
  select(-n) %>%
  rename(kinase=gene)


#' shuffle the kinases
ka_protreg_shuff <- ka_protreg %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(index = sample(seq_len(nrow(.)), nrow(.), replace = FALSE)) %>%
  mutate(kinase = kinase[index]) %>%
  select(-index) %>%
  unnest()

ka_kinPsite_protreg_shuff <- ka_kinPsite_protreg %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(index = sample(seq_len(nrow(.)), nrow(.), replace = FALSE)) %>%
  mutate(kinase = kinase[index]) %>%
  select(-index) %>%
  unnest()


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


#' shuffle phosphorylation data
phospho_shuff <- phospho %>%
  group_by(gene, position, residue) %>%
  nest() %>%
  ungroup() %>%
  mutate(index = sample(seq_len(nrow(.)), nrow(.), replace = FALSE)) %>%
  mutate(gene = gene[index], position = position[index], residue = residue[index]) %>%
  select(-index) %>%
  unnest()


#' load all and PSP regulatory phosphosites
regPsite <- read_tsv(file = "./output/files/cptac_kin_psites_funscoR.txt.gz")


#' correlate kinase activity with regulatory phosphosites

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
  mutate(class1 = "True pairs", class2="Substrate-based KA")


KA_RegPsite_shuff_cor <- ka_protreg_shuff %>%
  inner_join(regPsite, by = c("kinase" = "gene")) %>%
  select(sample, kinase, activity=log10P, position, residue, funscoR=probabilities, PSP_reg_status) %>%
  inner_join(phospho_shuff, by = c("sample", "kinase" = "gene", "position", "residue")) %>%
  select(sample, kinase, activity, position, residue, log2fc, everything()) %>%
  group_by(kinase, position, residue, PSP_reg_status, funscoR) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(cor = map(.x = data, .f = corr, x="activity", y = "log2fc", method = "pearson")) %>%
  select(-data) %>%
  unnest() %>%
  rename(pearson_r=estimate) %>%
  mutate(class1 = "Random pairs", class2="Substrate-based KA")


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
  mutate(class1 = "True pairs", class2="RegPsites-based KA")


KA_KinPsite_RegPsite_shuff_cor <- ka_kinPsite_protreg_shuff %>%
  inner_join(regPsite, by = c("kinase" = "gene")) %>%
  select(sample, kinase, activity=log10P, position, residue, funscoR=probabilities, PSP_reg_status) %>%
  inner_join(phospho_shuff, by = c("sample", "kinase" = "gene", "position", "residue")) %>%
  select(sample, kinase, activity, position, residue, log2fc, everything()) %>%
  group_by(kinase, position, residue, PSP_reg_status, funscoR) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(cor = map(.x = data, .f = corr, x="activity", y = "log2fc", method = "pearson")) %>%
  select(-data) %>%
  unnest() %>%
  rename(pearson_r=estimate) %>%
  mutate(class1 = "Random pairs", class2="RegPsites-based KA")


all_cor <- bind_rows(KA_RegPsite_cor, KA_RegPsite_shuff_cor, KA_KinPsite_RegPsite_cor, KA_KinPsite_RegPsite_shuff_cor)


#' plot the distributions
KA_RegPsite_cor_distribution <- all_cor %>%
  filter(class2 == "Substrate-based KA") %>%
  ggplot(mapping = aes(x = class1, y = pearson_r, fill = PSP_reg_status)) +
  geom_boxplot(notch = T, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.1, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("#018571", "#a6611a")) +
  theme(
    axis.title = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Pearson's r (kinase activity vs phosphosite)", fill = "PSP status")

#+ fig.width=9, fig.height=3
KA_RegPsite_cor_distribution

ggsave(filename = "KA_RegPsite_cor_distribution.png", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 3, width = 9)
unlink("KA_RegPsite_cor_distribution.png")


KA_RegPsite_cor_distribution <- all_cor %>%
  filter(class2 == "Substrate-based KA" & PSP_reg_status == "known") %>%
  ggplot(mapping = aes(x = class1, y = pearson_r, fill = class1)) +
  geom_boxplot(notch = T, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("#f7f7f7", "#67a9cf"), guide=F) +
  stat_compare_means(method = "wilcox.test", label.x = 2.3, label.y = -0.7, size = 3.5) +
  scale_y_continuous(limits = c(-0.7, 0.7)) +
  theme(
    axis.title = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Pearson's r (kinase activity vs PSP regulatory phosphosites)")

#+ fig.width=8, fig.height=2
KA_RegPsite_cor_distribution

ggsave(filename = "KA_RegPsite_cor_distribution_known_psp.png", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 2, width = 8)
unlink("KA_RegPsite_cor_distribution_known_psp.png")



#' load PhosphoSitePlus (PSP) regulatory sites
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
  mutate(position = as.numeric(position))



#' check the phosphosite regulation status (activation and/or inhibition)

# function to evaluate that
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

psp_reg <- psp_reg %>%
  mutate(reg_status = map_chr(.x = reg_fun, .f = status)) %>%
  select(-reg_fun)



#' plot the distribution
KA_RegPsite_cor_distribution <- all_cor %>%
  filter(class2 == "Substrate-based KA" & PSP_reg_status == "known") %>%
  inner_join(psp_reg, by = c("kinase" = "gene", "residue", "position")) %>%
  filter(!(reg_status == "unknown" | reg_status == "both")) %>%
  ggplot(mapping = aes(x = class1, y = pearson_r, fill = reg_status)) +
  geom_boxplot(notch = T, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.4, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("#99d594", "#fc8d59")) +
  theme(
    axis.title = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Pearson's r (kinase activity vs PSP regulatory phosphosite)", fill = "Regulatory status")

#+ fig.width=9, fig.height=3
KA_RegPsite_cor_distribution

ggsave(filename = "KA_RegPsite_cor_distribution_direction.png", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 3, width = 9)
unlink("KA_RegPsite_cor_distribution_direction.png")


KA_RegPsite_cor_distribution <- all_cor %>%
  filter(class2 == "Substrate-based KA" & PSP_reg_status == "known") %>%
  inner_join(psp_reg, by = c("kinase" = "gene", "residue", "position")) %>%
  filter(reg_status == "activation") %>%
  ggplot(mapping = aes(x = class1, y = pearson_r, fill = class1)) +
  geom_boxplot(notch = T, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("#f7f7f7", "#67a9cf"), guide=F) +
  stat_compare_means(method = "wilcox.test", label.x = 2.3, label.y = -0.7, size = 3.5) +
  scale_y_continuous(limits = c(-0.7, 0.7)) +
  theme(
    axis.title = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 10)) +
  labs(x = "", y = "Pearson's r (kinase activity vs PSP regulatory phosphosites)")

#+ fig.width=8, fig.height=2
KA_RegPsite_cor_distribution

ggsave(filename = "KA_RegPsite_cor_distribution_direction_activation.png", plot = KA_RegPsite_cor_distribution, path = "./output/plots/kinase_reg_psites_correlation/", height = 2, width = 8)
unlink("KA_RegPsite_cor_distribution_direction_activation.png")
