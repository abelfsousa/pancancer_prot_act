# Exploratory analysis of the sources of kinase and TF regulation


library(tidyverse)
library(viridis)
library(UpSetR)
library(RColorBrewer)



# load kinase activity data
k_subN <- 3
kin_act <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-n, -source_type) %>%
  rename(kin_activity = log10P) %>%
  select(kinase, everything())

# kin_act_imp <- data.table::fread(file = "./data/Danish/kinaseActMatImputed.tsv") %>%
#   as_tibble() %>%
#   rename(kinase=V1) %>%
#   pivot_longer(-kinase, names_to = "sample", values_to = "kin_activity")

kin_act_not_kin_sites <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf_notKinSites.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-n, -source_type) %>%
  rename(kin_activity = log10P) %>%
  select(kinase, everything())

kin_act_prot_not_reg <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-n, -source_type) %>%
  rename(kin_activity = log10P) %>%
  select(kinase, everything())


# load TF activity data
tf_act <- data.table::fread("./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  as_tibble() %>%
  rename(tf=V1) %>%
  pivot_longer(-tf, names_to = "sample", values_to = "tf_activity")


#' load CNV data
cnv <- read_tsv(file = "./output/files/cnv.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "cnv")


# load gene expression data (log2fc)
rna <- read_tsv(file = "./output/files/transcriptomics_log2fc.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "rna_log2fc") %>%
  filter(!is.na(rna_log2fc))


# load protein abundance data (log2fc)
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "prot_log2fc") %>%
  filter(!is.na(prot_log2fc))


# load phosphorylation data (log2fc)
# protein regressed-out quantile normalized phosphorylation data
phospho <- read_tsv(file = "./output/files/phosphoproteomicsQ_Protreg_allSamp_withZtransf.txt.gz") %>%
  pivot_longer(-c(gene, psite, psites), names_to = "sample", values_to = "phos_log2fc") %>%
  filter(!is.na(phos_log2fc)) %>%
  select(gene, psite, psites, sample, phos_log2fc) %>%
  filter(psites == 1) %>%
  select(-gene, -psites) %>%
  separate(col = "psite", into = c("gene", "psite"), sep = "_") %>%
  mutate(residue = str_extract(psite, "^[STY]{1}")) %>%
  mutate(position = as.numeric(str_extract(psite, "[0-9]+"))) %>%
  select(-psite) %>%
  select(gene, position, residue, everything())

# quantile normalized phosphorylation data
phospho_prot_not_reg <- read_tsv(file = "./output/files/phosphoproteomicsQ.txt.gz") %>%
  pivot_longer(-c(gene, psite, psites), names_to = "sample", values_to = "phos_log2fc") %>%
  filter(!is.na(phos_log2fc)) %>%
  select(gene, psite, psites, sample, phos_log2fc) %>%
  filter(psites == 1) %>%
  select(-gene, -psites) %>%
  separate(col = "psite", into = c("gene", "psite"), sep = "_") %>%
  mutate(residue = str_extract(psite, "^[STY]{1}")) %>%
  mutate(position = as.numeric(str_extract(psite, "[0-9]+"))) %>%
  select(-psite) %>%
  select(gene, position, residue, everything())


# load kinase phosphosites
kin_psites <- read_tsv(file = "./output/files/cptac_kin_psites_funscoR.txt.gz")


# load TF phosphosites
tf_psites <- read_tsv(file = "./output/files/cptac_tf_psites_funscoR.txt.gz")


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


# load PhosphoSitePlus regulatory sites annotation
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


# set up a correlation function
corr <- function(df, x, y, method){
  a = df[[x]]
  b = df[[y]]
  cor <- broom::tidy(cor.test(a, b, method = method))
  
  res <- cor %>% select(estimate, p.value)
  
  res
}



# -- kinases

# quantify kinase variation

# use 1.75 as activity cutoff
# 1.75 corresponds to the ~96.7th percentil
#quantile(kin_act %>% pull(kin_activity), c(0.967))

kin_act_var <- kin_act %>%
  group_by(kinase) %>%
  summarise(n = n(), act_sd = sd(kin_activity), reg_samples = sum(abs(kin_activity) > 1.75)) %>%
  ungroup() %>%
  mutate(var_act_cutoff = as.numeric(reg_samples > n*0.05), var_act_sd = as.numeric(act_sd > quantile(act_sd, 0.5))) %>%
  rename(protein = kinase) %>%
  mutate(protein_type = "Kinases") %>%
  select(protein_type, everything())

kin_var_comparison_plot <- kin_act_var %>%
  pivot_longer(-c(protein_type, protein, n, act_sd, reg_samples), names_to = "variation_criteria", values_to = "variable") %>%
  ggplot(mapping = aes(x = variation_criteria, fill = as.character(variable))) +
  geom_bar(position = "dodge") +
  scale_fill_discrete(labels = c("0" = "Non variable", "1" = "Variable"), name = "Variation") +
  scale_x_discrete(labels = c("var_act_cutoff" = "sum(act > |1.75|) > 5%", "var_act_sd" = "SD > median(SD)"), name = "Variation criterion")

ggsave(filename = "kin_variation_criterion.png", plot = kin_var_comparison_plot, path = "./output/plots/kin_tf_regulation/", height = 4, width = 5)
ggsave(filename = "kin_variation_criterion.pdf", plot = kin_var_comparison_plot, path = "./output/plots/kin_tf_regulation/", height = 4, width = 5)


# correlate kinase activities with phosphosites mapping to the same kinase
kin_phos_cor <- kin_act %>%
  #inner_join(kin_psites[, -c(4:5)], by = c("kinase" = "gene")) %>%
  #inner_join(phospho, by = c("kinase" = "gene", "sample", "position", "residue")) %>%
  inner_join(phospho, by = c("kinase" = "gene", "sample")) %>%
  select(kinase, position, residue, everything()) %>%
  group_by(kinase, position, residue) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(n = map_dbl(data, nrow)) %>%
  mutate(corr = map(.x = data, .f = corr, x="kin_activity", y = "phos_log2fc", method = "pearson")) %>%
  select(-data) %>%
  unnest(cols = corr) %>%
  rename(pearson_r=estimate)


# correlate kinase activities with TF activities
kin_tf_cor <- kin_act %>%
  inner_join(tf_act, by = c("sample")) %>%
  select(kinase, tf, everything()) %>%
  group_by(kinase, tf) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(n = map_dbl(data, nrow)) %>%
  mutate(corr = map(.x = data, .f = corr, x="kin_activity", y = "tf_activity", method = "pearson")) %>%
  select(-data) %>%
  unnest(cols = corr) %>%
  rename(pearson_r=estimate)


# correlate kinase activities with respective protein abundance
kin_prot_cor <- kin_act %>%
  inner_join(protein, by = c("kinase" = "gene", "sample")) %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(n = map_dbl(data, nrow)) %>%
  mutate(corr = map(.x = data, .f = corr, x="kin_activity", y = "prot_log2fc", method = "pearson")) %>%
  select(-data) %>%
  unnest(cols = corr) %>%
  rename(pearson_r=estimate)


# correlate kinase activities with respective gene expression
kin_rna_cor <- kin_act %>%
  inner_join(rna, by = c("kinase" = "gene", "sample")) %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(n = map_dbl(data, nrow)) %>%
  mutate(corr = map(.x = data, .f = corr, x="kin_activity", y = "rna_log2fc", method = "pearson")) %>%
  select(-data) %>%
  unnest(cols = corr) %>%
  rename(pearson_r=estimate)


# correlate kinase activities with respective CNV scores
kin_cnv_cor <- kin_act %>%
  inner_join(cnv, by = c("kinase" = "gene", "sample")) %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(n = map_dbl(data, nrow)) %>%
  mutate(corr = map(.x = data, .f = corr, x="kin_activity", y = "cnv", method = "pearson")) %>%
  select(-data) %>%
  unnest(cols = corr) %>%
  rename(pearson_r=estimate)

# merge kinase correlations
kinase_correlations <- bind_rows(list(kin_phos_cor, kin_prot_cor, kin_rna_cor, kin_cnv_cor), .id = "data") %>%
  mutate(data = map_chr(.x = data, .f = ~ if(.x == "1"){"phosphorylation"}else if(.x == "2"){"protein"}else if(.x == "3"){"RNA"}else{"CNV"})) %>%
  mutate(protein_type = "kinase") %>%
  select(protein_type, correlated_data = data, protein = kinase, psite_position = position, psite_residue = residue, n, pearson_r, pearson_r_pval = p.value)


# -- TFs

# quantify TF variation

# use 3.89 as activity cutoff
# 3.89 corresponds to the ~96.7th percentil
#quantile(tf_act %>% pull(tf_activity), c(0.967))

tf_act_var <- tf_act %>%
  group_by(tf) %>%
  summarise(n = n(), act_sd = sd(tf_activity), reg_samples = sum(abs(tf_activity) > 3.89)) %>%
  ungroup() %>%
  mutate(var_act_cutoff = as.numeric(reg_samples > n*0.05), var_act_sd = as.numeric(act_sd > quantile(act_sd, 0.5))) %>%
  rename(protein = tf) %>%
  mutate(protein_type = "TFs") %>%
  select(protein_type, everything())

tf_var_comparison_plot <- tf_act_var %>%
  pivot_longer(-c(protein_type, protein, n, act_sd, reg_samples), names_to = "variation_criteria", values_to = "variable") %>%
  ggplot(mapping = aes(x = variation_criteria, fill = as.character(variable))) +
  geom_bar(position = "dodge") +
  scale_fill_discrete(labels = c("0" = "Non variable", "1" = "Variable"), name = "Variation") +
  scale_x_discrete(labels = c("var_act_cutoff" = "sum(act > |3.89|) > 5%", "var_act_sd" = "SD > median(SD)"), name = "Variation criterion")

ggsave(filename = "tf_variation_criterion.png", plot = tf_var_comparison_plot, path = "./output/plots/kin_tf_regulation/", height = 4, width = 5)
ggsave(filename = "tf_variation_criterion.pdf", plot = tf_var_comparison_plot, path = "./output/plots/kin_tf_regulation/", height = 4, width = 5)


# correlate TF activities with phosphosites mapping to the same TF
tf_phos_cor <- tf_act %>%
  inner_join(phospho, by = c("tf" = "gene", "sample")) %>%
  select(tf, position, residue, everything()) %>%
  group_by(tf, position, residue) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(n = map_dbl(data, nrow)) %>%
  mutate(corr = map(.x = data, .f = corr, x="tf_activity", y = "phos_log2fc", method = "pearson")) %>%
  select(-data) %>%
  unnest(cols = corr) %>%
  rename(pearson_r=estimate)


# correlate TF activities with kinase activities
tf_kin_cor <- tf_act %>%
  inner_join(kin_act, by = c("sample")) %>%
  select(tf, kinase, everything()) %>%
  group_by(tf, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(n = map_dbl(data, nrow)) %>%
  mutate(corr = map(.x = data, .f = corr, x="tf_activity", y = "kin_activity", method = "pearson")) %>%
  select(-data) %>%
  unnest(cols = corr) %>%
  rename(pearson_r=estimate)


# correlate TF activities with respective protein abundance
tf_prot_cor <- tf_act %>%
  inner_join(protein, by = c("tf" = "gene", "sample")) %>%
  group_by(tf) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(n = map_dbl(data, nrow)) %>%
  mutate(corr = map(.x = data, .f = corr, x="tf_activity", y = "prot_log2fc", method = "pearson")) %>%
  select(-data) %>%
  unnest(cols = corr) %>%
  rename(pearson_r=estimate)


# correlate TF activities with respective gene expression
tf_rna_cor <- tf_act %>%
  inner_join(rna, by = c("tf" = "gene", "sample")) %>%
  group_by(tf) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(n = map_dbl(data, nrow)) %>%
  mutate(corr = map(.x = data, .f = corr, x="tf_activity", y = "rna_log2fc", method = "pearson")) %>%
  select(-data) %>%
  unnest(cols = corr) %>%
  rename(pearson_r=estimate)


# correlate TF activities with respective CNV scores
tf_cnv_cor <- tf_act %>%
  inner_join(cnv, by = c("tf" = "gene", "sample")) %>%
  group_by(tf) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(n = map_dbl(data, nrow)) %>%
  mutate(corr = map(.x = data, .f = corr, x="tf_activity", y = "cnv", method = "pearson")) %>%
  select(-data) %>%
  unnest(cols = corr) %>%
  rename(pearson_r=estimate)

# merge TF correlations
tf_correlations <- bind_rows(list(tf_phos_cor, tf_prot_cor, tf_rna_cor, tf_cnv_cor), .id = "data") %>%
  mutate(data = map_chr(.x = data, .f = ~ if(.x == "1"){"phosphorylation"}else if(.x == "2"){"protein"}else if(.x == "3"){"RNA"}else{"CNV"})) %>%
  mutate(protein_type = "TF") %>%
  select(protein_type, correlated_data = data, protein = tf, psite_position = position, psite_residue = residue, n, pearson_r, pearson_r_pval = p.value)


# export correlation data (significant correlations)
correlations <- kinase_correlations %>%
  bind_rows(tf_correlations) %>%
  filter(pearson_r_pval < 0.05) %>%
  mutate(correlated_data = fct_relevel(correlated_data, "phosphorylation", "protein", "RNA", "CNV")) %>%
  arrange(protein_type, correlated_data, pearson_r_pval)

write_tsv(correlations, "./output/files/prot_act_corr_data_types.txt")


# plot percentage of kinases and TFs classified as variable
kin_tf_act_var <- kin_act_var %>%
  bind_rows(tf_act_var) %>%
  select(protein_type, protein, var_act_cutoff) %>%
  group_by(protein_type, var_act_cutoff) %>%
  tally() %>%
  ungroup() %>%
  group_by(protein_type) %>%
  mutate(p = n/sum(n)) %>%
  ungroup()

kin_tf_act_var2 <- kin_tf_act_var %>%
  mutate(p2 = map2_dbl(.x=var_act_cutoff, .y=p, .f = ~ if(.x == 0){1}else{.y}))

variation_plot <- kin_tf_act_var %>%
  ggplot(mapping = aes(x = protein_type, y = p, fill = as.character(var_act_cutoff))) +
  geom_col() +
  geom_text(data = kin_tf_act_var2, mapping = aes(x = protein_type, y = p2-0.03, label = n), size = 5) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, color = "black"),
    legend.position = "bottom") +
  scale_fill_manual(name = "Variation", labels = c("0" = "Non-variable", "1" = "Variable"), values = c("0"="grey", "1"="#74add1")) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"), name = "Percentage") +
  scale_x_discrete(name = "Protein type")

ggsave(filename = "kin_tf_variation_plot.png", plot = variation_plot, path = "./output/plots/kin_tf_regulation/", height = 4, width = 3)
ggsave(filename = "kin_tf_variation_plot.pdf", plot = variation_plot, path = "./output/plots/kin_tf_regulation/", height = 4, width = 3)


# plot percentage of kinases and TFs regulated by phospho, protein, RNA and CNV

# kinases
kin_phos_reg <- kin_phos_cor %>%
  filter(kinase %in% kin_act_var[kin_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  #mutate(phospho = if_else(abs(pearson_r) > 0.4, T, F)) %>%
  #mutate(phospho = if_else(abs(pearson_r) > 0.4 & p.value < 0.05, T, F)) %>%
  #mutate(phospho = if_else(p.value < 0.05, T, F)) %>%
  #mutate(phospho = if_else(p.value < 0.05 & n >= quantile(n,0.75), T, F)) %>%
  #mutate(phospho = if_else(p.value < 0.05 & (n >= 50 & n <= 200), T, F)) %>%
  #mutate(phospho = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05){0}else if(.x < 0.05 & .y < 0.2){1}else if(.x < 0.05 & (.y >= 0.2 & .y < 0.4)){2}else if(.x < 0.05 & .y >= 0.4){3}else{stop(0)})) %>%
  mutate(phospho = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  group_by(kinase) %>%
  #summarise(phospho = any(phospho)) %>%
  summarise(phospho = max(phospho)) %>%
  ungroup()
  #mutate(phospho = as.numeric(phospho))

kin_tf_reg <- kin_tf_cor %>%
  filter(kinase %in% kin_act_var[kin_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  filter(tf %in% tf_act_var[tf_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  #mutate(tf = if_else(abs(pearson_r) > 0.4, T, F)) %>%
  #mutate(tf = if_else(abs(pearson_r) > 0.4 & p.value < 0.05, T, F)) %>%
  #mutate(tf = if_else(p.value < 0.05, T, F)) %>%
  #mutate(tf = if_else(p.value < 0.05 & n >= quantile(n,0.75), T, F)) %>%
  mutate(tf = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  group_by(kinase) %>%
  #summarise(tf = any(tf)) %>%
  summarise(tf = max(tf)) %>%
  ungroup()
  #mutate(tf = as.numeric(tf))

kin_prot_reg <- kin_prot_cor %>%
  filter(kinase %in% kin_act_var[kin_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  #mutate(protein = if_else(abs(pearson_r) > 0.4, 1, 0)) %>%
  #mutate(protein = if_else(abs(pearson_r) > 0.4 & p.value < 0.05, 1, 0)) %>%
  #mutate(protein = if_else(p.value < 0.05, 1, 0)) %>%
  #mutate(protein = if_else(p.value < 0.05 & n >= quantile(n,0.75), 1, 0)) %>%
  #mutate(protein = if_else(p.value < 0.05 & (n >= 300 & n <= 750), 1, 0)) %>%
  mutate(protein = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  select(-pearson_r, -p.value, -n)

kin_rna_reg <- kin_rna_cor %>%
  filter(kinase %in% kin_act_var[kin_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  #mutate(rna = if_else(abs(pearson_r) > 0.4, 1, 0)) %>%
  #mutate(rna = if_else(abs(pearson_r) > 0.4 & p.value < 0.05, 1, 0)) %>%
  #mutate(rna = if_else(p.value < 0.05, 1, 0)) %>%
  #mutate(rna = if_else(p.value < 0.05 & n >= quantile(n,0.75), 1, 0)) %>%
  mutate(rna = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  select(-pearson_r, -p.value, -n)

kin_cnv_reg <- kin_cnv_cor %>%
  filter(kinase %in% kin_act_var[kin_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  #mutate(cnv = if_else(abs(pearson_r) > 0.4, 1, 0)) %>%
  #mutate(cnv = if_else(abs(pearson_r) > 0.4 & p.value < 0.05, 1, 0)) %>%
  #mutate(cnv = if_else(p.value < 0.05, 1, 0)) %>%
  #mutate(cnv = if_else(p.value < 0.05 & n >= quantile(n,0.75), 1, 0)) %>%
  mutate(cnv = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  select(-pearson_r, -p.value, -n)

kinase_regulation <- reduce(
  .x = list(kin_phos_reg, kin_tf_reg, kin_prot_reg, kin_rna_reg, kin_cnv_reg),
  .f = ~ full_join(.x, .y, by = "kinase")) %>%
  pivot_longer(-kinase, names_to = "data", values_to = "regulated") %>%
  arrange(data) %>%
  mutate(regulated = replace_na(regulated, 4)) %>%
  mutate(prot_type = "kinase") %>%
  rename(protein = kinase) %>%
  select(prot_type, everything())

#TFs
tf_phos_reg <- tf_phos_cor %>%
  filter(tf %in% tf_act_var[tf_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  #mutate(phospho = if_else(abs(pearson_r) > 0.4, T, F)) %>%
  #mutate(phospho = if_else(abs(pearson_r) > 0.4 & p.value < 0.05, T, F)) %>%
  #mutate(phospho = if_else(p.value < 0.05, T, F)) %>%
  #mutate(phospho = if_else(p.value < 0.05 & n >= quantile(n,0.75), T, F)) %>%
  #mutate(phospho = if_else(p.value < 0.05 & (n >= 50 & n <= 200), T, F)) %>%
  mutate(phospho = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  group_by(tf) %>%
  #summarise(phospho = any(phospho)) %>%
  summarise(phospho = max(phospho)) %>%
  ungroup()
  #mutate(phospho = as.numeric(phospho))

tf_kin_reg <- tf_kin_cor %>%
  filter(tf %in% tf_act_var[tf_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  filter(kinase %in% kin_act_var[kin_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  #mutate(kinase = if_else(abs(pearson_r) > 0.4, T, F)) %>%
  #mutate(kinase = if_else(abs(pearson_r) > 0.4 & p.value < 0.05, T, F)) %>%
  #mutate(kinase = if_else(p.value < 0.05, T, F)) %>%
  #mutate(kinase = if_else(p.value < 0.05 & n >= quantile(n,0.75), T, F)) %>%
  mutate(kinase = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  group_by(tf) %>%
  #summarise(kinase = any(kinase)) %>%
  summarise(kinase = max(kinase)) %>%
  ungroup()
  #mutate(kinase = as.numeric(kinase))

tf_prot_reg <- tf_prot_cor %>%
  filter(tf %in% tf_act_var[tf_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  #mutate(protein = if_else(abs(pearson_r) > 0.4, 1, 0)) %>%
  #mutate(protein = if_else(abs(pearson_r) > 0.4 & p.value < 0.05, 1, 0)) %>%
  #mutate(protein = if_else(p.value < 0.05, 1, 0)) %>%
  #mutate(protein = if_else(p.value < 0.05 & n >= quantile(n,0.75), 1, 0)) %>%
  #mutate(protein = if_else(p.value < 0.05 & (n >= 300 & n <= 750), 1, 0)) %>%
  mutate(protein = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  select(-pearson_r, -p.value, -n)

tf_rna_reg <- tf_rna_cor %>%
  filter(tf %in% tf_act_var[tf_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  #mutate(rna = if_else(abs(pearson_r) > 0.4, 1, 0)) %>%
  #mutate(rna = if_else(abs(pearson_r) > 0.4 & p.value < 0.05, 1, 0)) %>%
  #mutate(rna = if_else(p.value < 0.05, 1, 0)) %>%
  #mutate(rna = if_else(p.value < 0.05 & n >= quantile(n,0.75), 1, 0)) %>%
  mutate(rna = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  select(-pearson_r, -p.value, -n)

tf_cnv_reg <- tf_cnv_cor %>%
  filter(tf %in% tf_act_var[tf_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  #mutate(cnv = if_else(abs(pearson_r) > 0.4, 1, 0)) %>%
  #mutate(cnv = if_else(abs(pearson_r) > 0.4 & p.value < 0.05, 1, 0)) %>%
  #mutate(cnv = if_else(p.value < 0.05, 1, 0)) %>%
  #mutate(cnv = if_else(p.value < 0.05 & n >= quantile(n,0.75), 1, 0)) %>%
  mutate(cnv = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  select(-pearson_r, -p.value, -n)

tf_regulation <- reduce(
  .x = list(tf_phos_reg, tf_kin_reg, tf_prot_reg, tf_rna_reg, tf_cnv_reg),
  .f = ~ full_join(.x, .y, by = "tf")) %>%
  pivot_longer(-tf, names_to = "data", values_to = "regulated") %>%
  arrange(data) %>%
  mutate(regulated = replace_na(regulated, 4)) %>%
  mutate(prot_type = "tf") %>%
  rename(protein = tf) %>%
  select(prot_type, everything())


regulation_plot <- tf_regulation %>% 
  bind_rows(kinase_regulation) %>%
  mutate(regulated = as.character(regulated)) %>%
  filter(!data %in% c("tf", "kinase")) %>%
  group_by(prot_type, data, regulated) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(prot_type, data) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  mutate(data = fct_relevel(data, "cnv", "rna", "protein", "phospho")) %>%
  ggplot(mapping = aes(x = data, y = p, fill = regulated)) +
  geom_col(position = "stack", color = "black") + 
  facet_wrap(~ prot_type, scales = "free_x", labeller = labeller(prot_type = c("kinase" = "Kinases", "tf" = "TFs"))) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 20, color = "black"),
    axis.text = element_text(size = 18, color = "black"),
    axis.title = element_text(size = 20, color = "black"),
    legend.text = element_text(size = 18, color = "black"),
    legend.title = element_blank(),
    legend.position = "bottom") +
  #scale_fill_viridis(discrete = T, name = "Regulated", labels = c("0" = "Non-correlated (P > 5%)", "1" = "Correlated (P < 5% and r < 0.2)", "2" = "Correlated (0.2 <= r < 0.4)", "3" = "Correlated (r >= 0.4)", "4" = "Unknown")) +
  scale_fill_brewer(type = "qual", palette = "Set3", name = "Regulated", labels = c("0" = "Non-correlated (P > 5%)", "1" = "Correlated (P < 5% and r < 0.2)", "2" = "Correlated (0.2 <= r < 0.4)", "3" = "Correlated (r >= 0.4)", "4" = "Unknown")) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"), name = "Percentage") +
  scale_x_discrete(name = "Data", labels = c("cnv" = "CNV", "phospho" = "Phospho", "protein" = "Protein", "rna" = "RNA", "kinase" = "Kinase")) +
  guides(fill = guide_legend(nrow = 2))

ggsave(filename = "kin_tf_regulation_perct.png", plot = regulation_plot, path = "./output/plots/kin_tf_regulation/", height = 6, width = 9)
ggsave(filename = "kin_tf_regulation_perct.pdf", plot = regulation_plot, path = "./output/plots/kin_tf_regulation/", height = 6, width = 9)


# plot the intersection between the sets of kinases/TFs regulated by phospho, protein, RNA and CNV

kin_listInput <- list(
  TF = kin_tf_reg %>% filter(tf == 1) %>% pull(kinase),
  Phosphorylation = kin_phos_reg %>% filter(phospho == 1) %>% pull(kinase),
  Protein = kin_prot_reg %>% filter(protein == 1) %>% pull(kinase),
  RNA = kin_rna_reg %>% filter(rna == 1) %>% pull(kinase),
  CNV = kin_cnv_reg %>% filter(cnv == 1) %>% pull(kinase))

kin_regulation_sets <- upset(fromList(kin_listInput), order.by = "freq")

ggsave(filename = "kin_regulation_sets.png", plot = print(kin_regulation_sets), path = "./output/plots/kin_tf_regulation/", height = 4, width = 4)
ggsave(filename = "kin_regulation_sets.pdf", plot = print(kin_regulation_sets), path = "./output/plots/kin_tf_regulation/", height = 4, width = 4)


tf_listInput <- list(
  Kinase = tf_kin_reg %>% filter(kinase == 1) %>% pull(tf),
  Phosphorylation = tf_phos_reg %>% filter(phospho == 1) %>% pull(tf),
  Protein = tf_prot_reg %>% filter(protein == 1) %>% pull(tf),
  RNA = tf_rna_reg %>% filter(rna == 1) %>% pull(tf),
  CNV = tf_cnv_reg %>% filter(cnv == 1) %>% pull(tf))

tf_regulation_sets <- upset(fromList(tf_listInput), order.by = "freq")

ggsave(filename = "tf_regulation_sets.png", plot = print(tf_regulation_sets), path = "./output/plots/kin_tf_regulation/", height = 4, width = 4)
ggsave(filename = "tf_regulation_sets.pdf", plot = print(tf_regulation_sets), path = "./output/plots/kin_tf_regulation/", height = 4, width = 4)

