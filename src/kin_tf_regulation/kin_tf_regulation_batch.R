# Exploratory analysis of the sources of kinase and TF regulation


library(tidyverse)
library(viridis)
library(UpSetR)
library(RColorBrewer)



# load kinase activity data
k_subN <- 3
kin_act <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  #kin_act <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf_notKinSites.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-n, -source_type) %>%
  rename(kin_activity = log10P) %>%
  select(kinase, everything())

# kin_act <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ.txt.gz") %>%
#   filter(source_type == "DB_text-mining") %>%
#   filter(n >= k_subN) %>%
#   select(-n, -source_type) %>%
#   rename(kin_activity = log10P) %>%
#   select(kinase, everything())

kin_act_imp <- data.table::fread(file = "./data/Danish/kinaseActMatImputed.tsv") %>%
  as_tibble() %>%
  rename(kinase=V1) %>%
  pivot_longer(-kinase, names_to = "sample", values_to = "kin_activity")

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
phospho2 <- read_tsv(file = "./output/files/phosphoproteomicsQ.txt.gz") %>%
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


# load metadata
getSamples <- function(all_samples, data_types){
  
  l = length(data_types)
  
  samples <- all_samples %>%
    filter(data %in% data_types) %>%
    rename(info = data) %>%
    group_by(info, batch, cancer, tissue) %>%
    summarise(sample = list(sample)) %>%
    ungroup() %>%
    group_by(batch, cancer, tissue) %>%
    summarise(sample = list(sample)) %>%
    ungroup() %>%
    filter(map_dbl(.x=sample, .f=length) == l) %>%
    mutate(overlap = map(.x=sample, .f = ~ reduce(.x=.x, .f=intersect))) %>%
    select(-sample) %>%
    unnest(cols = overlap) %>%
    rename(sample = overlap)
  
  return(samples)
  
}

samples_metadata <- read_tsv("./output/files/all_samples_annotation.txt")
tf_samples <- getSamples(samples_metadata, data_types = c("mRNA", "protein")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue)

kin_samples <- getSamples(samples_metadata, data_types = c("phosphorylation")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue)


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

#ggsave(filename = "kin_variation_criterion.png", plot = kin_var_comparison_plot, path = "./output/plots/kin_tf_regulation/", height = 4, width = 5)
#ggsave(filename = "kin_variation_criterion.pdf", plot = kin_var_comparison_plot, path = "./output/plots/kin_tf_regulation/", height = 4, width = 5)


# correlate kinase activities with phosphosites mapping to the same kinase
kin_phos_cor <- kin_act %>%
  inner_join(phospho, by = c("kinase" = "gene", "sample")) %>%
  inner_join(kin_samples, by = "sample") %>%
  #filter(batch == "hcc-proteogenomics") %>%
  select(batch, kinase, position, residue, everything()) %>%
  group_by(batch, kinase, position, residue) %>%
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
  inner_join(kin_samples, by = "sample") %>%
  #filter(batch == "hcc-proteogenomics") %>%
  select(batch, kinase, tf, everything()) %>%
  group_by(batch, kinase, tf) %>%
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
  inner_join(kin_samples, by = "sample") %>%
  #filter(batch == "hcc-proteogenomics") %>%
  group_by(batch, kinase) %>%
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
  inner_join(kin_samples, by = "sample") %>%
  #filter(batch == "hcc-proteogenomics") %>%
  group_by(batch, kinase) %>%
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
  inner_join(kin_samples, by = "sample") %>%
  #filter(batch == "hcc-proteogenomics") %>%
  group_by(batch, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(n = map_dbl(data, nrow)) %>%
  mutate(corr = map(.x = data, .f = corr, x="kin_activity", y = "cnv", method = "pearson")) %>%
  select(-data) %>%
  unnest(cols = corr) %>%
  rename(pearson_r=estimate)


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

#ggsave(filename = "tf_variation_criterion.png", plot = tf_var_comparison_plot, path = "./output/plots/kin_tf_regulation/", height = 4, width = 5)
#ggsave(filename = "tf_variation_criterion.pdf", plot = tf_var_comparison_plot, path = "./output/plots/kin_tf_regulation/", height = 4, width = 5)


# correlate TF activities with phosphosites mapping to the same TF
tf_phos_cor <- tf_act %>%
  inner_join(phospho, by = c("tf" = "gene", "sample")) %>%
  inner_join(tf_samples, by = "sample") %>%
  #filter(batch == "hcc-proteogenomics") %>%
  select(batch, tf, position, residue, everything()) %>%
  group_by(batch, tf, position, residue) %>%
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
  inner_join(tf_samples, by = "sample") %>%
  #filter(batch == "hcc-proteogenomics") %>%
  select(batch, tf, kinase, everything()) %>%
  group_by(batch, tf, kinase) %>%
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
  inner_join(tf_samples, by = "sample") %>%
  #filter(batch == "hcc-proteogenomics") %>%
  group_by(batch, tf) %>%
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
  inner_join(tf_samples, by = "sample") %>%
  #filter(batch == "hcc-proteogenomics") %>%
  group_by(batch, tf) %>%
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
  inner_join(tf_samples, by = "sample") %>%
  #filter(batch == "hcc-proteogenomics") %>%
  group_by(batch, tf) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(n = map_dbl(data, nrow)) %>%
  mutate(corr = map(.x = data, .f = corr, x="tf_activity", y = "cnv", method = "pearson")) %>%
  select(-data) %>%
  unnest(cols = corr) %>%
  rename(pearson_r=estimate)


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
  geom_text(data = kin_tf_act_var2, mapping = aes(x = protein_type, y = p2-0.03, label = n)) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = "black"),
    legend.position = "bottom") +
  scale_fill_manual(name = "Variation", labels = c("0" = "Non-variable", "1" = "Variable"), values = c("0"="grey", "1"="#74add1")) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"), name = "Percentage") +
  scale_x_discrete(name = "Protein type")

#ggsave(filename = "kin_tf_variation_plot.png", plot = variation_plot, path = "./output/plots/kin_tf_regulation/", height = 4, width = 2)
#ggsave(filename = "kin_tf_variation_plot.pdf", plot = variation_plot, path = "./output/plots/kin_tf_regulation/", height = 4, width = 2)


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
  mutate(phospho = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05 | is.na(.x)){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  group_by(batch, kinase) %>%
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
  mutate(tf = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05 | is.na(.x)){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  group_by(batch, kinase) %>%
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
  mutate(protein = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05 | is.na(.x)){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  select(-pearson_r, -p.value, -n)

kin_rna_reg <- kin_rna_cor %>%
  filter(kinase %in% kin_act_var[kin_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  #mutate(rna = if_else(abs(pearson_r) > 0.4, 1, 0)) %>%
  #mutate(rna = if_else(abs(pearson_r) > 0.4 & p.value < 0.05, 1, 0)) %>%
  #mutate(rna = if_else(p.value < 0.05, 1, 0)) %>%
  #mutate(rna = if_else(p.value < 0.05 & n >= quantile(n,0.75), 1, 0)) %>%
  mutate(rna = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05 | is.na(.x)){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  select(-pearson_r, -p.value, -n)

kin_cnv_reg <- kin_cnv_cor %>%
  filter(kinase %in% kin_act_var[kin_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  #mutate(cnv = if_else(abs(pearson_r) > 0.4, 1, 0)) %>%
  #mutate(cnv = if_else(abs(pearson_r) > 0.4 & p.value < 0.05, 1, 0)) %>%
  #mutate(cnv = if_else(p.value < 0.05, 1, 0)) %>%
  #mutate(cnv = if_else(p.value < 0.05 & n >= quantile(n,0.75), 1, 0)) %>%
  mutate(cnv = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05 | is.na(.x)){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  select(-pearson_r, -p.value, -n)

kinase_regulation <- reduce(
  .x = list(kin_phos_reg, kin_tf_reg, kin_prot_reg, kin_rna_reg, kin_cnv_reg),
  .f = ~ full_join(.x, .y, by = c("batch", "kinase"))) %>%
  pivot_longer(-c(batch, kinase), names_to = "data", values_to = "regulated") %>%
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
  mutate(phospho = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05 | is.na(.x)){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  group_by(batch, tf) %>%
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
  mutate(kinase = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05 | is.na(.x)){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  group_by(batch, tf) %>%
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
  mutate(protein = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05 | is.na(.x)){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  select(-pearson_r, -p.value, -n)

tf_rna_reg <- tf_rna_cor %>%
  filter(tf %in% tf_act_var[tf_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  #mutate(rna = if_else(abs(pearson_r) > 0.4, 1, 0)) %>%
  #mutate(rna = if_else(abs(pearson_r) > 0.4 & p.value < 0.05, 1, 0)) %>%
  #mutate(rna = if_else(p.value < 0.05, 1, 0)) %>%
  #mutate(rna = if_else(p.value < 0.05 & n >= quantile(n,0.75), 1, 0)) %>%
  mutate(rna = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05 | is.na(.x)){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  select(-pearson_r, -p.value, -n)

tf_cnv_reg <- tf_cnv_cor %>%
  filter(tf %in% tf_act_var[tf_act_var$var_act_cutoff == 1, "protein", drop = T]) %>%
  #mutate(cnv = if_else(abs(pearson_r) > 0.4, 1, 0)) %>%
  #mutate(cnv = if_else(abs(pearson_r) > 0.4 & p.value < 0.05, 1, 0)) %>%
  #mutate(cnv = if_else(p.value < 0.05, 1, 0)) %>%
  #mutate(cnv = if_else(p.value < 0.05 & n >= quantile(n,0.75), 1, 0)) %>%
  mutate(cnv = map2_dbl(.x = p.value, .y = pearson_r, .f = ~ if(.x > 0.05 | is.na(.x)){0}else{if(.y < 0.2){1}else if(.y >= 0.2 & .y < 0.4){2}else if(.y >= 0.4){3}else{stop(0)}})) %>%
  select(-pearson_r, -p.value, -n)

tf_regulation <- reduce(
  .x = list(tf_phos_reg, tf_kin_reg, tf_prot_reg, tf_rna_reg, tf_cnv_reg),
  .f = ~ full_join(.x, .y, by = c("batch","tf"))) %>%
  pivot_longer(-c(batch, tf), names_to = "data", values_to = "regulated") %>%
  arrange(data) %>%
  mutate(regulated = replace_na(regulated, 4)) %>%
  mutate(prot_type = "tf") %>%
  rename(protein = tf) %>%
  select(prot_type, everything())


regulation_plot <- tf_regulation %>% 
  bind_rows(kinase_regulation) %>%
  mutate(regulated = as.character(regulated)) %>%
  filter(!data %in% c("tf", "kinase")) %>%
  group_by(batch, prot_type, data, regulated) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(batch, prot_type, data) %>%
  mutate(p = n/sum(n)) %>%
  ungroup() %>%
  mutate(data = fct_relevel(data, "cnv", "rna", "protein", "phospho")) %>%
  ggplot(mapping = aes(x = data, y = p, fill = regulated)) +
  geom_col(position = "stack") + 
  facet_grid(prot_type ~ batch,
             labeller = labeller(
               prot_type = c("kinase" = "Kinases", "tf" = "TFs"),
               batch = c("cbttc" = "Brain", "ccle-breast" = "Breast (ccle)","ccle-colorectal" = "Colorectal (ccle)","colon-opportunities" = "Colon", "discovery-ccrcc" = "Kidney", "discovery-luad"="Lung", "discovery-ucec" = "Uterus", "eogc-proteogenomics" = "Stomach", "hcc-proteogenomics" = "Liver", "tcga-brca" = "Breast (tcga)", "tcga-ov" = "Ovary", "tcga-coread" = "Colorectal (tcga)"))) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 18, color = "black"),
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18, color = "black"),
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_blank(),
    legend.position = "bottom") +
  #scale_fill_viridis(discrete = T, name = "Regulated", labels = c("0" = "Non-correlated (P > 5%)", "1" = "Correlated (P < 5% and r < 0.2)", "2" = "Correlated (0.2 <= r < 0.4)", "3" = "Correlated (r >= 0.4)", "4" = "Unknown")) +
  scale_fill_brewer(type = "qual", palette = "Set3", name = "Regulated", labels = c("0" = "Non-correlated (P > 5%)", "1" = "Correlated (P < 5% and r < 0.2)", "2" = "Correlated (0.2 <= r < 0.4)", "3" = "Correlated (r >= 0.4)", "4" = "Unknown")) +
  scale_y_continuous(labels = function(x) paste0(x*100, "%"), name = "Percentage") +
  scale_x_discrete(name = "Data", labels = c("cnv" = "C", "phospho" = "Pho", "protein" = "Pro", "rna" = "R")) +
  guides(fill = guide_legend(nrow = 2))

ggsave(filename = "kin_tf_regulation_perct_batch.png", plot = regulation_plot, path = "./output/plots/kin_tf_regulation/", height = 6, width = 20)
ggsave(filename = "kin_tf_regulation_perct_batch.pdf", plot = regulation_plot, path = "./output/plots/kin_tf_regulation/", height = 6, width = 20)

