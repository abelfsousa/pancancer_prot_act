library(tidyverse)



# load RPPA phosphorylation data
load("./data/protein/tcga/rppa/pan_can_tcga_rppa_mat_updated.Rdata")
pan_can_tcga_rppa_mat <- pan_can_tcga_rppa_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  as_tibble() %>%
  #select_at(.vars = c(1, which(str_sub(colnames(.),14,15) == "01"))) %>%
  select(feature, ends_with("-01")) %>%
  pivot_longer(-feature, names_to = "sample", values_to = "expr") %>%
  mutate(sample = str_replace_all(str_sub(sample, 1, 12), "-", "."))


load("./data/protein/tcga/rppa/matched_prot_phosphoprot_list_updated.Rdata")
matched_prot_phosphoprot_list <- matched_prot_phosphoprot_list %>%
  as_tibble()

rppa_pho <- pan_can_tcga_rppa_mat %>%
  inner_join(matched_prot_phosphoprot_list, by = c("feature" = "fil_rppa_phospho_ab")) %>%
  select(psite = feature, protein = fil_rppa_prot_ab, gene = fil_gene_prot_list, sample, psite_expr = expr) %>%
  mutate(across(.cols = where(is.factor), .fns = as.character))


# load kinase-activity inference data
# (quantile-normalized protein regressed-out phosphorylation data)
# select quantifications with more than 3 substrates
k_subN <- 3
kins <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  select(kinase, sample, kin_activity=log10P)

source("src/utils/getSamples.R")
metadata <- read_tsv(file = "./output/files/all_samples_annotation.txt")
metadata <- getSamples(metadata, c("phosphorylation", "protein"))


sample_rppa <- rppa_pho %>%
  select(sample) %>%
  distinct()

sample_cptac <- kins %>%
  select(sample) %>%
  distinct()

sample_cptac %>%
  inner_join(metadata, by = "sample") %>%
  filter(sample %in% sample_rppa$sample) %>%
  #semi_join(sample_rppa, by = "sample") %>%
  group_by(batch, cancer, tissue) %>%
  summarise(n = n(), s = list(sample)) %>%
  ungroup()

sample_rppa %>%
  summarise(x = n(), y = sum(sample %in% sample_cptac$sample)) %>%
  mutate(z = (y/x))


load("data/Danish/tissue_type_sel.Rdata")
tissue_type_sel %>%
  names() %>%
  str_replace_all("-", ".") %>%
  str_sub(1, 12) %>%
  unique() %>%
  as_tibble() %>%
  rename(sample = value) %>%
  summarise(x = n(), y = sum(sample %in% sample_cptac$sample)) %>%
  mutate(z = (y/x))

#https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations
tissue_type_sel %>%
  unname() %>%
  as_tibble() %>%
  distinct() %>%
  mutate(
    tissue = c("bladder", "breast", "uterus", "colorectal", "head_neck", "kidney", "kidney", "brain", "liver", "lung", "lung", "ovary", "prostate", "bone", "skin", "stomach", "thyroid")) %>%
  filter(tissue %in% unique(metadata$tissue))

tissue_type_sel %>%
  unname() %>%
  as_tibble() %>%
  distinct() %>%
  mutate(
    tissue = c("bladder", "breast", "uterus", "colorectal", "head_neck", "kidney", "kidney", "brain", "liver", "lung", "lung", "ovary", "prostate", "bone", "skin", "stomach", "thyroid")) %>%
  filter(!tissue %in% unique(metadata$tissue))
