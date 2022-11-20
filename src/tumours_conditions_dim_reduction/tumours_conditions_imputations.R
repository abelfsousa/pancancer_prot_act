# load R packages
library(missForest)
library(tidyverse)


# load kinase activities from tumour samples
kin_act_tumors <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= 3) %>%
  select(-source_type, -n, activity = log10P, protein = kinase)


# load kinase activities from human pertubations
kin_act_pertub <- read_tsv(file = "./output/files/KA_esetNR.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= 3) %>%
  select(-source_type, -n, activity = log10P, protein = kinase)


# load metadata of human perturbations
metadata_pertub <- read_tsv("./output/files/esetNR_cond_anno.txt")


# load metadata of tumor samples
metadata_tumors <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
  filter(data == "phosphorylation") %>%
  select(-data)


# set up a function to extract the kinases and samples quantified in at least X% and X% of samples
get_features <- function(df, perct_kinases, perct_samples){
  
  new_df <- df %>%
    pivot_wider(names_from = "protein", values_from = "activity") %>%
    pivot_longer(cols = -sample, names_to = "protein", values_to = "activity")
  
  kinases <- new_df %>%
    group_by(protein) %>%
    summarise(selected = if_else(sum(!is.na(activity)) >= n()*perct_kinases, 1, 0)) %>%
    ungroup() %>%
    filter(selected == 1) %>%
    pull(protein)
  
  samples <- new_df %>%
    filter(protein %in% kinases) %>%
    group_by(sample) %>%
    summarise(selected = if_else(sum(!is.na(activity)) >= n()*perct_samples, 1, 0)) %>%
    ungroup() %>%
    filter(selected == 1) %>%
    pull(sample)
  
  l <- list(kinases, samples)
  
  return(l)
}


# select kinases and samples quantified in at least 60% and 70% of the samples and kinases, respectively
# impute missing values using a random forest

tumors_selection <- get_features(kin_act_tumors, 0.6, 0.7)
pertub_selection <- get_features(kin_act_pertub, 0.6, 0.7)

tumors_samples <- tumors_selection[[2]]
pertub_samples <- pertub_selection[[2]]
common_kinases <- intersect(tumors_selection[[1]], pertub_selection[[1]])

# -- strategy 1: impute missing values by dataset
kin_act_tumors_sel <- kin_act_tumors %>%
  filter(protein %in% common_kinases & sample %in% tumors_samples) %>%
  pivot_wider(names_from = "protein", values_from = "activity") %>%
  column_to_rownames(var = "sample") %>%
  as.matrix()

kin_act_pertub_sel <- kin_act_pertub %>%
  filter(protein %in% common_kinases & sample %in% pertub_samples) %>%
  pivot_wider(names_from = "protein", values_from = "activity") %>%
  column_to_rownames(var = "sample") %>%
  as.matrix()

kin_act_tumors_imp <- missForest(kin_act_tumors_sel)
kin_act_pertub_imp <- missForest(kin_act_pertub_sel)

# merge matrices and export
kin_act_tumors_imp <- kin_act_tumors_imp$ximp %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(-sample, names_to = "protein", values_to = "activity")

kin_act_pertub_imp <- kin_act_pertub_imp$ximp %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(-sample, names_to = "protein", values_to = "activity")

kin_act_merged <- bind_rows(kin_act_tumors_imp, kin_act_pertub_imp)
write_tsv(kin_act_merged, "./output/files/kin_act_tumours_perturb_imputed_v1.txt")


# -- strategy 2: impute missing across both datasets
kin_act_tumors_pertub_sel <- bind_rows(kin_act_tumors, kin_act_pertub) %>%
  filter(protein %in% common_kinases & (sample %in% pertub_samples | sample %in% tumors_samples)) %>%
  pivot_wider(names_from = "protein", values_from = "activity") %>%
  column_to_rownames(var = "sample") %>%
  as.matrix()

kin_act_tumors_pertub_imp <- missForest(kin_act_tumors_pertub_sel)

kin_act_tumors_pertub_imp <- kin_act_tumors_pertub_imp$ximp %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(-sample, names_to = "protein", values_to = "activity")

write_tsv(kin_act_tumors_pertub_imp, "./output/files/kin_act_tumours_perturb_imputed_v2.txt")

