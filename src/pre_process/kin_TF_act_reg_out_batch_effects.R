library(tidyverse)

source("./src/utils/getSamples.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load samples annotation
samples_annotation <- read_tsv(file = "./output/files/all_samples_annotation.txt")
samples_annotation <- getSamples(samples_annotation, c("protein"))

samples_annotation <- samples_annotation %>%
  mutate(batch = str_c(batch, tissue, sep = "-"))


# load transcription factor activities
#tfs <- data.table::fread("./data/progeny/TF_activity_log2FC.csv") %>%
tfs <- data.table::fread("./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  as_tibble() %>%
  rename(tf=V1) %>%
  pivot_longer(-tf, names_to = "sample", values_to = "tf_activity")


# load kinase activities and imputed values
# matrix used for PCA analysis and UMAP
kin_activities <- data.table::fread("./data/Danish/kinaseActMatImputed.tsv") %>%
  as_tibble() %>%
  rename(kinase = V1) %>%
  pivot_longer(-kinase, names_to = "sample", values_to = "kin_activity")


# regress-out batch effects from TFs and kinase activities
get_residuals <- function(df){
  m <- lm(activity ~ batch, data = df)
  resd <- residuals(m)
  
  resd
}


# TFs
tfs_regout <- tfs %>%
  rename(activity = tf_activity) %>%
  inner_join(samples_annotation[, c("sample", "batch")], by = "sample") %>%
  group_by(tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(activity_res = map(.x = data, .f = get_residuals)) %>%
  unnest() %>%
  select(-activity, -batch) %>%
  pivot_wider(names_from = "sample", values_from = "activity_res")

write_tsv(tfs_regout, "./output/files/TF_activity_log2FC_batch_regOut.txt")


# kinases
kin_regout <- kin_activities %>%
  rename(activity = kin_activity) %>%
  inner_join(samples_annotation[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(activity_res = map(.x = data, .f = get_residuals)) %>%
  unnest() %>%
  select(-activity, -batch) %>%
  pivot_wider(names_from = "sample", values_from = "activity_res")

write_tsv(kin_regout, "./output/files/kinaseActMatImputed_batch_regOut.txt")
