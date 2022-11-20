# Estimate the association between TF activity and the mutational status of all genes with mutations


library(tidyverse)

source("./src/utils/TfA_mut.R")
source("./src/utils/getSamples.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load samples metadata
metadata <- read_tsv(file = "./output/files/all_samples_annotation.txt")
metadata <- getSamples(metadata, c("protein", "mRNA", "mutation"))


# load mutation matrix
mut <- read_tsv(file = "./output/files/mutations_matrix.txt.gz")


# load TF activity data
tf <- read_csv(file = "./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  rename(tf = X1) %>%
  pivot_longer(-tf, names_to = "sample", values_to = "activity")


# prepare batch covariate
batch <- metadata %>%
  mutate(batch = str_c(batch, tissue, sep = "_")) %>%
  select(sample, batch)


# perform the associations
associations <- TfA_mut(mutMat = mut, tfA = tf, covars = batch)


# write associations
write_tsv(associations, "./output/files/tf_mutStatus_allGenes.txt.gz")
