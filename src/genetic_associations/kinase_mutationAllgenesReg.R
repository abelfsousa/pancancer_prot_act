# Estimate the association between kinase-activity and the mutational status of all genes with mutations


library(tidyverse)

source("./src/utils/KA_mut.R")
source("./src/utils/getSamples.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load samples metadata
metadata <- read_tsv(file = "./output/files/all_samples_annotation.txt")
metadata <- getSamples(metadata, "phosphorylation")


# load mutation matrix
mut <- read_tsv(file = "./output/files/mutations_matrix.txt.gz")


# load kinase-activity inference data
# (quantile-normalized protein regressed-out phosphorylation data)
ka <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz")


# prepare batch covariate
batch <- metadata %>%
  mutate(batch = str_c(batch, tissue, sep = "_")) %>%
  select(sample, batch)


# perform the associations
associations <- KA_mut(mutMat = mut, ka = ka, covars = batch)


# write associations
write_tsv(associations, "./output/files/ka_mutStatus_allGenes.txt.gz")
