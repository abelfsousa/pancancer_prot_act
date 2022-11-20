# Estimate the association between kinase-activity principal components and the mutational status of all genes with mutations


library(tidyverse)
library(parallel)

source("./src/utils/KA_mut.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy



# load mutation matrix
mut <- read_tsv(file = "./output/files/mutations_matrix.txt.gz") %>%
  filter(n >= 50) %>%
  select(-n) %>%
  pivot_longer(cols = -c(gene), names_to = "sample", values_to = "mutated")


# load kinase-activity principal components
kaPCs <- data.table::fread(file = "./data/Danish/kinasePCAMat.tsv") %>%
  as_tibble() %>%
  rename(sample = V1) %>%
  select(1:11) %>%
  pivot_longer(-sample, names_to = "prin_comp", values_to = "value")


# load metadata
metadata <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
  filter(data == "phosphorylation") %>%
  select(-data) %>%
  select(sample, batch, tissue)


# prepare data to model
dataMod <- kaPCs %>%
  inner_join(mut, by = c("sample")) %>%
  select(prin_comp, gene, sample, value, mutated) %>%
  rename(log10P = value) %>%
  group_by(prin_comp, gene) %>%
  nest() %>%
  ungroup()



# perform the associations

# single core run
# associations <- dataMod %>%
#   mutate(model = map(.x = data, .f = linear_model, covs = metadata, mutN = 5)) %>%
#   select(-data) %>%
#   unnest() %>%
#   filter(!is.na(estimate)) %>%
#   mutate(p.adjust = p.adjust(p.value, "BH"))

# multi core run
cores = detectCores()-1

models <- mclapply(
  X = dataMod$data,
  FUN = linear_model,
  covs = metadata[, c("sample", "batch")],
  mutN = 5,
  mc.cores = cores)

associations <- dataMod %>%
  mutate(model = models) %>%
  select(-data) %>%
  unnest() %>%
  filter(!is.na(estimate)) %>%
  mutate(p.adjust = p.adjust(p.value, "BH"))


# export data
write_tsv(x = associations, path = "./output/files/kaPCs_mutStatus_allGenes.txt.gz")
