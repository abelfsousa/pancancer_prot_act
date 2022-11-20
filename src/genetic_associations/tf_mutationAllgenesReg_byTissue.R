# Estimate the association between TF activity and the mutational status of all genes with mutations
# By tissue


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


# prepare the data to calculate the associations
filter_mut <- function(x, mut_data){
  new_mut <- mut_data %>%
    select_if(.predicate = colnames(.) %in% c("n", "gene", x)) %>%
    pivot_longer(cols = -c(n, gene), names_to = "sample", values_to = "value") %>%
    group_by(gene) %>%
    mutate(n = sum(value == 1)) %>%
    ungroup() %>%
    pivot_wider(names_from = "sample", values_from = "value")
  
  return(new_mut)
}

associations_data <- metadata %>%
  select(-cancer) %>%
  group_by(tissue, batch) %>%
  summarise(sample = list(sample)) %>%
  ungroup() %>%
  mutate(mut = map(.x = sample, .f = filter_mut, mut_data = mut)) %>%
  #mutate(ka = map(.x = sample, .f = ~ ka %>% filter(sample %in% .x)))
  mutate( tf = map( .x = sample, .f = ~ .y %>% filter(sample %in% .x), .y = tf ))


# perform the associations
associations1 <- associations_data %>%
  mutate(res = map2(.x = mut, .y = tf, TfA_mut, samp_mut = 5, mut_tfN = 5)) %>%
  select(-sample, -mut, -tf) %>%
  unnest()

associations2 <- associations_data %>%
  mutate(res = map2(.x = mut, .y = tf, TfA_mut, samp_mut = 5, mut_tfN = 4)) %>%
  select(-sample, -mut, -tf) %>%
  unnest()

associations3 <- associations_data %>%
  mutate(res = map2(.x = mut, .y = tf, TfA_mut, samp_mut = 5, mut_tfN = 3)) %>%
  select(-sample, -mut, -tf) %>%
  unnest()

associations4 <- associations_data %>%
  mutate(res = map2(.x = mut, .y = tf, TfA_mut, samp_mut = 5, mut_tfN = 10)) %>%
  select(-sample, -mut, -tf) %>%
  unnest()

# write associations
write_tsv(associations1, "./output/files/tf_mutStatus_allGenes_byTissue1.txt.gz")
write_tsv(associations2, "./output/files/tf_mutStatus_allGenes_byTissue2.txt.gz")
write_tsv(associations3, "./output/files/tf_mutStatus_allGenes_byTissue3.txt.gz")
write_tsv(associations4, "./output/files/tf_mutStatus_allGenes_byTissue4.txt.gz")
