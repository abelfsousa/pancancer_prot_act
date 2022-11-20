library(tidyverse)
library(caret)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load CPTAC/CCLE samples
source("./src/utils/getSamples.R")
cptac_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")
cptac_samples <- getSamples(cptac_samples, c("protein","mRNA"))

cell_lines_batch <- read_tsv(file = "./output/files/proteomics_samples.txt")


# load NCI60/CRC65 cell lines
nci60_crc65 <- read_tsv(file = "./output/files/nci60_crc65_cell_lines.txt")

shared_cell_lines <- nci60_crc65 %>%
  group_by(batch) %>%
  summarise(cell_line = list(cell_line)) %>%
  pull(cell_line) %>%
  reduce(intersect)

alternative_names <- nci60_crc65 %>%
  filter(!is.na(alternative_names)) %>%
  mutate(alternative_names = str_split(alternative_names, ";")) %>%
  unnest() %>%
  distinct()


# load TF activity inference data from NCI60/CRC65 dataset
crc65_TF <- read_csv(file = "./data/dorothea/nci60_crc65/TF_activity_log2FC_allsamples_CRC65.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "tf_activity") %>%
  rename(tf=X1) %>%
  filter(!sample %in% shared_cell_lines) %>%
  select(sample, tf, tf_activity) %>%
  mutate(batch = "CRC65")

nci60_TF <- read_csv(file = "./data/dorothea/nci60_crc65/TF_activity_log2FC_allsamples_NCI60.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "tf_activity") %>%
  rename(tf=X1) %>%
  select(sample, tf, tf_activity) %>%
  mutate(batch = "NCI60")

nci60_crc65_TF <- bind_rows(crc65_TF, nci60_TF)


# load TF activity inference data from CPATC/CCLE dataset
ccle_TF <- read_csv(file = "./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "tf_activity") %>%
  rename(tf=X1) %>%
  select(sample, tf, tf_activity) %>%
  inner_join(cptac_samples[, c("sample", "batch")], by = "sample") %>%
  filter(batch == "ccle") %>%
  anti_join(nci60_crc65_TF, by = "sample") %>%
  select(-batch) %>%
  inner_join(cell_lines_batch[, c("sample", "batch")], by = "sample")
#anti_join(alternative_names, by = c("sample" = "alternative_names"))


# bind the TF activities of both datasets
tf_activity <- bind_rows(ccle_TF, nci60_crc65_TF)

batch <- tf_activity %>%
  select(batch, sample) %>%
  distinct()

tf_mat <- tf_activity %>%
  select(-batch) %>%
  pivot_wider(names_from = "sample", values_from = "tf_activity") %>%
  as.data.frame() %>%
  column_to_rownames(var = "tf") %>%
  as.matrix() %>%
  t()

pca <- prcomp(tf_mat)

plot <- pca %>%
  pluck("x") %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  as_tibble() %>%
  select(1:11) %>%
  inner_join(batch, by = "sample") %>%
  ggplot(mapping = aes(x = PC1, y = PC2, color = batch)) +
  geom_point()
plot  
