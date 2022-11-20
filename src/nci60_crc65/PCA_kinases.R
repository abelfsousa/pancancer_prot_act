library(tidyverse)
library(caret)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load CPTAC/CCLE samples
cptac_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
  filter(data == "phosphorylation")


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


# load kinase activity inference data from NCI60/CRC65 dataset
# select quantifications with more than 3 substrates
nci60_crc65_KA <- read_tsv(file = "./output/files/nci60_crc65_kinase_activities.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= 3) %>%
  select(-source_type, -n) %>%
  filter(!(batch == "CRC65" & sample %in% shared_cell_lines)) %>%
  select(batch, sample, kinase, kin_activity=log10P)


# load kinase activity inference data without imputations from CPATC/CCLE dataset
# (quantile-normalized protein regressed-out phosphorylation data)
# select quantifications with more than 3 substrates
ccle_KA <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= 3) %>%
  select(-source_type, -n) %>%
  select(sample, kinase, kin_activity=log10P) %>%
  inner_join(cptac_samples[, c("sample", "batch")], by = "sample") %>%
  filter(batch == "ccle") %>%
  anti_join(nci60_crc65_KA, by = "sample") %>%
  mutate(batch = "Roumeliotis") #%>%
  #anti_join(alternative_names, by = c("sample" = "alternative_names"))


# bind the kinase activities of both datasets
kin_activity <- bind_rows(ccle_KA, nci60_crc65_KA)

batch <- kin_activity %>%
  select(batch, sample) %>%
  distinct()

kin_mat <- kin_activity %>%
  group_by(kinase) %>%
  filter(n() >= 145*0.8) %>%
  ungroup() %>%
  select(-batch) %>%
  pivot_wider(names_from = "sample", values_from = "kin_activity") %>%
  as.data.frame() %>%
  column_to_rownames(var = "kinase") %>%
  as.matrix() %>%
  t()

pre_process <- preProcess(kin_mat, method = c("bagImpute"))
kin_mat <- predict(pre_process, kin_mat)

pca <- prcomp(kin_mat)

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
  