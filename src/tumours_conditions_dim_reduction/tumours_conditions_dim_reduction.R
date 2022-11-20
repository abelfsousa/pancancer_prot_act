# load R packages
library(tidyverse)
library(umap)


# load kinase activities with imputed values
kin_act_imp_v1 <- read_tsv("./output/files/kin_act_tumours_perturb_imputed_v1.txt")
kin_act_imp_v2 <- read_tsv("./output/files/kin_act_tumours_perturb_imputed_v2.txt")


# load metadata of human perturbations
metadata_pertub <- read_tsv("./output/files/esetNR_cond_anno.txt")


# load metadata of tumor samples
metadata_tumors <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
  filter(data == "phosphorylation") %>%
  select(-data) %>%
  #mutate(batch = if_else(batch == "ccle", "ccle-colorectal", batch))
  mutate(batch = map2_chr(batch, tissue, ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x}))


# sample types
perturbations <- metadata_pertub %>%
  select(sample, description) %>%
  mutate(type = "perturbations")

tumours <- metadata_tumors %>%
  select(sample, description = batch) %>%
  mutate(type = "tumours")

sample_type <- bind_rows(perturbations, tumours)


# set up a matrix for dimensionality reduction
# remove PGE2 perturbations
kin_act_matrix <- kin_act_imp_v1 %>%
  inner_join(sample_type, by = "sample") %>%
  filter(description != "PGE2") %>%
  select(-description, -type) %>%
  pivot_wider(names_from = "protein", values_from = "activity") %>%
  column_to_rownames(var = "sample") %>%
  as.matrix()


# perform UMAP analysis
kin_act_umap <- umap(kin_act_matrix)


# perform PCA analysis
kin_act_pca <- prcomp(kin_act_matrix)


# plot UMAP
change_description <- function(x){
  label <- NA
  if(str_detect(x, "DNA damage")){
    label <- "DNA damage"
  } else if(str_detect(x, "Mitosis")){
    label <- "Mitosis"
  } else if(str_detect(x, "EGF")){
    label <- "EGF"
  }
  return(label)
}

umap_reduction <- kin_act_umap %>%
  pluck("layout") %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  as_tibble() %>%
  rename(UMAP1 = V1, UMAP2 = V2) %>%
  inner_join(sample_type, by = "sample") %>%
  mutate(description2 = map_chr(description, change_description)) %>%
  mutate(description2 = str_replace_na(description2, "Other"))

umap_plot <- ggplot(data = umap_reduction, mapping = aes(x = UMAP1, y = UMAP2, color = description2)) +
  geom_point() +
  theme_classic() +
  theme(
    legend.position = "right") +
  scale_color_manual(name = "Description", values = c("Other" = "grey", "DNA damage" = "red", "EGF" = "chartreuse3", "Mitosis" = "blue"))

ggsave(filename = "./output/plots/tumours_conditions_dim_reduction/tumours_conditions_umap.png", plot = umap_plot, width = 8, height = 6)


# correlate UMAP projections with kinase activities
kinase_umap_corr <- umap_reduction %>%
  select(sample, UMAP1, UMAP2) %>%
  pivot_longer(-sample, names_to = "UMAP", values_to = "value") %>%
  inner_join(kin_act_imp_v1, by = "sample") %>%
  group_by(protein, UMAP) %>%
  summarise(corr = cor(activity, value)) %>%
  #nest() %>%
  ungroup() %>%
  #mutate(corr = map_dbl(data, ~ cor.test(.x$value, .x$activity)$estimate)) %>%
  #select(-data) %>%
  pivot_wider(names_from = "UMAP", values_from = "corr") %>%
  arrange(protein)

# kinase_umap_corr <- cor(kin_act_matrix, as.matrix(setNames(as.data.frame(kin_act_umap$layout), c("UMAP1", "UMAP2")))) %>%
#   as.data.frame() %>%
#   rownames_to_column("protein") %>%
#   as_tibble() %>%
#   arrange(protein)


# plot PCA
var_expl <- kin_act_pca %>%
  pluck("sdev") %>%
  as_tibble() %>%
  rename(sdev = value) %>%
  mutate(PC = str_c("PC", seq_len(nrow(.)))) %>%
  select(PC, sdev) %>%
  mutate(var = sdev^2) %>%
  mutate(var_perc = (var/sum(var))*100)

pca_reduction <- kin_act_pca %>%
  pluck("x") %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  as_tibble() %>%
  inner_join(sample_type, by = "sample")

pca_plot1 <- ggplot(data = pca_reduction, mapping = aes(x = PC1, y = PC2, color = type)) +
  geom_point()
