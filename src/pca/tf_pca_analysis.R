# PCA on the TF activities

# load R packages
library(RColorBrewer)
library(tidyverse)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load samples annotation
source("./src/utils/getSamples.R")
samples_annotation <- read_tsv(file = "./output/files/all_samples_annotation.txt")
samples_annotation <- getSamples(samples_annotation, c("protein", "mRNA")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue) %>%
  mutate(batch = str_replace(batch, "tcga-brca", "CPTAC-breast\n(TCGA)")) %>%
  mutate(batch = str_replace(batch, "tcga-ov", "CPTAC-ovarian\n(TCGA)")) %>%
  mutate(batch = str_replace(batch, "tcga-coread", "CPTAC-colorectal\n(TCGA)")) %>%
  mutate(batch = str_replace(batch, "cbttc", "CPTAC-brain")) %>%
  mutate(batch = str_replace(batch, "discovery-ccrcc", "CPTAC-kidney\n(discovery)")) %>%
  mutate(batch = str_replace(batch, "discovery-luad", "CPTAC-lung\n(discovery)")) %>%
  mutate(batch = str_replace(batch, "discovery-ucec", "CPTAC-uterus\n(discovery)")) %>%
  mutate(batch = str_replace(batch, "colon-opportunities", "CPTAC-colon\n(confirmatory)")) %>%
  mutate(batch = str_replace(batch, "eogc-proteogenomics", "CPTAC-stomach")) %>%
  mutate(batch = str_replace(batch, "hcc-proteogenomics", "CPTAC-liver")) %>%
  mutate(batch = str_replace(batch, "ccle-colorectal", "CCLE-colorectal")) %>%
  mutate(batch = str_replace(batch, "ccle-breast", "CCLE-breast"))


# load TF activities
tf_matrix <- read_csv(file = "./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  as.data.frame() %>%
  column_to_rownames(var = "X1") %>%
  t()


# perform PCA
pca <- prcomp(tf_matrix, scale. = T, center = T)

tf_matrix_pcs <- pca %>% 
  pluck("x") %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  as_tibble()

var_expl <- pca %>%
  pluck("sdev") %>%
  as_tibble() %>%
  rename(sdev = value) %>%
  mutate(PC = str_c("PC", seq_len(nrow(.)))) %>%
  select(PC, sdev) %>%
  mutate(var = sdev^2) %>%
  mutate(var_perc = (var/sum(var))*100)

var_barplot <- var_expl %>%
  mutate(PC = fct_reorder(PC, var_perc, function(x) x)) %>%
  slice(1:20) %>%
  ggplot(mapping = aes(x = PC, y = var_perc)) +
  geom_col(fill = "#2171b5") +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title.x = element_text(colour = "black", size = 18),
    axis.text.x = element_text(colour = "black", size = 14),
    axis.title.y = element_blank(),
    axis.text.y = element_text(colour = "black", size = 16)) +
  labs(x = "Principal components", y = "% of variance")

ggsave(filename = "tf_pcs_var_explained.png", plot = var_barplot, path = "./output/plots/pca/", width = 2, height = 6)
ggsave(filename = "tf_pcs_var_explained.pdf", plot = var_barplot, path = "./output/plots/pca/", width = 2, height = 6)


pca_plot <- tf_matrix_pcs %>%
  inner_join(samples_annotation, by = "sample") %>%
  ggplot(mapping = aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_brewer(type = "qual", palette = "Paired") +
  theme_classic() +
  theme(
    #plot.margin = unit(c(0.1, 0.1, 5, 0.1), "cm"),
    plot.title = element_text(colour = "black", size = 18, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 18),
    axis.text = element_text(colour = "black", size = 16),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 14),
    #legend.position = c(0,0),
    #legend.justification = c(0.05,1.3),
    #legend.margin = margin(c(-0.2,0,0,0), unit="cm"),
    #legend.text = element_text(colour = "black", size = 12, margin = margin(t = 5))
    legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1), ncol = 3)) +
  labs(x = "PC1", y = "PC2", title = "TF activity PCA", color = "Study")

ggsave(filename = "tf_pca_sample_projection.png", plot = pca_plot, path = "./output/plots/pca/", width = 6, height = 6)
ggsave(filename = "tf_pca_sample_projection.pdf", plot = pca_plot, path = "./output/plots/pca/", width = 6, height = 6)
