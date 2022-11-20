#' ---
#' title: "UMAP analysis on the TF activities"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=12, fig.height=6)


#' load R packages
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(ggrepel)
library(umap)

source("./src/utils/KA_mut.R")

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load CPTAC TF activities
tf_matrix <- read_csv(file = "./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  as.data.frame() %>%
  column_to_rownames(var = "X1") %>%
  t()


#' load samples annotation
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


#' ## create UMAP projection (TFs as variables)
set.seed(123)
tf_umap1 = umap(tf_matrix)

projections <- tf_umap1 %>%
  pluck("layout") %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  as_tibble() %>%
  rename(UMAP1 = V1, UMAP2 = V2) %>%
  arrange(sample)
write_tsv(projections, "./output/files/tfActivities_umap_projections.txt")


#' ## plot UMAP1 vs UMAP2 colored by batch
umap_plot <- projections %>%
  inner_join(samples_annotation, by = "sample") %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, color = batch)) +
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
  labs(x = "UMAP1", y = "UMAP2", title = "TF activity UMAP", color = "Study")

#+ fig.width=4, fig.height=6
umap_plot

ggsave(filename = "tf_umap_sample_projection.png", plot = umap_plot, path = "./output/plots/umap/", width = 6, height = 6)
ggsave(filename = "tf_umap_sample_projection.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 6, height = 6)


#' ## correlate each TF with the UMAP projections
tfs <- tf_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  as_tibble() %>%
  arrange(sample)

cor_mat <- cor(tfs[, -c(1)], projections[, -c(1)])


#' ## remove TF redundancies
source(file = "./data/Danish/filterRedundantTerms.R")
tf_targets <- read_csv(file = "./data/dorothea/dorothea_ABC.csv") %>%
  select(-confidence, -mor) %>%
  group_by(tf) %>%
  summarise(target = list(target), n = n()) %>%
  ungroup() %>%
  filter(n >= 3) %>%
  select(-n) %>%
  unnest() %>%
  rename(ont=tf, gene=target) %>%
  as.matrix()

# select non-redundant TFs
nr_tfs <- pruneGOTerms(tf_targets, category = "broad", numTerms = 1, treeHeight = 0.99)

sel <- cor_mat %>%
  abs() %>%
  as.data.frame() %>%
  rownames_to_column(var = "tf") %>%
  as_tibble() %>%
  filter(tf %in% nr_tfs) %>%
  filter(UMAP1 > 0.45 | UMAP2 > 0.45)

cor_mat2 <- cor_mat[rownames(cor_mat) %in% sel$tf, ]


#' ## heatmap of the correlations
heatmap <- Heatmap(
  column_title = "UMAP",
  column_title_gp = gpar(fontsize = 16),
  col = c("#d8b365", "#f5f5f5", "#5ab4ac"),
  show_heatmap_legend = T,
  border = T,
  column_labels = c("1", "2"),
  column_names_rot = 0,
  rect_gp = gpar(col = "white", lwd = 0.5),
  column_names_centered = T,
  matrix = cor_mat2,
  name = "Pearson's r",
  heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter", labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize = 17)),
  cluster_columns = F,
  row_names_gp = gpar(fontsize = 11),
  column_names_gp = gpar(fontsize = 16))

#+ fig.width=2, fig.height=6
draw(heatmap, heatmap_legend_side = "bottom")

pdf(file = "./output/plots/umap/umap_tfs_cor.pdf", height = 7, width = 3)
draw(heatmap, heatmap_legend_side = "bottom")
dev.off()

png(file="./output/plots/umap/umap_tfs_cor.png", units = "in", res = 300, height=7, width=3)
draw(heatmap, heatmap_legend_side = "bottom")
dev.off()


#' ## plot UMAP1 vs UMAP2 colored by TFs activity
umap_plot <- projections %>%
  inner_join(tfs, by = "sample") %>%
  select_at(.vars = c("sample", "UMAP1", "UMAP2", sel$tf)) %>%
  pivot_longer(cols = -c(sample, UMAP1, UMAP2), names_to = "tf", values_to = "activity") %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, color = activity)) +
  geom_point(size = 0.2) +
  facet_wrap(~ tf, ncol = 10) +
  scale_colour_gradient2(low = "#ca0020", mid = "white", high = "#0571b0", na.value = "white", name = "TF activity", limits = c(-6, 6)) +
  #scale_colour_gradient(low = "#de2d26", high = "#2c7fb8") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 11),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14),
    legend.position = "bottom",
    legend.key.width = unit(0.7, "cm"),
    legend.title = element_text(colour = "black", size = 16),
    legend.text = element_text(colour = "black", size = 14))

#+ fig.width=10, fig.height=6
umap_plot

ggsave(filename = "umap_tfs_activity.png", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
ggsave(filename = "umap_tfs_activity.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)


#' ## create UMAP projection (samples as variables)
tf_umap2 = umap(t(tf_matrix))


#' ## plot UMAP1 vs UMAP2
umap_plot <- tf_umap2 %>%
  pluck("layout") %>%
  as.data.frame() %>%
  rownames_to_column(var = "tf") %>%
  as_tibble() %>%
  rename(UMAP1 = V1, UMAP2 = V2) %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, label = tf)) +
  geom_point(size = 1, color = "darkblue") +
  geom_text_repel(size = 2, segment.size = 0.2) +
  theme_minimal()

#+ fig.width=10, fig.height=10
umap_plot

ggsave(filename = "umap_tf_projection.png", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 10)
ggsave(filename = "umap_tf_projection.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 10)
unlink("umap_tf_projection.png")
unlink("umap_tf_projection.pdf")

