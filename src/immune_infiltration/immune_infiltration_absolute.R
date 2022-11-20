#' ---
#' title: "Immune infiltration  (using absolute cibersort estimates)"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=12, fig.height=6)


#' load R packages
library(ComplexHeatmap)
library(tidyverse)


source("./src/utils/getSamples.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load samples annotation
samples_annotation <- read_tsv(file = "./output/files/all_samples_annotation.txt")
metadata <- getSamples(samples_annotation, c("mRNA", "protein"))


#' load immune cells estimates (absolute values)\
#' estimated by Danish using CIBERSORT
cibersort <- read_csv(file = "./data/cibersort/CIBERSORT.CPTAC.absolute.csv") %>%
  rename(sample = X1) %>%
  #filter(str_detect(sample, "^X{1}"))
  mutate(sample = str_replace(sample, "^X{1}", "")) %>%
  select(-`P-value`, -Correlation, -RMSE, -`Absolute score (no.sumto1)`) %>%
  pivot_longer(-sample, names_to = "cell", values_to = "absolute")


#' load kinase activities and imputed values\
#' matrix used for PCA analysis and UMAP
kin_activities <- data.table::fread("./data/Danish/kinaseActMatImputed.tsv") %>%
  as.data.frame() %>%
  column_to_rownames(var = "V1") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  as_tibble() %>%
  pivot_longer(-sample, names_to = "kinase", values_to = "activity")


#' load kinase-activity principal components\
#' select 1:10 principal components
pcs <- data.table::fread(file = "./data/Danish/kinasePCAMat.tsv") %>%
  as_tibble() %>%
  rename(sample = V1) %>%
  select(1:11) %>%
  pivot_longer(-sample, names_to = "pc", values_to = "value")


#' load UMAP projections
umap <- read_tsv("./output/files/kinActivities_umap_projections.txt") %>%
  pivot_longer(-sample, names_to = "umap", values_to = "value")


#' # correlate cibersort relative values with UMAP projections
cor_mat <- cibersort %>%
  inner_join(umap, by = "sample") %>%
  group_by(cell, umap) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x=data, .f = ~ broom::tidy(cor.test(.x$absolute, .x$value)) %>% select(estimate, p.value))) %>%
  select(-data) %>%
  unnest()

correlations <- cor_mat %>%
  select(-p.value) %>%
  pivot_wider(names_from = "umap", values_from = "estimate") %>%
  as.data.frame() %>%
  column_to_rownames(var = "cell") %>%
  as.matrix()

pvalues <- cor_mat %>%
  select(-estimate) %>%
  pivot_wider(names_from = "umap", values_from = "p.value") %>%
  as.data.frame() %>%
  column_to_rownames(var = "cell") %>%
  as.matrix()

heatmap <- Heatmap(
  matrix = correlations,
  name = "Pearson's r",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(pvalues[i, j] < 0.05)
      grid.text("*", x, y, gp = gpar(fontsize = 10))
  })

#+ fig.width=4, fig.height=4
heatmap

pdf(file = "./output/plots/immune_infiltration/umap_immune_absolute_cor.pdf", width = 4, height = 4)
heatmap
dev.off()


#' ## plot UMAP1 vs UMAP2 colored by immune cell estimates
sel_cells <- cor_mat %>%
  filter(abs(estimate) > 0.1 & p.value < 0.05) %>%
  pull(cell) %>%
  unique()

umap_plot <- umap %>%
  pivot_wider(names_from = "umap", values_from = "value") %>%
  inner_join(cibersort, by = "sample") %>%
  filter(cell %in% sel_cells) %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, color = absolute)) +
  geom_point(size = 0.5) +
  facet_wrap(~ cell) +
  theme_minimal() +
  scale_color_gradient(low = "white", high = "red", name = "estimate")

#+ fig.width=8, fig.height=4
umap_plot

ggsave(filename = "umap_projection_immune_absolute_cor.pdf", plot = umap_plot, path = "./output/plots/immune_infiltration/", width = 8, height = 4)
unlink("umap_projection_immune_absolute_cor.pdf")


#' # correlate cibersort relative values with principal components
cor_mat <- cibersort %>%
  inner_join(pcs, by = "sample") %>%
  group_by(cell, pc) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x=data, .f = ~ broom::tidy(cor.test(.x$absolute, .x$value)) %>% select(estimate, p.value))) %>%
  select(-data) %>%
  unnest()

correlations <- cor_mat %>%
  select(-p.value) %>%
  pivot_wider(names_from = "pc", values_from = "estimate") %>%
  as.data.frame() %>%
  column_to_rownames(var = "cell") %>%
  as.matrix()

pvalues <- cor_mat %>%
  select(-estimate) %>%
  pivot_wider(names_from = "pc", values_from = "p.value") %>%
  as.data.frame() %>%
  column_to_rownames(var = "cell") %>%
  as.matrix()

heatmap <- Heatmap(
  matrix = correlations,
  name = "Pearson's r",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(pvalues[i, j] < 0.05)
      grid.text("*", x, y, gp = gpar(fontsize = 10))
  })

#+ fig.width=6, fig.height=4
heatmap

pdf(file = "./output/plots/immune_infiltration/pcs_immune_absolute_cor.pdf", width = 6, height = 4)
heatmap
dev.off()


#' # correlate cibersort relative values with kinase activities\
cor_mat <- cibersort %>%
  inner_join(kin_activities, by = "sample") %>%
  group_by(cell, kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x=data, .f = ~ broom::tidy(cor.test(.x$absolute, .x$activity)) %>% select(estimate, p.value))) %>%
  select(-data) %>%
  unnest()

correlations <- cor_mat %>%
  select(-p.value) %>%
  pivot_wider(names_from = "kinase", values_from = "estimate") %>%
  as.data.frame() %>%
  column_to_rownames(var = "cell") %>%
  as.matrix() %>%
  t()

pvalues <- cor_mat %>%
  select(-estimate) %>%
  pivot_wider(names_from = "kinase", values_from = "p.value") %>%
  as.data.frame() %>%
  column_to_rownames(var = "cell") %>%
  as.matrix() %>%
  t()


#' ## all kinases (90)
heatmap <- Heatmap(
  matrix = correlations,
  name = "Pearson's r",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(pvalues[i, j] < 0.05)
      grid.text("*", x, y, gp = gpar(fontsize = 10))
  })

#+ fig.width=8, fig.height=12
heatmap

pdf(file = "./output/plots/immune_infiltration/kin_act_immune_absolute_cor_all_kinases.pdf", width = 8, height = 12)
heatmap
dev.off()


#' ## select most variable kinases (52; UMAP absolute correlation > 0.4)
kin_umap_cor <- kin_activities %>%
  inner_join(umap, by = "sample") %>%
  group_by(kinase, umap) %>%
  summarise(r = cor(value, activity)) %>%
  ungroup() %>%
  pivot_wider(names_from = "umap", values_from = "r") %>%
  filter(abs(UMAP1) > 0.4 | abs(UMAP2) > 0.4)

correlations <- correlations[rownames(correlations) %in% kin_umap_cor$kinase, ]
pvalues <- pvalues[rownames(correlations), ]

heatmap <- Heatmap(
  matrix = correlations,
  name = "Pearson's r",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(pvalues[i, j] < 0.05)
      grid.text("*", x, y, gp = gpar(fontsize = 10))
  })

#+ fig.width=8, fig.height=12
heatmap

pdf(file = "./output/plots/immune_infiltration/kin_act_immune_absolute_cor_most_variable_kinases.pdf", width = 8, height = 12)
heatmap
dev.off()
