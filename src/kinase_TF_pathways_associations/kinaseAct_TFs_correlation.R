#' ---
#' title: "Kinase to transcription factors/pathways correlation"
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
metadata <- getSamples(samples_annotation, c("protein"))


#' load transcription factor activities
#tfs <- data.table::fread("./data/progeny/TF_activity_FPKM.csv") %>%
tfs <- data.table::fread("./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  as_tibble() %>%
  rename(tf=V1) %>%
  pivot_longer(-tf, names_to = "sample", values_to = "tf_activity")


#' load pathway activities
pathways <- data.table::fread("./data/progeny/pathway_prog_log2FC.csv") %>%
  as_tibble() %>%
  rename(sample=V1) %>%
  pivot_longer(-sample, names_to = "pathway", values_to = "pathway_activity") %>%
  mutate(sample = str_replace(sample, "^X{1}", ""))


#' load kinase-activity inference data without imputations\
#' (quantile-normalized protein regressed-out phosphorylation data)\
#' select quantifications with more than 3 substrates
k_subN <- 3
ka <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  select(sample, kinase, kin_activity=log10P)


#' load kinase activities and imputed values\
#' matrix used for PCA analysis and UMAP
kin_activities <- data.table::fread("./data/Danish/kinaseActMatImputed.tsv") %>%
  as_tibble() %>%
  rename(kinase = V1) %>%
  pivot_longer(-kinase, names_to = "sample", values_to = "kin_activity")


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


#' # correlate UMAP projections with transcription factor activities
cor_mat <- umap %>%
  inner_join(tfs, by = "sample") %>%
  group_by(umap, tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x=data, .f = ~ broom::tidy(cor.test(.x$value, .x$tf_activity)) %>% select(estimate, p.value))) %>%
  select(-data) %>%
  unnest()

correlations <- cor_mat %>%
  select(-p.value) %>%
  pivot_wider(names_from = "umap", values_from = "estimate") %>%
  as.data.frame() %>%
  column_to_rownames(var = "tf") %>%
  as.matrix() %>%
  t()

pvalues <- cor_mat %>%
  select(-estimate) %>%
  pivot_wider(names_from = "umap", values_from = "p.value") %>%
  as.data.frame() %>%
  column_to_rownames(var = "tf") %>%
  as.matrix() %>%
  t()

heatmap <- Heatmap(
  matrix = correlations,
  name = "Pearson's r",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 4),
  row_title = "UMAP projections",
  column_title = "Transcription factors",
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(pvalues[i, j] < 0.05)
      grid.text("*", x, y, gp = gpar(fontsize = 8))
})

#+ fig.width=16, fig.height=2
heatmap

pdf(file = "./output/plots/kinase_TF_pathways_associations/umap_tf_activity_correlation_heatmap.pdf", width = 16, height = 2)
heatmap
dev.off()


#' # correlate UMAP projections with pathway activities
cor_mat <- umap %>%
  inner_join(pathways, by = "sample") %>%
  group_by(umap, pathway) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x=data, .f = ~ broom::tidy(cor.test(.x$value, .x$pathway_activity)) %>% select(estimate, p.value))) %>%
  select(-data) %>%
  unnest()

correlations <- cor_mat %>%
  select(-p.value) %>%
  pivot_wider(names_from = "umap", values_from = "estimate") %>%
  as.data.frame() %>%
  column_to_rownames(var = "pathway") %>%
  as.matrix() %>%
  t()

pvalues <- cor_mat %>%
  select(-estimate) %>%
  pivot_wider(names_from = "umap", values_from = "p.value") %>%
  as.data.frame() %>%
  column_to_rownames(var = "pathway") %>%
  as.matrix() %>%
  t()

heatmap <- Heatmap(
  matrix = correlations,
  name = "Pearson's r",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  row_title = "UMAP projections",
  column_title = "Pathways",
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(pvalues[i, j] < 0.05)
      grid.text("*", x, y, gp = gpar(fontsize = 8))
  })

#+ fig.width=4, fig.height=2
heatmap

pdf(file = "./output/plots/kinase_TF_pathways_associations/umap_pathway_activity_correlation_heatmap.pdf", width = 4, height = 2)
heatmap
dev.off()


#' # correlate principal components with transcription factor activities
cor_mat <- pcs %>%
  inner_join(tfs, by = "sample") %>%
  group_by(pc, tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x=data, .f = ~ broom::tidy(cor.test(.x$value, .x$tf_activity)) %>% select(estimate, p.value))) %>%
  select(-data) %>%
  unnest()

correlations <- cor_mat %>%
  select(-p.value) %>%
  pivot_wider(names_from = "pc", values_from = "estimate") %>%
  as.data.frame() %>%
  column_to_rownames(var = "tf") %>%
  as.matrix() %>%
  t()

pvalues <- cor_mat %>%
  select(-estimate) %>%
  pivot_wider(names_from = "pc", values_from = "p.value") %>%
  as.data.frame() %>%
  column_to_rownames(var = "tf") %>%
  as.matrix() %>%
  t()

heatmap <- Heatmap(
  matrix = correlations,
  name = "Pearson's r",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 4),
  row_title = "Kinase principal components",
  column_title = "Transcription factors")

#+ fig.width=16, fig.height=4
heatmap

pdf(file = "./output/plots/kinase_TF_pathways_associations/pcs_tf_activity_correlation_heatmap.pdf", width = 16, height = 4)
heatmap
dev.off()


#' # correlate principal components with pathway activities
cor_mat <- pcs %>%
  inner_join(pathways, by = "sample") %>%
  group_by(pc, pathway) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x=data, .f = ~ broom::tidy(cor.test(.x$value, .x$pathway_activity)) %>% select(estimate, p.value))) %>%
  select(-data) %>%
  unnest()

correlations <- cor_mat %>%
  select(-p.value) %>%
  pivot_wider(names_from = "pc", values_from = "estimate") %>%
  as.data.frame() %>%
  column_to_rownames(var = "pathway") %>%
  as.matrix() %>%
  t()

pvalues <- cor_mat %>%
  select(-estimate) %>%
  pivot_wider(names_from = "pc", values_from = "p.value") %>%
  as.data.frame() %>%
  column_to_rownames(var = "pathway") %>%
  as.matrix() %>%
  t()

heatmap <- Heatmap(
  matrix = correlations,
  name = "Pearson's r",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  row_title = "Kinase principal components",
  column_title = "Pathways")

#+ fig.width=4, fig.height=4
heatmap

pdf(file = "./output/plots/kinase_TF_pathways_associations/pcs_pathway_activity_correlation_heatmap.pdf", width = 4, height = 4)
heatmap
dev.off()


#' # correlate kinase and transcription factor activities\
#' ## use the kinase activities without imputations\
#' select kinase-TF pairs with more than 10 samples
cor_mat <- ka %>%
  inner_join(tfs, by = "sample") %>%
  group_by(kinase, tf) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, ~ nrow(.x) > 10)) %>%
  mutate(cor = map(.x=data, .f = ~ broom::tidy(cor.test(.x$kin_activity, .x$tf_activity)) %>% select(estimate, p.value))) %>%
  select(-data) %>%
  unnest()
write_tsv(cor_mat, "./output/files/kinaseNotImputed_TF_activity_correlation.txt.gz")

correlations <- cor_mat %>%
  select(-p.value) %>%
  pivot_wider(names_from = "kinase", values_from = "estimate") %>%
  as.data.frame() %>%
  column_to_rownames(var = "tf") %>%
  as.matrix() %>%
  t()

pvalues <- cor_mat %>%
  select(-estimate) %>%
  pivot_wider(names_from = "kinase", values_from = "p.value") %>%
  as.data.frame() %>%
  column_to_rownames(var = "tf") %>%
  as.matrix() %>%
  t()

heatmap <- Heatmap(
  matrix = correlations,
  name = "Pearson's r",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 4),
  row_title = "Kinases",
  column_title = "Transcription factors")

#+ fig.width=17, fig.height=10
heatmap

pdf(file = "./output/plots/kinase_TF_pathways_associations/kin_tf_activity_correlation_heatmap_all_kinases.pdf", width = 17, height = 10)
heatmap
dev.off()


#' # correlate kinase and transcription factor activities\
#' ## use the kinase activities with imputations (matrix used for PCA and UMAP)
cor_mat <- kin_activities %>%
  inner_join(tfs, by = "sample") %>%
  group_by(kinase, tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x=data, .f = ~ broom::tidy(cor.test(.x$kin_activity, .x$tf_activity)) %>% select(estimate, p.value))) %>%
  select(-data) %>%
  unnest()
write_tsv(cor_mat, "./output/files/kinaseImputed_TF_activity_correlation.txt.gz")

correlations <- cor_mat %>%
  select(-p.value) %>%
  pivot_wider(names_from = "kinase", values_from = "estimate") %>%
  as.data.frame() %>%
  column_to_rownames(var = "tf") %>%
  as.matrix() %>%
  t()

pvalues <- cor_mat %>%
  select(-estimate) %>%
  pivot_wider(names_from = "kinase", values_from = "p.value") %>%
  as.data.frame() %>%
  column_to_rownames(var = "tf") %>%
  as.matrix() %>%
  t()

heatmap <- Heatmap(
  matrix = correlations,
  name = "Pearson's r",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 4),
  row_title = "Kinases",
  column_title = "Transcription factors")

#+ fig.width=17, fig.height=7
heatmap

pdf(file = "./output/plots/kinase_TF_pathways_associations/kin_tf_activity_correlation_heatmap_imputed_kinases.pdf", width = 17, height = 7)
heatmap
dev.off()


#' # correlate kinase and pathway activities\
#' ## use the kinase activities with imputations (matrix used for PCA and UMAP)
cor_mat <- kin_activities %>%
  inner_join(pathways, by = "sample") %>%
  group_by(kinase, pathway) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x=data, .f = ~ broom::tidy(cor.test(.x$kin_activity, .x$pathway_activity)) %>% select(estimate, p.value))) %>%
  select(-data) %>%
  unnest()

correlations <- cor_mat %>%
  select(-p.value) %>%
  pivot_wider(names_from = "kinase", values_from = "estimate") %>%
  as.data.frame() %>%
  column_to_rownames(var = "pathway") %>%
  as.matrix() %>%
  t()

pvalues <- cor_mat %>%
  select(-estimate) %>%
  pivot_wider(names_from = "kinase", values_from = "p.value") %>%
  as.data.frame() %>%
  column_to_rownames(var = "pathway") %>%
  as.matrix() %>%
  t()

heatmap <- Heatmap(
  matrix = correlations,
  name = "Pearson's r",
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 6),
  row_title = "Kinases",
  column_title = "Pathways")

#+ fig.width=3, fig.height=6
heatmap

pdf(file = "./output/plots/kinase_TF_pathways_associations/kin_pathway_activity_correlation_heatmap_imputed_kinases.pdf", width = 3, height = 6)
heatmap
dev.off()

