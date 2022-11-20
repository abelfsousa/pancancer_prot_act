library(tidyverse)
library(ComplexHeatmap)


# load kinase activities and imputed values
kin_matrix <- data.table::fread("./data/Danish/kinaseActMatImputed.tsv") %>%
  as.data.frame() %>%
  column_to_rownames(var = "V1") %>%
  t()


# perform PCA
kin_pca <- prcomp(kin_matrix, scale. = T, center = T)

correlations <- cor(kin_matrix, kin_pca$x[,1:20])

cor_hm <- Heatmap(
  correlations, 
  name = "Person's r")

pdf("./output/plots/cell_cycle/kinases_pcs_corr.pdf", width = 10, height = 20)
cor_hm
dev.off()


kin_matrix2 <- kin_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  as_tibble() %>%
  pivot_longer(cols = -sample, names_to = "kinase", values_to = "activity")

kin_pca2 <- kin_pca %>%
  pluck("x") %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  as_tibble() %>%
  pivot_longer(cols = -sample, names_to = "pc", values_to = "pc_value") %>%
  filter(pc %in% paste0(c("PC"), 1:20))

correlations2 <- kin_matrix2 %>%
  inner_join(kin_pca2, by = c("sample"))
