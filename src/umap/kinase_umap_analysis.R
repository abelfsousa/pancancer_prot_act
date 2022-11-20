#' ---
#' title: "UMAP analysis on the kinase activities"
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


#' load kinase activities and imputed values
kin_matrix <- data.table::fread("./data/Danish/kinaseActMatImputed.tsv") %>%
  as.data.frame() %>%
  column_to_rownames(var = "V1") %>%
  t()

#kin_matrix <- read_tsv("./output/files/kinaseActMatImputed_batch_regOut.txt") %>%
#  column_to_rownames(var = "kinase") %>%
#  t()

kin_matrix_var <- kin_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(-sample, names_to = "kinase", values_to = "activity") %>%
  group_by(kinase) %>%
  summarise(act_sd = sd(activity)) %>%
  ungroup() %>%
  arrange(desc(act_sd))


#' load samples annotation
source("./src/utils/getSamples.R")
samples_annotation <- read_tsv(file = "./output/files/all_samples_annotation.txt")
samples_annotation <- getSamples(samples_annotation, c("phosphorylation")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue) %>%
  mutate(batch = str_replace(batch, "tcga-brca", "CPTAC-breast\n(TCGA)")) %>%
  mutate(batch = str_replace(batch, "tcga-ov", "CPTAC-ovarian\n(TCGA)")) %>%
  mutate(batch = str_replace(batch, "cbttc", "CPTAC-brain")) %>%
  mutate(batch = str_replace(batch, "discovery-ccrcc", "CPTAC-kidney\n(discovery)")) %>%
  mutate(batch = str_replace(batch, "discovery-luad", "CPTAC-lung\n(discovery)")) %>%
  mutate(batch = str_replace(batch, "discovery-ucec", "CPTAC-uterus\n(discovery)")) %>%
  mutate(batch = str_replace(batch, "colon-opportunities", "CPTAC-colon\n(confirmatory)")) %>%
  mutate(batch = str_replace(batch, "eogc-proteogenomics", "CPTAC-stomach")) %>%
  mutate(batch = str_replace(batch, "hcc-proteogenomics", "CPTAC-liver")) %>%
  mutate(batch = str_replace(batch, "ccle-colorectal", "CCLE-colorectal"))


#' ## create UMAP projection (kinases as variables)
set.seed(123)
kin_umap1 = umap(kin_matrix)

projections <- kin_umap1 %>%
  pluck("layout") %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  as_tibble() %>%
  rename(UMAP1 = V1, UMAP2 = V2) %>%
  arrange(sample)
write_tsv(projections, "./output/files/kinActivities_umap_projections.txt")


#' ## correlate UMAP projections with the batch covariate
batch <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
  getSamples("phosphorylation") %>%
  mutate(batch = map2_chr(batch, tissue, ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(sample, batch) %>%
  mutate(batch = as.factor(batch))

batch_matrix <- model.matrix(object = as.formula(~ batch + 0), data = batch, contrasts.arg=list(batch=contrasts(batch$batch, contrasts=F))) %>%
  as.data.frame() %>%
  as_tibble() %>%
  rename_all(.funs = ~ str_replace(.x, "batch", ""))

batch <- batch %>%
  bind_cols(batch_matrix) %>%
  select(-batch)

projections_batch <- projections %>%
  inner_join(batch, by = "sample") %>%
  select(-sample, -cbttc)

cor_mat <- t(cor(projections_batch[, 1:2], projections_batch[, 3:ncol(projections_batch)]))

heatmap <- Heatmap(
  column_title = "UMAP",
  column_title_gp = gpar(fontsize = 14),
  show_heatmap_legend = T,
  border = T,
  column_labels = c("1", "2"),
  column_names_rot = 0,
  rect_gp = gpar(col = "white", lwd = 0.5),
  column_names_centered = T,
  matrix = cor_mat,
  name = "Pearson's r",
  heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter", labels_gp = gpar(fontsize = 12), title_gp = gpar(fontsize = 14)),
  cluster_columns = F,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 12))

#+ fig.width=3, fig.height=4
draw(heatmap, heatmap_legend_side = "bottom")

pdf(file = "./output/plots/umap/kinase_umap_batch_cor.pdf", height = 4, width = 3)
draw(heatmap, heatmap_legend_side = "bottom")
dev.off()

png(file="./output/plots/umap/kinase_umap_batch_cor.png", units = "in", res = 300, height=4, width=3)
draw(heatmap, heatmap_legend_side = "bottom")
dev.off()


#' ## plot UMAP1 vs UMAP2 colored by batch
umap_plot <- projections %>%
  inner_join(samples_annotation, by = "sample") %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, color = batch)) +
  geom_point(size = 0.5, alpha = 1) +
  scale_color_brewer(type = "qual", palette = "Paired") +
  theme_classic() +
  theme(
    plot.title = element_text(colour = "black", size = 16, hjust = 0.5),
    axis.title = element_text(colour = "black", size = 16),
    axis.text = element_text(colour = "black", size = 14),
    #legend.title = element_text(colour = "black", size = 12, hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", size = 11, margin = margin(t = 4)),
    legend.position = "right") +
  labs(x = "UMAP1", y = "UMAP2", title = "Kinase activity UMAP", color = "Study") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))

#+ fig.width=4, fig.height=6
umap_plot

ggsave(filename = "kinase_umap_sample_projection.png", plot = umap_plot, path = "./output/plots/umap/", width = 4, height = 6)
ggsave(filename = "kinase_umap_sample_projection.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 4, height = 6)


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
  labs(x = "UMAP1", y = "UMAP2", title = "Kinase activity UMAP", color = "Study")

#+ fig.width=4, fig.height=6
umap_plot

ggsave(filename = "kinase_umap_sample_projection.png", plot = umap_plot, path = "./output/plots/umap/", width = 6, height = 6)
ggsave(filename = "kinase_umap_sample_projection.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 6, height = 6)


#' ## correlate each kinase with the UMAP projections
kinases <- kin_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  as_tibble() %>%
  arrange(sample)

cor_mat <- cor(kinases[, -c(1)], projections[, -c(1)])


#' ## heatmap of the correlations
heatmap <- Heatmap(
  matrix = cor_mat,
  name = "Pearson's r",
  cluster_columns = F,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8))

#+ fig.width=3, fig.height=10
heatmap

pdf(file="./output/plots/umap/umap_kinases_cor.pdf", height=10, width=3)
print(heatmap)
dev.off()

png(file="./output/plots/umap/umap_kinases_cor.png", units = "in", res = 300, height=10, width=3)
print(heatmap)
dev.off()


#' ## plot UMAP1 vs UMAP2 colored by kinases activity
sel <- cor_mat %>%
  abs() %>%
  as.data.frame() %>%
  rownames_to_column(var = "kinase") %>%
  as_tibble() %>%
  filter(UMAP1 > 0.55 | UMAP2 > 0.55)

umap_plot <- projections %>%
  inner_join(kinases, by = "sample") %>%
  #select(sample, UMAP1, UMAP2, CSNK2A1, AKT1) %>%
  select_at(.vars = c("sample", "UMAP1", "UMAP2", sel$kinase)) %>%
  pivot_longer(cols = -c(sample, UMAP1, UMAP2), names_to = "kinase", values_to = "activity") %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, color = activity)) +
  geom_point(size = 0.2) +
  facet_wrap(~ kinase) +
  scale_colour_gradient2(low = "#e31a1c", mid = "#d9d9d9", high = "darkblue") +
  theme_minimal()

#+ fig.width=6, fig.height=4
umap_plot

ggsave(filename = "umap_kinases_activity.png", plot = umap_plot, path = "./output/plots/umap/", width = 6, height = 4)
ggsave(filename = "umap_kinases_activity.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 6, height = 4)


#' ## remove kinase redundancies
source(file = "./data/Danish/filterRedundantTerms.R")
ks <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz") %>%
  select(-pair, -source) %>%
  distinct() %>%
  mutate(source_type = str_replace(source_type, "database|text-mining", "DB_text-mining")) %>%
  distinct() %>%
  filter(source_type == "DB_text-mining") %>%
  select(kinase, substrate) %>%
  distinct() %>%
  group_by(kinase) %>%
  summarise(substrate = list(substrate), n = n()) %>%
  ungroup() %>%
  filter(n >= 3) %>%
  select(-n) %>%
  unnest() %>%
  rename(ont=kinase, gene=substrate) %>%
  as.matrix()

# select non-redundant kinases
nr_kinases <- pruneGOTerms(ks, category = "broad", numTerms = 1, treeHeight = 0.85)

sel <- cor_mat %>%
  abs() %>%
  as.data.frame() %>%
  rownames_to_column(var = "kinase") %>%
  as_tibble() %>%
  filter(kinase %in% nr_kinases) %>%
  filter(UMAP1 > 0.45 | UMAP2 > 0.45)

most_var_kinases <- kin_matrix_var %>%
  filter(act_sd > quantile(act_sd, 0.5)) %>%
  pull(kinase)

sel <- cor_mat %>%
  #abs() %>%
  as.data.frame() %>%
  rownames_to_column(var = "kinase") %>%
  as_tibble() %>%
  filter(kinase %in% nr_kinases) %>%
  #filter(kinase %in% slice_max(kin_matrix_var, order_by = act_sd, n=40)$kinase)
  filter(kinase %in% most_var_kinases)


cor_mat2 <- cor_mat[rownames(cor_mat) %in% sel$kinase, ]

col_fun = circlize::colorRamp2(c(min(cor_mat2), (min(cor_mat2)+max(cor_mat2))/2, max(cor_mat2)), c("#d8b365", "#f5f5f5", "#5ab4ac"))
lgd = Legend(col_fun = col_fun, title = "Pearson'r", at = seq(-1,1,0.5), direction = "horizontal", title_position = "topcenter", legend_width = unit(3, "cm"))

heatmap <- Heatmap(
  col = col_fun,
  show_heatmap_legend = F,
  border = T,
  column_labels = c("Coord.1", "Coord.2"),
  column_names_rot = 22.5,
  rect_gp = gpar(col = "white", lwd = 0.5),
  column_names_centered = T,
  matrix = cor_mat2,
  name = "Pearson's r",
  cluster_columns = F,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8))

heatmap <- Heatmap(
  column_title = "UMAP",
  column_title_gp = gpar(fontsize = 14),
  #col = c("#d8b365", "#f5f5f5", "#5ab4ac"),
  col = c("#ca0020", "white", "#0571b0"),
  #clustering_method_rows = "average",
  show_heatmap_legend = T,
  border = T,
  column_labels = c("1", "2"),
  column_names_rot = 0,
  rect_gp = gpar(col = "white", lwd = 0.6),
  column_names_centered = T,
  matrix = cor_mat2,
  name = "Pearson's r",
  heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter", labels_gp = gpar(fontsize = 14), title_gp = gpar(fontsize = 16), at = c(-1, -0.5, 0, 0.5, 1), legend_width = unit(4, "cm")),
  cluster_columns = F,
  row_names_gp = gpar(fontsize = 11),
  column_names_gp = gpar(fontsize = 14))

#+ fig.width=2, fig.height=6
draw(heatmap, heatmap_legend_side = "bottom")

pdf(file = "./output/plots/umap/umap_kinases_cor2.pdf", height = 6, width = 2)
draw(heatmap, heatmap_legend_side = "bottom")
dev.off()

png(file="./output/plots/umap/umap_kinases_cor2.png", units = "in", res = 300, height = 6, width=2)
draw(heatmap, heatmap_legend_side = "bottom")
dev.off()


enr_kinases <- clusterProfiler::bitr(sel$kinase, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)[,2]
enr_universe <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz") %>%
  filter(source_type %in% c("database", "text-mining")) %>%
  select(kinase, substrate) %>%
  distinct() %>%
  group_by(kinase) %>%
  tally() %>%
  #filter(n >= 3) %>%
  pull(kinase)
enr_universe <- clusterProfiler::bitr(enr_universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db::org.Hs.eg.db)[,2]
  
enr_object <- clusterProfiler::enrichKEGG(gene = enr_kinases, organism = "hsa", keyType = "ncbi-geneid", use_internal_data = T, universe = enr_universe)

enr_table <- enr_object@result %>%
  as_tibble() %>%
  filter(p.adjust < 0.05)

categories <- enr_table %>%
  filter(str_detect(Description, c("MAPK|ErbB|cAMP|Ras|PI3K|mTOR|VEGF|Wnt"))) %>%
  mutate(log10_p = -log10(p.adjust)) %>%
  mutate(Description = fct_reorder(.f=Description, .x=log10_p, .fun = function(x) x))

categories_barplot <- ggplot(data = categories, mapping = aes(x = Description, y = log10_p, fill = Count)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  scale_fill_gradient( high = "#132B43", low = "#56B1F7") +
  theme(
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_blank()) +
  labs(x = "KEGG pathway", y = "Adjusted P-value (-log10)")
categories_barplot

ggsave(filename = "~/Desktop/categories_barplot.pdf", plot = categories_barplot, height = 2, width = 6)


umap_plot <- projections %>%
  inner_join(kinases, by = "sample") %>%
  #select_at(.vars = c("sample", "UMAP1", "UMAP2", sel$kinase)) %>%
  select_at(.vars = c("sample", "UMAP1", "UMAP2", "CSNK2A1", "CDK1", "PRKCB", "BRAF", "MAP2K1", "MAPK1")) %>%
  pivot_longer(cols = -c(sample, UMAP1, UMAP2), names_to = "kinase", values_to = "activity") %>%
  mutate(kinase = fct_relevel(kinase, "CSNK2A1", "CDK1", "PRKCB", "BRAF", "MAP2K1", "MAPK1")) %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, color = activity)) +
  geom_point(size = 1.5) +
  facet_wrap(~ kinase, nrow = 2) +
  scale_colour_gradient2(low = "#ca0020", mid = "white", high = "#0571b0", na.value = "white", name = "Kinase activity", breaks = seq(-2,2,1), limits = c(-2,2)) +
  #scale_colour_gradient(low = "#de2d26", high = "#2c7fb8") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 16),
    axis.text = element_text(colour = "black", size = 14),
    axis.title = element_text(colour = "black", size = 16),
    legend.position = "bottom",
    legend.key.width = unit(0.7, "cm"),
    legend.title = element_text(colour = "black", size = 18),
    legend.text = element_text(colour = "black", size = 16))

#+ fig.width=7, fig.height=6
umap_plot

ggsave(filename = "umap_kinases_activity2.png", plot = umap_plot, path = "./output/plots/umap/", width = 6, height = 6)
ggsave(filename = "umap_kinases_activity2.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 6, height = 6)


#' load all mutation to principal component associations
pc_mut_associations <- read_tsv("./output/files/kaPCs_mutStatus_allGenes.txt.gz") %>%
  filter(p.adjust < 0.05)


#' load mutation matrix\
#' filter mutation matrix by genes with mutations in at least 50 samples
mut <- read_tsv(file = "./output/files/mutations_matrix.txt.gz") %>%
  filter(n >= 50) %>%
  select(-n) %>%
  pivot_longer(cols = -c(gene), names_to = "sample", values_to = "mutated")


#' ## plot UMAP1 vs UMAP2 colored mutation state (PCs associations)
umap_plot <- projections %>%
  inner_join(mut, by = "sample") %>%
  filter(gene %in% unique(pc_mut_associations$gene)) %>%
  mutate(mutated = as.character(mutated)) %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, color = mutated)) +
  geom_point(size = 1) +
  facet_wrap(~ gene) +
  theme_minimal()

#+ fig.width=6, fig.height=3
umap_plot

ggsave(filename = "kinase_umap_pcs_mutation_association.png", plot = umap_plot, path = "./output/plots/umap/", width = 6, height = 3)
ggsave(filename = "kinase_umap_pcs_mutation_association.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 6, height = 3)
unlink("kinase_umap_pcs_mutation_association.png")
unlink("kinase_umap_pcs_mutation_association.pdf")


#' ## perform associations between the UMAP projections and the mutational state of samples\

#' prepare data to model
dataMod <- projections %>%
  pivot_longer(-sample, names_to = "UMAP", values_to = "value") %>%
  inner_join(mut, by = c("sample")) %>%
  select(UMAP, gene, sample, value, mutated) %>%
  rename(log10P = value) %>%
  group_by(UMAP, gene) %>%
  nest() %>%
  ungroup()

metadata <- samples_annotation %>%
  select(sample, batch)


#' run
associations1 <- dataMod %>%
  mutate(model = map(.x = data, .f = lr_test, covs = metadata, dep = "log10P", ind = c("batch", "mutated"), ind_drop = "mutated")) %>%
  select(-data) %>%
  unnest() %>%
  mutate(p.adjust = p.adjust(p.value, "BH"))

associations2 <- dataMod %>%
  mutate(model = map(.x = data, .f = linear_model, covs = metadata, mutN = 5)) %>%
  select(-data) %>%
  unnest() %>%
  mutate(p.adjust = p.adjust(p.value, "BH"))


#' ## plot UMAP1 vs UMAP2 colored by kinases activity\
#' ## highlight samples with mutations
umap_space <- projections %>%
  inner_join(kinases, by = "sample") %>%
  select_at(.vars = c("sample", "UMAP1", "UMAP2", sel$kinase)) %>%
  pivot_longer(cols = -c(sample, UMAP1, UMAP2), names_to = "kinase", values_to = "activity")

#' KRAS mutation
umap_plot <- umap_space %>%
  inner_join(mut[mut$gene == "KRAS", ], by = "sample") %>%
  mutate(mutated = as.character(mutated)) %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, fill = activity, color = mutated)) +
  geom_point(shape = 21, size = 1, stroke = 0.5) +
  facet_wrap(~ kinase) +
  scale_fill_gradient2(low = "#e31a1c", mid = "#d9d9d9", high = "darkblue") +
  scale_color_manual(values = c("#00000000", "#000000")) +
  theme_minimal() +
  labs(color = "KRAS\nmutation")

#+ fig.width=10, fig.height=6
umap_plot

ggsave(filename = "umap_kinases_activity_KRAS.png", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
ggsave(filename = "umap_kinases_activity_KRAS.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
unlink("umap_kinases_activity_KRAS.png")
unlink("umap_kinases_activity_KRAS.pdf")

#' EGFR mutation
umap_plot <- umap_space %>%
  inner_join(mut[mut$gene == "EGFR", ], by = "sample") %>%
  mutate(mutated = as.character(mutated)) %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, fill = activity, color = mutated)) +
  geom_point(shape = 21, size = 1, stroke = 0.5) +
  facet_wrap(~ kinase) +
  scale_fill_gradient2(low = "#e31a1c", mid = "#d9d9d9", high = "darkblue") +
  scale_color_manual(values = c("#00000000", "#000000")) +
  theme_minimal() +
  labs(color = "EGFR\nmutation")

#+ fig.width=10, fig.height=6
umap_plot

ggsave(filename = "umap_kinases_activity_EGFR.png", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
ggsave(filename = "umap_kinases_activity_EGFR.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
unlink("umap_kinases_activity_EGFR.png")
unlink("umap_kinases_activity_EGFR.pdf")

#' TP53 mutation
umap_plot <- umap_space %>%
  inner_join(mut[mut$gene == "TP53", ], by = "sample") %>%
  mutate(mutated = as.character(mutated)) %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, fill = activity, color = mutated)) +
  geom_point(shape = 21, size = 1, stroke = 0.5) +
  facet_wrap(~ kinase) +
  scale_fill_gradient2(low = "#e31a1c", mid = "#d9d9d9", high = "darkblue") +
  scale_color_manual(values = c("#00000000", "#000000")) +
  theme_minimal() +
  labs(color = "TP53\nmutation")

#+ fig.width=10, fig.height=6
umap_plot

ggsave(filename = "umap_kinases_activity_TP53.png", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
ggsave(filename = "umap_kinases_activity_TP53.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
unlink("umap_kinases_activity_TP53.png")
unlink("umap_kinases_activity_TP53.pdf")

#' VHL mutation
umap_plot <- umap_space %>%
  inner_join(mut[mut$gene == "VHL", ], by = "sample") %>%
  mutate(mutated = as.character(mutated)) %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, fill = activity, color = mutated)) +
  geom_point(shape = 21, size = 1, stroke = 0.5) +
  facet_wrap(~ kinase) +
  scale_fill_gradient2(low = "#e31a1c", mid = "#d9d9d9", high = "darkblue") +
  scale_color_manual(values = c("#00000000", "#000000")) +
  theme_minimal() +
  labs(color = "VHL\nmutation")

#+ fig.width=10, fig.height=6
umap_plot

ggsave(filename = "umap_kinases_activity_VHL.png", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
ggsave(filename = "umap_kinases_activity_VHL.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
unlink("umap_kinases_activity_VHL.png")
unlink("umap_kinases_activity_VHL.pdf")

#' ARID1A mutation
umap_plot <- umap_space %>%
  inner_join(mut[mut$gene == "ARID1A", ], by = "sample") %>%
  mutate(mutated = as.character(mutated)) %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, fill = activity, color = mutated)) +
  geom_point(shape = 21, size = 1, stroke = 0.5) +
  facet_wrap(~ kinase) +
  scale_fill_gradient2(low = "#e31a1c", mid = "#d9d9d9", high = "darkblue") +
  scale_color_manual(values = c("#00000000", "#000000")) +
  theme_minimal() +
  labs(color = "ARID1A\nmutation")

#+ fig.width=10, fig.height=6
umap_plot

ggsave(filename = "umap_kinases_activity_ARID1A.png", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
ggsave(filename = "umap_kinases_activity_ARID1A.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
unlink("umap_kinases_activity_ARID1A.png")
unlink("umap_kinases_activity_ARID1A.pdf")


#' MUC5B mutation
umap_plot <- umap_space %>%
  inner_join(mut[mut$gene == "MUC5B", ], by = "sample") %>%
  mutate(mutated = as.character(mutated)) %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, fill = activity, color = mutated)) +
  geom_point(shape = 21, size = 1, stroke = 0.5) +
  facet_wrap(~ kinase) +
  scale_fill_gradient2(low = "#e31a1c", mid = "#d9d9d9", high = "darkblue") +
  scale_color_manual(values = c("#00000000", "#000000")) +
  theme_minimal() +
  labs(color = "MUC5B\nmutation")

#+ fig.width=10, fig.height=6
umap_plot

ggsave(filename = "umap_kinases_activity_MUC5B.png", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
ggsave(filename = "umap_kinases_activity_MUC5B.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
unlink("umap_kinases_activity_MUC5B.png")
unlink("umap_kinases_activity_MUC5B.pdf")


#' PLEC mutation
umap_plot <- umap_space %>%
  inner_join(mut[mut$gene == "PLEC", ], by = "sample") %>%
  mutate(mutated = as.character(mutated)) %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, fill = activity, color = mutated)) +
  geom_point(shape = 21, size = 1, stroke = 0.5) +
  facet_wrap(~ kinase) +
  scale_fill_gradient2(low = "#e31a1c", mid = "#d9d9d9", high = "darkblue") +
  scale_color_manual(values = c("#00000000", "#000000")) +
  theme_minimal() +
  labs(color = "PLEC\nmutation")

#+ fig.width=10, fig.height=6
umap_plot

ggsave(filename = "umap_kinases_activity_PLEC.png", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
ggsave(filename = "umap_kinases_activity_PLEC.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
unlink("umap_kinases_activity_PLEC.png")
unlink("umap_kinases_activity_PLEC.pdf")


#' ZNF469 mutation
umap_plot <- umap_space %>%
  inner_join(mut[mut$gene == "ZNF469", ], by = "sample") %>%
  mutate(mutated = as.character(mutated)) %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, fill = activity, color = mutated)) +
  geom_point(shape = 21, size = 1, stroke = 0.5) +
  facet_wrap(~ kinase) +
  scale_fill_gradient2(low = "#e31a1c", mid = "#d9d9d9", high = "darkblue") +
  scale_color_manual(values = c("#00000000", "#000000")) +
  theme_minimal() +
  labs(color = "ZNF469\nmutation")

#+ fig.width=10, fig.height=6
umap_plot

ggsave(filename = "umap_kinases_activity_ZNF469.png", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
ggsave(filename = "umap_kinases_activity_ZNF469.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
unlink("umap_kinases_activity_ZNF469.png")
unlink("umap_kinases_activity_ZNF469.pdf")


#' CHD5 mutation
umap_plot <- umap_space %>%
  inner_join(mut[mut$gene == "CHD5", ], by = "sample") %>%
  mutate(mutated = as.character(mutated)) %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, fill = activity, color = mutated)) +
  geom_point(shape = 21, size = 1, stroke = 0.5) +
  facet_wrap(~ kinase) +
  scale_fill_gradient2(low = "#e31a1c", mid = "#d9d9d9", high = "darkblue") +
  scale_color_manual(values = c("#00000000", "#000000")) +
  theme_minimal() +
  labs(color = "CHD5\nmutation")

#+ fig.width=10, fig.height=6
umap_plot

ggsave(filename = "umap_kinases_activity_CHD5.png", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
ggsave(filename = "umap_kinases_activity_CHD5.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 10, height = 6)
unlink("umap_kinases_activity_CHD5.png")
unlink("umap_kinases_activity_CHD5.pdf")


#' ## plot UMAP1 vs UMAP2 colored by sample purity
purity <- read_tsv("./output/files/samples_purity.txt") %>%
  select(-metric)

umap_plot <- projections %>%
  inner_join(purity, by = "sample") %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, color = score)) +
  geom_point() +
  scale_colour_gradient2() +
  theme_minimal() +
  labs(color = "purity")

#+ fig.width=6, fig.height=3
umap_plot

ggsave(filename = "kinase_umap_samples_purity.png", plot = umap_plot, path = "./output/plots/umap/", width = 6, height = 3)
ggsave(filename = "kinase_umap_samples_purity.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 6, height = 3)
unlink("kinase_umap_samples_purity.png")
unlink("kinase_umap_samples_purity.pdf")


#' ## create UMAP projection (samples as variables)
set.seed(1234)
kin_umap2 = umap(t(kin_matrix))


#' ## plot UMAP1 vs UMAP2
umap_plot <- kin_umap2 %>%
  pluck("layout") %>%
  as.data.frame() %>%
  rownames_to_column(var = "kinase") %>%
  as_tibble() %>%
  rename(UMAP1 = V1, UMAP2 = V2) %>%
  ggplot(mapping = aes(x = UMAP1, y = UMAP2, label = kinase)) +
  geom_point(size = 1, color = "darkblue") +
  geom_text_repel(size = 2) +
  theme_minimal()

#+ fig.width=7, fig.height=4
umap_plot

ggsave(filename = "umap_kinase_projection.png", plot = umap_plot, path = "./output/plots/umap/", width = 7, height = 4)
ggsave(filename = "umap_kinase_projection.pdf", plot = umap_plot, path = "./output/plots/umap/", width = 7, height = 4)
unlink("umap_kinase_projection.png")
unlink("umap_kinase_projection.pdf")

