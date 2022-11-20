# load R packages
library(tidyverse)
library(ggrepel)

source("./src/utils/getSamples.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load samples annotation
samples_annotation <- read_tsv(file = "./output/files/all_samples_annotation.txt")
samples_annotation <- getSamples(samples_annotation, c("protein")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue)


# load immune cells estimates (absolute values)
cibersort <- read_csv(file = "./data/cibersort/CIBERSORT.CPTAC.absolute.csv") %>%
  rename(sample = X1) %>%
  #filter(str_detect(sample, "^X{1}"))
  mutate(sample = str_replace(sample, "^X{1}", "")) %>%
  select(-`P-value`, -Correlation, -RMSE, -`Absolute score (no.sumto1)`) %>%
  pivot_longer(-sample, names_to = "cell", values_to = "estimate")

#cibersort <- read_csv(file = "./data/cibersort/CIBERSORT.CPTAC.relative.csv") %>%
#  rename(sample = `Input Sample`) %>%
#  #filter(str_detect(sample, "^X{1}"))
#  mutate(sample = str_replace(sample, "^X{1}", "")) %>%
#  select(-`P-value`, -`Pearson Correlation`, -`RMSE`) %>%
#  pivot_longer(-sample, names_to = "cell", values_to = "estimate")



# load kinase activities and imputed values
# matrix used for PCA analysis and UMAP
kin_activities <- data.table::fread("./data/Danish/kinaseActMatImputed.tsv") %>%
  as_tibble() %>%
  rename(kinase = V1) %>%
  pivot_longer(-kinase, names_to = "sample", values_to = "kin_activity")


# load transcription factor activities
#tfs <- data.table::fread("./data/progeny/TF_activity_log2FC.csv") %>%
tfs <- data.table::fread("./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  as_tibble() %>%
  rename(tf=V1) %>%
  pivot_longer(-tf, names_to = "sample", values_to = "tf_activity")


# load string network (vertices names sorted for each edge; duplicates removed)
string_net <- read_tsv("./output/files/string_network_sorted_pairs.txt.gz")


# load kinase-TFs associations
sortPrt <- function(x, y){
  sorted <- sort(c(x,y))
  tibble(a=sorted[1],b=sorted[2])
}

kin_tf <- read_tsv(file = "./output/files/kinaseImputed_TF_activity_associations.txt.gz") %>%
  mutate(x = map2(kinase, tf, sortPrt)) %>%
  unnest() %>%
  select(a, b, everything()) %>%
  inner_join(string_net, by = c("a", "b")) %>%
  filter(padj < 0.05) %>%
  select(kinase, tf, kin_tf_beta = kin_beta, kin_tf_pval = kin_pval, kin_tf_padj = padj, string=weight)


# function to fit a multiple linear model
reg <- function(df, metadata, covars, coeff){
  dat <- df %>%
    inner_join(metadata, by = "sample")
  
  # fit the models
  f <- as.formula(str_c(covars[1], "~", str_c(covars[-1], collapse = "+")))
  mod <- lm(f, data = dat)
  
  # compute and extract multiple statistical values
  mod <- broom::tidy(mod)
  
  beta <- mod[mod$term == coeff, "estimate", drop=T]
  pval <- mod[mod$term == coeff, "p.value", drop=T]
  
  
  res <- tibble(beta=beta, pval=pval)
  
  res
}


# associate kinase and immune cells using linear models
# use the kinase activities with imputations (matrix used for PCA and UMAP)
# use experimental batch as covariate
associations <- kin_activities %>%
  inner_join(cibersort, by = "sample") %>%
  group_by(kinase, cell) %>%
  nest() %>%
  ungroup() %>%
  mutate(models = map(.x=data, .f = reg, metadata = samples_annotation, covars = c("estimate", "batch", "kin_activity"), coeff = "kin_activity")) %>%
  select(-data) %>%
  unnest() %>%
  rename(kin_cell_beta = beta, kin_cell_pval = pval) %>%
  mutate(kin_cell_padj = p.adjust(kin_cell_pval, method = "BH"))

write_tsv(associations, "./output/files/kinaseImputed_immune_cells_associations.txt.gz")

kin_cell <- associations %>%
  filter(kin_cell_padj < 0.05)

# volcano plot 
vplot <- associations %>%
  mutate(l = pmap_chr(.l = ., .f = ~ if(..5 < 0.05){str_c(..1, ..2, sep="--")}else{NA})) %>%
  mutate(c = map_chr(.x = kin_cell_padj, .f = ~ if(.x < 0.05){"1"}else{"2"})) %>%
  ggplot(mapping = aes(x = kin_cell_beta, y = -log10(kin_cell_padj), label = l, color = c)) +
  geom_point(size = 0.5, show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, color = "grey60") + 
  geom_text_repel(size = 1.5, segment.size = 0.1, color = "black") +
  scale_color_manual(values = c("red", "grey60")) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 12),
    legend.key.size = unit(0.4, "cm")) +
  labs(x = expression(beta ~ "(kinase activity)"), y = "Adjusted P-value (-log10)")

ggsave(filename = "kinase_immune_cells_associations_volcano_plot.png", plot = vplot, path = "./output/plots/immune_infiltration/", width = 4, height = 4)
ggsave(filename = "kinase_immune_cells_associations_volcano_plot.pdf", plot = vplot, path = "./output/plots/immune_infiltration/", width = 4, height = 4)
unlink("kinase_immune_cells_associations_volcano_plot.png")
unlink("kinase_immune_cells_associations_volcano_plot.pdf")


# associate TFs and immune cells using linear models
# use experimental batch as covariate
associations <- tfs %>%
  inner_join(cibersort, by = "sample") %>%
  group_by(tf, cell) %>%
  nest() %>%
  ungroup() %>%
  mutate(models = map(.x=data, .f = reg, metadata = samples_annotation, covars = c("estimate", "batch", "tf_activity"), coeff = "tf_activity")) %>%
  select(-data) %>%
  unnest() %>%
  rename(tf_cell_beta = beta, tf_cell_pval = pval) %>%
  mutate(tf_cell_padj = p.adjust(tf_cell_pval, method = "BH"))

write_tsv(associations, "./output/files/tf_immune_cells_associations.txt.gz")

tf_cell <- associations %>%
  filter(tf_cell_padj < 0.05)


# volcano plot 
vplot <- associations %>%
  mutate(l = pmap_chr(.l = ., .f = ~ if(-log10(..5) > 50){str_c(..1, ..2, sep="--")}else{NA})) %>%
  mutate(c = map_chr(.x = tf_cell_padj, .f = ~ if(.x < 0.05){"1"}else{"2"})) %>%
  ggplot(mapping = aes(x = tf_cell_beta, y = -log10(tf_cell_padj), label = l, color = c)) +
  geom_point(size = 0.5, show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, color = "grey60") + 
  geom_text_repel(size = 1.5, segment.size = 0.1, color = "black") +
  scale_color_manual(values = c("blue", "grey60")) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 12),
    legend.key.size = unit(0.4, "cm")) +
  labs(x = expression(beta ~ "(TF activity)"), y = "Adjusted P-value (-log10)")

ggsave(filename = "tf_immune_cells_associations_volcano_plot.png", plot = vplot, path = "./output/plots/immune_infiltration/", width = 4, height = 4)
ggsave(filename = "tf_immune_cells_associations_volcano_plot.pdf", plot = vplot, path = "./output/plots/immune_infiltration/", width = 4, height = 4)
unlink("tf_immune_cells_associations_volcano_plot.png")
unlink("tf_immune_cells_associations_volcano_plot.pdf")


# cases of immune cells associated with the same kinase and tf, and kinase and tf associated with each other
immune_associations <- kin_cell %>%
  inner_join(tf_cell, by = "cell") %>%
  select(cell, kinase, tf, everything()) %>%
  inner_join(kin_tf, by = c("kinase", "tf")) %>%
  arrange(desc(string))

write_tsv(immune_associations, "./output/files/kinase_tf_immune_cells_associations.txt")


