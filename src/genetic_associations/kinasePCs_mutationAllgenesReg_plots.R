#' ---
#' title: "Kinase principal components to mutation association - plots"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=12, fig.height=6)


#' load R packages
library(RColorBrewer)
library(tidyverse)
library(viridis)
library(cowplot)
library(ggrepel)

source("./src/utils/cor_functions.R")
source("./src/utils/kinPCA_mut_plot.R")

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load samples metadata
metadata <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
  filter(data == "phosphorylation")


#' load kinase principal components
prin_comps <- data.table::fread("./data/Danish/kinasePCAMat.tsv") %>%
  as_tibble() %>%
  rename(sample = V1) %>%
  select(1:11) %>%
  pivot_longer(-sample, names_to = "pc", values_to = "value")


#' load PCA loadings
loadings <- data.table::fread("./data/Danish/kinasePCALoadings.tsv") %>%
  as_tibble() %>%
  rename(gene = V1) %>%
  select(1:11) %>%
  pivot_longer(-gene, names_to = "pc", values_to = "value")


#' load kinase activities and imputed values
kin_matrix <- data.table::fread("./data/Danish/kinaseActMatImputed.tsv") %>%
  as_tibble() %>%
  rename(gene = V1) %>%
  pivot_longer(-gene, names_to = "sample", values_to = "activity")


#' correlate kinases (variables) to principal components
pc_kinaseCor <- kin_matrix %>%
  inner_join(prin_comps, by = "sample") %>%
  select(pc, gene, sample, value, activity) %>%
  group_by(pc, gene) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = corr, x = "value", y = "activity")) %>%
  select(-data) %>%
  unnest()


#' load mutation matrix
mut <- read_tsv(file = "./output/files/mutations_matrix.txt.gz")

# filter mutation matrix
mut <- mut %>%
  filter(n >= 50) %>%
  select(-n) %>%
  pivot_longer(cols = -c(gene), names_to = "sample", values_to = "mutated")


#' load all mutation to principal component associations
pc_mut_associations <- read_tsv("./output/files/kaPCs_mutStatus_allGenes.txt.gz") %>%
  filter(p.adjust < 0.05)


#' load all mutation to kinase activity associations
kinase_mut_associations <- read_tsv("./output/files/ka_mutStatus_allGenes.txt.gz") %>%
  filter(p.adjust < 0.05)


#' get all combinations of PCs
# 1/2 ways
combinations <- pc_mut_associations %>%
  select(x=prin_comp) %>%
  distinct() %>%
  mutate(y=x) %>%
  expand(x, y) %>%
  filter(x != y) %>%
  mutate(id = map2_chr(.x = x, .y = y, ~ str_c(sort(c(.x,.y)), collapse = "_"))) %>%
  distinct(id, .keep_all = T) %>%
  select(-id)

# 2/2 ways
conct <- function(x, y) str_c(x, 1:y, sep = "_")
combinations <- pc_mut_associations %>%
  select(x=prin_comp) %>%
  distinct() %>%
  mutate(y=x) %>%
  mutate_if(is.character, .funs = conct, y = nrow(.)) %>%
  expand(x, y) %>%
  filter(x != y) %>%
  separate(col = x, into = c("x", "id_x"), sep = "_") %>%
  separate(col = y, into = c("y", "id_y"), sep = "_") %>%
  select(x, y, id_x, id_y) %>%
  mutate_at(3:4, as.numeric) %>%
  mutate(id = map2_chr(.x = id_x, .y = id_y, ~ str_c(sort(c(.x,.y)), collapse = "_"))) %>%
  distinct(id, .keep_all = T) %>%
  select(-id, -id_x, -id_y)


#' plot kinase principal components to mutation association
plot <- pc_mut_associations %>%
  arrange((estimate)) %>%
  mutate(gene = fct_inorder(gene)) %>%
  mutate(prin_comp = fct_inorder(prin_comp)) %>%
  ggplot() +
  geom_point(mapping = aes(x = prin_comp, y = gene, color = estimate, size = p.adjust)) +
  scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
  theme(
    axis.text = element_text(color = "black", size = 8),
    axis.title = element_text(color = "black", size = 10),
    legend.text = element_text(color = "black", size = 6),
    legend.title = element_text(color = "black", size = 8),
    legend.key.width=unit(0.2,"cm"),
    legend.key.height=unit(0.2,"cm")) +
  #guides(color = guide_legend(override.aes = list(size=2)), size = guide_legend(override.aes = list(size=2))) +
  labs(x = "Kinase principal components", y = "Mutated genes", color = "Effect size", size = "Adjusted\nP-value")

#+ fig.width=3, fig.height=3
plot

ggsave(filename = "kinPCs_mutation_associations1.png", plot = plot, path = "./output/plots/genetic_associations/", width = 3, height = 3)
unlink("kinPCs_mutation_associations1.png")


plot <- pc_mut_associations %>%
  arrange((estimate)) %>%
  select(prin_comp, gene) %>%
  mutate(pair = str_c(str_c("Kin", prin_comp, sep = ""), gene, sep = " ~ ")) %>%
  inner_join(mut, by = "gene") %>%
  inner_join(prin_comps, by = c("prin_comp" = "pc", "sample")) %>%
  mutate(pair = fct_inorder(pair)) %>%
  ggplot(mapping = aes(x = as.character(mutated), y = value)) +
  geom_boxplot(mapping = aes(fill = mutated), show.legend = F, alpha = 0.5, outlier.shape = NA, lwd = 0.2) +
  geom_jitter(mapping = aes(color = mutated), show.legend = F, alpha = 0.5, width = 0.1, size = 0.1) +
  #facet_wrap(~ pair, nrow = 1, scales = "free") +
  facet_wrap(facets = vars(prin_comp, gene), nrow = 1, scales = "free") +
  theme(
    axis.text = element_text(color = "black", size = 5),
    axis.title = element_text(color = "black", size = 7),
    axis.ticks = element_line(color = "black", size = 0.2),
    strip.text = element_text(color = "black", size = 5)) +
  labs(x = "Mutated", y = "Principal component")

#+ fig.width=5, fig.height=2
plot

ggsave(filename = "kinPCs_mutation_associations2.png", plot = plot, path = "./output/plots/genetic_associations/", width = 5, height = 2)
unlink("kinPCs_mutation_associations2.png")


#' PCA space scatterplots
#plots <- tibble(
#  x = c("PC1", "PC1", "PC1", "PC1", "PC1", "PC1", "PC1", "PC1", "PC1", "PC2", "PC3"),
#  y = c("PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC6", "PC7", "PC10", "PC3", "PC4"),
#  gx = c("VHL", "VHL", "VHL", "VHL", "VHL", "VHL", "VHL", "VHL", "VHL", "EGFR", "USH2A"),
#  gY = c("EGFR", "USH2A", "PTEN", "KRAS", "EGFR", "SALL1", "VWF", "KRAS", "PTEN", "USH2A", "PTEN"))

plots <- combinations %>%
  inner_join(pc_mut_associations[, c("prin_comp", "gene")], by = c("x" = "prin_comp")) %>%
  inner_join(pc_mut_associations[, c("prin_comp", "gene")], by = c("y" = "prin_comp")) %>%
  rename(gx = gene.x, gy = gene.y)

plotMapper <- function(x, y, gx, gy){
  
  scatterplot <- kinPCA_mut_plot(prin_comps, loadings, mut, metadata, kinase_mut_associations, x, y, gx, gy, 11)
  
  file_name = str_c("kin_pca_space_", tolower(x), "_", gx, "_", tolower(y), "_", gy, ".png")
  
  ggsave(filename = file_name, plot = scatterplot, path = "./output/plots/genetic_associations/", width = 12, height = 6)
  unlink(file_name)
  
  scatterplot
}

plots <- pmap(.l = plots, .f = ~ plotMapper(x = ..1, y = ..2, gx = ..3, gy = ..4))



plots[[1]]

plots[[2]]

plots[[3]]

plots[[4]]

plots[[5]]

plots[[6]]

plots[[7]]

plots[[8]]

plots[[9]]






#' PC3 to PC5 scatterplot\
#' highlight samples with mutations on KRAS (PC5 associated) and EGFR (PC5 associated)

lim = 11

# first scatterplot
df1 <- prin_comps %>%
  filter(pc %in% c("PC3", "PC5")) %>%
  pivot_wider(names_from = "pc", values_from = "value") %>%
  inner_join(mut[mut$gene == "KRAS", ], by = "sample") %>%
  mutate(mutated = as.character(mutated)) %>%
  inner_join(metadata[, c("sample", "batch", "cancer")], by = "sample")

scatterplot1 <- df1 %>%
  ggplot() +
  theme_classic() +
  geom_point(mapping = aes(x = PC3, y = PC5, color = mutated), alpha = 0.7) +
  geom_hline(yintercept = 0, lwd = .5, color = "grey") +
  geom_vline(xintercept = 0, lwd = .5, color = "grey") +
  geom_segment(aes(x = 10, y = -.2, xend = 10, yend = .2), lwd = .5, color = "grey") +
  geom_segment(aes(x = 5, y = -.2, xend = 5, yend = .2), lwd = .5, color = "grey") +
  geom_segment(aes(x = -10, y = -.2, xend = -10, yend = .2), lwd = .5, color = "grey") +
  geom_segment(aes(x = -5, y = -.2, xend = -5, yend = .2), lwd = .5, color = "grey") +
  geom_segment(aes(x = -.2, y = 5, xend = .2, yend = 5), lwd = .5, color = "grey") +
  geom_segment(aes(x = -.2, y = 10, xend = .2, yend = 10), lwd = .5, color = "grey") +
  geom_segment(aes(x = -.2, y = -5, xend = .2, yend = -5), lwd = .5, color = "grey") +
  geom_segment(aes(x = -.2, y = -10, xend = .2, yend = -10), lwd = .5, color = "grey") +
  scale_y_continuous(limits = c(-lim, lim)) +
  scale_x_continuous(limits = c(-lim, lim)) +
  scale_color_manual(values = c("#67a9cf", "red"), labels = c("No", "Yes")) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    plot.title = element_text(hjust = 0.5)) +
  labs(x = "PC3", y = "PC5", title = str_c("KRAS", " (", "PC5", " associated)"),  color = str_c("KRAS",  " LoF mutation"))

mutGx <- kinase_mut_associations[ kinase_mut_associations$gene == "KRAS", "kinase", drop=T]

df1_loadings <- loadings %>%
  filter(pc %in% c("PC3", "PC5")) %>%
  pivot_wider(names_from = "pc", values_from = "value") %>%
  rename(kinase = gene) %>%
  filter(kinase %in% mutGx) %>%
  mutate(exp = map2(.x = PC3, .y = PC5, .f = expand_coord, limit = lim)) %>%
  unnest()
  
scatterplot1 <- scatterplot1 +
  geom_segment(data = df1_loadings, aes(x = 0, y = 0, xend = x, yend = y), arrow = arrow(length = unit(1/2, "picas"))) +
  geom_text_repel(data = df1_loadings, mapping = aes(x = x, y = y, label = kinase), box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines"), size = 4, colour="black")

barplot1 <- df1 %>%
  mutate(cancer = fct_infreq(cancer)) %>%
  ggplot() +
  theme_classic() +
  coord_flip() +
  geom_bar(mapping = aes(x = cancer, fill = mutated), position = "fill") +
  scale_fill_manual(values = c("#67a9cf", "red"), labels = c("No", "Yes"), guide = F) +
  theme(axis.text = element_text(color = "black"), axis.title = element_text(color = "black")) +
  labs(x = "Cancer", y = "%", fill = "Mutated")

plot1 <- ggdraw(scatterplot1) +
  draw_plot(barplot1, .66, .65, .33, .3)



# second scatterplot
df2 <- prin_comps %>%
  filter(pc %in% c("PC3", "PC5")) %>%
  pivot_wider(names_from = "pc", values_from = "value") %>%
  inner_join(mut[mut$gene == "EGFR", ], by = "sample") %>%
  mutate(mutated = as.character(mutated)) %>%
  inner_join(metadata[, c("sample", "batch", "cancer")], by = "sample")

scatterplot2 <- df2 %>%
  ggplot() +
  theme_classic() +
  geom_point(mapping = aes(x = PC3, y = PC5, color = mutated), alpha = 0.7) +
  geom_hline(yintercept = 0, lwd = .5, color = "grey") +
  geom_vline(xintercept = 0, lwd = .5, color = "grey") +
  geom_segment(aes(x = 10, y = -.2, xend = 10, yend = .2), lwd = .5, color = "grey") +
  geom_segment(aes(x = 5, y = -.2, xend = 5, yend = .2), lwd = .5, color = "grey") +
  geom_segment(aes(x = -10, y = -.2, xend = -10, yend = .2), lwd = .5, color = "grey") +
  geom_segment(aes(x = -5, y = -.2, xend = -5, yend = .2), lwd = .5, color = "grey") +
  geom_segment(aes(x = -.2, y = 5, xend = .2, yend = 5), lwd = .5, color = "grey") +
  geom_segment(aes(x = -.2, y = 10, xend = .2, yend = 10), lwd = .5, color = "grey") +
  geom_segment(aes(x = -.2, y = -5, xend = .2, yend = -5), lwd = .5, color = "grey") +
  geom_segment(aes(x = -.2, y = -10, xend = .2, yend = -10), lwd = .5, color = "grey") +
  scale_y_continuous(limits = c(-lim, lim)) +
  scale_x_continuous(limits = c(-lim, lim)) +
  scale_color_manual(values = c("#67a9cf", "red"), labels = c("No", "Yes")) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    plot.title = element_text(hjust = 0.5)) +
  labs(x = "PC3", y = "PC5", title = str_c("EGFR", " (", "PC5", " associated)"),  color = str_c("EGFR",  " LoF mutation"))

mutGy <- kinase_mut_associations[ kinase_mut_associations$gene == "EGFR", "kinase", drop=T]

df2_loadings <- loadings %>%
  filter(pc %in% c("PC3", "PC5")) %>%
  pivot_wider(names_from = "pc", values_from = "value") %>%
  rename(kinase = gene) %>%
  filter(kinase %in% mutGy) %>%
  mutate(exp = map2(.x = PC3, .y = PC5, .f = expand_coord, limit = lim)) %>%
  unnest()
  
scatterplot2 <- scatterplot2 +
  geom_segment(data = df2_loadings, aes(x = 0, y = 0, xend = x, yend = y), arrow = arrow(length = unit(1/2, "picas"))) +
  geom_text_repel(data = df2_loadings, mapping = aes(x = x, y = y, label = kinase), box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines"), size = 4, colour="black")


barplot2 <- df2 %>%
  mutate(cancer = fct_infreq(cancer)) %>%
  ggplot() +
  theme_classic() +
  coord_flip() +
  geom_bar(mapping = aes(x = cancer, fill = mutated), position = "fill") +
  scale_fill_manual(values = c("#67a9cf", "red"), labels = c("No", "Yes"), guide = F) +
  theme(axis.text = element_text(color = "black"), axis.title = element_text(color = "black")) +
  labs(x = "Cancer", y = "%", fill = "Mutated")

plot2 <- ggdraw(scatterplot2) +
  draw_plot(barplot2, .66, .65, .33, .3)

scatterplots <- plot_grid(plot1, plot2)

scatterplots

ggsave(filename = "kin_pca_space_pc3_pc5_KRAS_EGFR.png", plot = scatterplots, path = "./output/plots/genetic_associations/", width = 12, height = 6)
unlink("kin_pca_space_pc3_pc5_KRAS_EGFR.png")
