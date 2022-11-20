#' ---
#' title: "Kinase activity to mutation association - plots"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=6, fig.height=6)



#' load R packages
library(tidyverse)
library(viridis)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(ggbeeswarm)


source("./src/utils/KA_mut.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load mutation matrix
samp_mut = 20
mut <- read_tsv(file = "./output/files/mutations_matrix.txt.gz") %>%
  filter(n >= samp_mut) %>%
  select(-n) %>%
  pivot_longer(cols = -c(gene), names_to = "sample", values_to = "mutated")


#' load kinase-activity inference data
k_subN = 3
k_sampN = 10
ka <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  group_by(kinase) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > k_sampN) %>%
  select(-n)


#' load string network (vertices names sorted for each edge; duplicates removed)
string_net <- read_tsv("./output/files/string_network_sorted_pairs.txt.gz")


#' load mutation-kinase activity associations

# function to sort the kinase-gene pairs alphabetically
sortPrt <- function(x, y){
  sorted <- sort(c(x,y))
  tibble(a=sorted[1],b=sorted[2])
}

associations <- read_tsv(file = "./output/files/ka_mutStatus_allGenes.txt.gz") %>%
  filter(kinase != gene) %>%
  mutate(p.adjust = p.adjust(p.value, "BH")) %>%
  mutate(x = map2(kinase, gene, sortPrt)) %>%
  unnest() %>%
  select(a, b, everything()) %>%
  left_join(string_net, by = c("a", "b")) %>%
  mutate(signf = map_dbl(.x = p.adjust, .f = ~ if(.x < 0.01){3}else if(.x > 0.01 & .x < 0.05){2}else{1}))


#' load shortest paths for multiple string network cutoffs (scores)
#files <- c("200_100", "400_100", "600_100", "800_100")
#net_paths <- map(.x = files, .f = ~ read_rds(path = paste0("./output/r_objects/mut_associations_shortestPaths_", .x, ".rds")) %>% bind_rows() %>% mutate(run = .x))

#net_paths <- bind_rows(net_paths) %>%
#  filter(set == "non_random", run == "200_100") %>%
#  mutate(n_edges = as.numeric(n_edges)) %>%
#  filter(n_edges == 1)


#' volcano plot
volcano <- associations %>%
  mutate(p.adjust = -log10(p.adjust)) %>%
  mutate(l = pmap_chr(.l = ., .f = ~ if(..7 > 2){str_c(..4, ..3, sep = "--")}else{NA})) %>%
  ggplot(mapping = aes(x = estimate, y = p.adjust, color = as.character(signf), label = l)) +
  geom_point(size = 0.5) +
  geom_text_repel(size = 2, segment.size = 0.3, color = "black") +
  geom_vline(xintercept = 0, color = "grey", linetype = "longdash", size = 0.3) +
  #geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "longdash") +
  #geom_hline(yintercept = -log10(0.01), color = "grey", linetype = "longdash") +
  scale_colour_manual(values = c("#969696", "#006d2c", "#662506"), guide = F) +
  labs(x = "Effect size (beta)", y = "Adjusted P-value (-log10)")

#+ fig.width=6, fig.height=6
volcano

ggsave(filename = "kinase_association_pairs_volcano_plot.png", plot = volcano, path = "./output/plots/genetic_associations/", width = 6, height = 6)
unlink("kinase_association_pairs_volcano_plot.png")


#' volcano plot - colored using string network

# function that select the associations to show in the plot
helper <- function(p, w, k, g, W, P){
  if(is.na(w)){
    NA
  } else if(w >= W & p < P){
    str_c(g, k, sep = "--")
  } else {
    NA
  }
}

volcano <- associations %>%
  filter(gene != kinase) %>%
  mutate(l = pmap_chr(.l = ., .f = ~ helper(..7, ..8, ..3, ..4, W = 400, P = 0.05))) %>%
  mutate(p.adjust = -log10(p.adjust)) %>%
  ggplot(mapping = aes(x = estimate, y = p.adjust, color = weight, label = l)) +
  geom_point(size = 0.5) +
  geom_text_repel(size = 2, segment.size = 0.3, color = "black") +
  geom_vline(xintercept = 0, color = "grey", linetype = "longdash", size = 0.3) +
  #geom_hline(yintercept = -log10(0.05), color = "grey", linetype = "longdash") +
  #geom_hline(yintercept = -log10(0.01), color = "grey", linetype = "longdash") +
  labs(x = "Effect size (beta)", y = "Adjusted P-value (-log10)")

#+ fig.width=6, fig.height=6
volcano

ggsave(filename = "kinase_association_pairs_volcano_plot2.png", plot = volcano, path = "./output/plots/genetic_associations/", width = 6, height = 6)
unlink("kinase_association_pairs_volcano_plot2.png")


#' volcano plot - colored using discrete string edge weight scale

# function to discretize the string edge weight scale
helper <- function(w){
  if(is.na(w)){
    out <- "Unknown"
  } else if(w > 0 & w < 0.4){
    out <- "(0,0.4)"
  } else if(w >= 0.4 & w < 0.8){
    out <- "[0.4,0.8)"
  } else if(w >= 0.8){
    out <- "[0.8,1)"
  } else {
    out <- NA
  }
  out
}

volcano <- associations %>%
  filter(gene != kinase) %>%
  mutate(pair = pmap_chr(.l = ., .f = ~ if(is.na(..8)){NA}else if(..7 < 0.001 & ..8 >= 800){str_c(..4, ..3, sep = "--")}else{NA})) %>%
  mutate(weight = weight/1000) %>%
  mutate(string = map_chr(.x = weight, .f = helper)) %>%
  mutate(string = fct_relevel(string, "Unknown", "(0,0.4)", "[0.4,0.8)", "[0.8,1)")) %>%
  mutate(p.adjust = -log10(p.adjust)) %>%
  ggplot(mapping = aes(x = estimate, y = p.adjust, color = string, label = pair)) +
  geom_point(mapping = aes(alpha = string), size = 1) +
  geom_text_repel(size = 3.5, segment.size = 0.3, color = "black", point.padding = 0.5, box.padding = 0.5) +
  geom_vline(xintercept = 0, color = "black", linetype = "longdash", size = 0.1) +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "longdash", size = 0.1) +
  geom_hline(yintercept = -log10(0.01), color = "black", linetype = "longdash", size = 0.1) +
  theme_classic() +
  theme(
    plot.title = element_text(colour = "black", size = 18, hjust = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.title = element_text(colour = "black", size = 17),
    axis.text = element_text(colour = "black", size = 15),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 15),
    legend.title = element_text(colour = "black", size = 17)) +
  scale_alpha_manual(values = c(0.1, 0.8, 0.9, 1), guide = F) +
  scale_color_manual(values = c("#636363", "#6baed6", "#2171b5", "#08306b")) +
  scale_x_continuous(limits = c(-3,3), breaks = seq(-3,3,1)) +
  labs(x = expression("Mutation" ~ beta ~ "(effect size)"), y = "Adjusted P-value (-log10)", color = "String\nweight", title = "Kinases") +
  guides(color = guide_legend(override.aes = list(alpha=1, size = 2), nrow = 2))

#+ fig.width=4, fig.height=6
volcano

ggsave(filename = "kinase_association_pairs_volcano_plot3.png", plot = volcano, path = "./output/plots/genetic_associations/", width = 4, height = 6, dpi = 900)
ggsave(filename = "kinase_association_pairs_volcano_plot3.pdf", plot = volcano, path = "./output/plots/genetic_associations/", width = 4, height = 6)
unlink("kinase_association_pairs_volcano_plot3.png")
unlink("kinase_association_pairs_volcano_plot3.pdf")

volcano <- associations %>%
  filter(gene != kinase) %>%
  mutate(pair = pmap_chr(.l = ., .f = ~ if(is.na(..8)){NA}else if(..7 < 0.001 & ..8 >= 800){str_c(..4, ..3, sep = "--")}else{NA})) %>%
  mutate(weight = weight/1000) %>%
  mutate(string = map_chr(.x = weight, .f = helper)) %>%
  mutate(string = fct_relevel(string, "Unknown", "(0,0.4)", "[0.4,0.8)", "[0.8,1)")) %>%
  mutate(p.adjust = -log10(p.adjust)) %>%
  ggplot(mapping = aes(x = estimate, y = p.adjust, color = string, label = pair)) +
  geom_point(mapping = aes(alpha = string), size = 1) +
  geom_text_repel(size = 3.5, segment.size = 0.3, color = "black", point.padding = 0.5, box.padding = 0.5, seed = 123) +
  geom_vline(xintercept = 0, color = "black", linetype = "longdash", size = 0.1) +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "longdash", size = 0.1) +
  geom_hline(yintercept = -log10(0.01), color = "black", linetype = "longdash", size = 0.1) +
  theme_classic() +
  theme(
    plot.title = element_text(colour = "black", size = 18, hjust = 0.5),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.title = element_text(colour = "black", size = 17),
    axis.text = element_text(colour = "black", size = 15),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 15),
    legend.title = element_text(colour = "black", size = 17)) +
  scale_alpha_manual(values = c(0.1, 0.8, 0.9, 1), guide = F) +
  scale_color_manual(values = c("#636363", "#6baed6", "#2171b5", "#08306b")) +
  scale_x_continuous(limits = c(-3,3), breaks = seq(-3,3,1)) +
  labs(x = expression("Mutation" ~ beta ~ "(effect size)"), y = "Adjusted P-value (-log10)", color = "String network\nedge weight", title = "Kinases") +
  guides(color = guide_legend(override.aes = list(alpha=1, size = 2), nrow = 2))

#+ fig.width=5, fig.height=5
volcano

ggsave(filename = "kinase_association_pairs_volcano_plot4.png", plot = volcano, path = "./output/plots/genetic_associations/", width = 5, height = 5, dpi = 900)
ggsave(filename = "kinase_association_pairs_volcano_plot4.pdf", plot = volcano, path = "./output/plots/genetic_associations/", width = 5, height = 5)
unlink("kinase_association_pairs_volcano_plot4.png")
unlink("kinase_association_pairs_volcano_plot4.pdf")


volcano <- associations %>%
  filter(gene != kinase) %>%
  mutate(pair = pmap_chr(.l = ., .f = ~ if(is.na(..8)){NA}else if(..7 < 0.001 & ..8 >= 800){str_c(..4, ..3, sep = "--")}else{NA})) %>%
  mutate(weight = weight/1000) %>%
  mutate(string = map_chr(.x = weight, .f = helper)) %>%
  mutate(string = fct_relevel(string, "Unknown", "(0,0.4)", "[0.4,0.8)", "[0.8,1)")) %>%
  mutate(p.adjust = -log10(p.adjust)) %>%
  ggplot(mapping = aes(x = estimate, y = p.adjust, color = string, label = pair)) +
  geom_point(mapping = aes(alpha = string), size = 2) +
  geom_text_repel(size = 4, segment.size = 0.3, color = "black", point.padding = 0.5, box.padding = 0.5, seed = 123) +
  geom_vline(xintercept = 0, color = "black", linetype = "longdash", size = 0.1) +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "longdash", size = 0.1) +
  geom_hline(yintercept = -log10(0.01), color = "black", linetype = "longdash", size = 0.1) +
  theme_classic() +
  theme(
    plot.title = element_text(colour = "black", size = 18, hjust = 0.5),
    #panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.title = element_text(colour = "black", size = 17),
    axis.text = element_text(colour = "black", size = 15),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 15),
    legend.title = element_text(colour = "black", size = 17)) +
  scale_alpha_manual(values = c(0.1, 0.8, 0.9, 1), guide = F) +
  scale_color_manual(values = c("#636363", "#6baed6", "#2171b5", "#08306b")) +
  scale_x_continuous(limits = c(-3,3), breaks = seq(-3,3,1)) +
  labs(x = expression("Mutation" ~ beta ~ "(effect size)"), y = "Adjusted P-value (-log10)", color = "String network\nedge weight", title = "Kinases") +
  guides(color = guide_legend(override.aes = list(alpha=1, size = 3), nrow = 2))

#+ fig.width=5, fig.height=5
volcano

ggsave(filename = "kinase_association_pairs_volcano_plot5.png", plot = volcano, path = "./output/plots/genetic_associations/", width = 5, height = 5, dpi = 900)
ggsave(filename = "kinase_association_pairs_volcano_plot5.pdf", plot = volcano, path = "./output/plots/genetic_associations/", width = 5, height = 5)


#' plot kinases and mutated genes involved in the associations
mut_genes <- associations %>%
  filter(p.adjust < 0.05) %>%
  group_by(gene) %>%
  summarise(n = n(), mean_estimate = mean(estimate), sd_estimate = sd(estimate), mean_pval = mean(p.adjust)) %>%
  ungroup() %>%
  filter(n > 2) %>%
  mutate(set = "mutated genes")

kinases <- associations %>%
  filter(p.adjust < 0.05) %>%
  group_by(kinase) %>%
  summarise(n = n(), mean_estimate = mean(estimate), sd_estimate = sd(estimate), mean_pval = mean(p.adjust)) %>%
  ungroup() %>%
  filter(n > 2) %>%
  rename(gene=kinase) %>%
  mutate(set = "kinases")


plot_kinases <- kinases %>%
  mutate(gene = fct_reorder(.f = gene, .x = mean_estimate, .fun = function(x) x)) %>%
  ggplot() +
  geom_point(mapping = aes(x = gene, y = mean_estimate, color = mean_pval, size = n)) +
  geom_linerange(mapping = aes(x = gene, ymin = mean_estimate-(sd_estimate*2), ymax = mean_estimate+(sd_estimate*2), color = mean_pval), size=0.3) +
  scale_color_viridis() +
  scale_size(limits = c(3, 20)) +
  coord_flip() +
  facet_wrap(~ set, scales = "free") +
  labs(x = "", y = "Mean effect size (mean Beta)", color = "Mean Adj. P-value", size = "Number of\nassociations")

plot_genes <- mut_genes %>%
  mutate(gene = fct_reorder(.f = gene, .x = mean_estimate, .fun = function(x) x)) %>%
  ggplot() +
  geom_point(mapping = aes(x = gene, y = mean_estimate, color = mean_pval, size = n)) +
  geom_linerange(mapping = aes(x = gene, ymin = mean_estimate-(sd_estimate*2), ymax = mean_estimate+(sd_estimate*2), color = mean_pval), size=0.3) +
  scale_color_viridis() +
  scale_size(limits = c(2, 8), breaks = c(4, 6, 8)) +
  coord_flip() +
  facet_wrap(~ set, scales = "free") +
  labs(x = "", y = "Mean effect size (mean Beta)", color = "Mean Adj. P-value", size = "Number of\nassociations")


plots <- plot_grid(plot_kinases, plot_genes)


#+ fig.width=8, fig.height=8
plots

ggsave(filename = "kinase_association_pairs_most_frequent_mut_genes_kinases.png", plot = plots, path = "./output/plots/genetic_associations/", width = 8, height = 8)
unlink("kinase_association_pairs_most_frequent_mut_genes_kinases.png")


#' compare distribution of edge weights between significant and non-significant associations
p_values <- associations %>%
  filter(!is.na(weight)) %>%
  mutate(signf = if_else(p.adjust < 0.05, "Adj.P < 5%", "Adj.P > 5%")) %>%
  select(signf, weight) %>%
  nest() %>%
  mutate(p_value = map_dbl(.x = data, .f = ~ wilcox.test(formula = weight ~ signf, data = .x)$p.value)) %>%
  select(-data) %>%
  mutate(y.position = 380) %>%
  mutate(group1 = 1.2, group2 = 5) %>%
  mutate(p_value = formatC(p_value, digits = 1, format = "e"))

plot <- associations %>%
  filter(!is.na(weight)) %>%
  mutate(signf = if_else(p.adjust < 0.05, "Adj.P < 5%", "Adj.P > 5%")) %>%
  mutate(signf = fct_relevel(signf, "Adj.P > 5%")) %>%
  ggplot(mapping = aes(x = signf, y = weight, fill = signf)) +
  geom_boxplot(alpha = 0.8, notch = F, outlier.size = 0.1, size = 0.2, outlier.shape = NA) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 16, colour = "black"),
    legend.text = element_text(size = 14, colour = "black"),
    legend.title = element_text(size = 16, colour = "black")) +
  stat_pvalue_manual(data = p_values, label = "Wilcoxon, p = {p_value}", xmax = NULL, size = 4, inherit.aes = F) +
  scale_y_continuous(limits = c(NA, 450), name = "String network edge weight") +
  scale_fill_manual(values = c("#bdbdbd", "#67a9cf"), name = "Associations") +
  coord_flip()

#+ fig.width=5, fig.height=1.5
plot

ggsave(filename = "kinase_association_pairs_weights_distribution.png", plot = plot, path = "./output/plots/genetic_associations/", width = 5, height = 1.5)
ggsave(filename = "kinase_association_pairs_weights_distribution.pdf", plot = plot, path = "./output/plots/genetic_associations/", width = 5, height = 1.5)
unlink("kinase_association_pairs_weights_distribution.png")
unlink("kinase_association_pairs_weights_distribution.pdf")


plot <- associations %>%
  mutate(signf = if_else(p.adjust < 0.05, "2", "1")) %>%
  ggplot(mapping = aes(x = signf, y = weight, fill = signf)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, width = 0.1) +
  theme_classic() +
  stat_compare_means(label.y = 600, label.x = 2.4, size = 3) +
  scale_x_discrete(labels = c("P > 5%", "P < 5%"), name = "Adjusted P-value") +
  scale_y_continuous(name = "String network edge weight") +
  scale_fill_manual(values = c("#4292c6", "#9ecae1"), guide = F) +
  coord_flip()

#+ fig.width=4, fig.height=2
plot

ggsave(filename = "kinase_association_pairs_weights_distribution2.png", plot = plot, path = "./output/plots/genetic_associations/", width = 4, height = 2)
unlink("kinase_association_pairs_weights_distribution2.png")


plot <- associations %>%
  mutate(signf = if_else(p.adjust < 0.05, "2", "1")) %>%
  ggplot(mapping = aes(x = weight)) +
  geom_histogram(mapping = aes(y = stat(density), fill = signf), bins = 30, position = "dodge", alpha = 1) +
  geom_density(mapping = aes(color = signf, y = stat(density)), size = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("#4292c6", "#9ecae1"), labels = c("Non-significant", "Significant"), name = "Associations") +
  scale_color_manual(values = c("#4292c6", "#9ecae1"), labels = c("Non-significant", "Significant"), name = "Associations") +
  labs(y = "Density")

#+ fig.width=6, fig.height=3
plot

ggsave(filename = "kinase_association_pairs_weights_distribution3.png", plot = plot, path = "./output/plots/genetic_associations/", width = 6, height = 3)
unlink("kinase_association_pairs_weights_distribution3.png")


p_value <- associations %>%
  filter(!is.na(weight)) %>%
  mutate(signf = if_else(p.adjust < 0.05, "Adj.P < 5%", "Adj.P > 5%")) %>%
  select(signf, weight)
p_value <- formatC(wilcox.test(formula = weight ~ signf, data = p_value)$p.value, digits = 1, format = "e")

means <- associations %>%
  filter(!is.na(weight)) %>%
  mutate(signf = if_else(p.adjust < 0.05, "Adj.P < 5%", "Adj.P > 5%")) %>%
  group_by(signf) %>%
  summarise(weight = mean(weight)) %>%
  ungroup()

plot <- associations %>%
  filter(!is.na(weight)) %>%
  mutate(signf = if_else(p.adjust < 0.05, "Adj.P < 5%", "Adj.P > 5%")) %>%
  mutate(signf = fct_relevel(signf, "Adj.P > 5%")) %>%
  ggplot() +
  #geom_histogram(mapping = aes(x = weight, y = stat(ndensity), fill = signf), bins = 30, position = "dodge", alpha = 0.8) +
  geom_density(mapping = aes(x = weight, color = signf, fill = signf, alpha = signf, y = stat(scaled))) +
  geom_vline(data = means, mapping = aes(xintercept = weight, color = signf), linetype = "dashed", show.legend = F, size = 1) +
  annotate("text", x = 700, y = 0.8, label = str_c("Wilcoxon, p = ", p_value, sep = ""), size = 4) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14, colour = "black"),
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(size = 14, colour = "black")) +
  #scale_fill_manual(values = c("#4292c6", "#9ecae1"), name = "Associations") +
  #scale_color_manual(values = c("#4292c6", "#9ecae1"), name = "Associations") +
  scale_fill_viridis(discrete = T, name = "Associations") +
  scale_color_viridis(discrete = T, name = "Associations") +
  scale_alpha_manual(values = c(0.8,0.5), guide = F) +
  labs(x = "String network edge weight", y = "Density")

#+ fig.width=5, fig.height=1.5
plot

ggsave(filename = "kinase_association_pairs_weights_distribution4.png", plot = plot, path = "./output/plots/genetic_associations/", width = 5, height = 1.5)
ggsave(filename = "kinase_association_pairs_weights_distribution4.pdf", plot = plot, path = "./output/plots/genetic_associations/", width = 5, height = 1.5)
unlink("kinase_association_pairs_weights_distribution4.png")
unlink("kinase_association_pairs_weights_distribution4.pdf")


#' plot examples of associations
examples <- associations %>%
  mutate(pair = pmap_chr(.l = ., .f = ~ if(is.na(..8)){NA}else if(..7 < 0.001 & ..8 >= 800){str_c(..4, ..3, sep = "--")}else{NA})) %>%
  filter(!is.na(pair)) %>%
  mutate(pair = fct_reorder(pair, p.adjust, median)) %>%
  #mutate(pair = fct_recode(pair, `STK11--\nPRKACA` = "STK11--PRKACA", `TP53--\nCDK1` = "TP53--CDK1", `PTEN--\nAKT3` = "PTEN--AKT3", `TP53--\nMAPK13` = "TP53--MAPK13")) %>%
  inner_join(ka, by = "kinase") %>%
  inner_join(mut, by = c("gene", "sample")) %>%
  mutate(mutated = fct_relevel(as.character(mutated), "0", "1"))

helper <- function(vect, lims){
  x <- vect
  x <- x[x >= min(lims) & x <= max(lims)]
  
  return(x)
}

examples2 <- examples %>%
  group_by(pair, mutated) %>%
  summarise(log10P = list(log10P)) %>%
  ungroup() %>%
  #mutate(limits = list(c(-3,3), c(-3,3), c(-3,3), c(-3,3), c(-3,3), c(-3,3), c(-1,1), c(-1,1))) %>%
  mutate(limits = map(.x = log10P, .f = ~ c(quantile(.x, 0.25)-(1.5*IQR(.x)), quantile(.x, 0.75)+(1.5*IQR(.x))))) %>%
  mutate(log10P = map2(.x = log10P, .y = limits, .f = helper)) %>%
  select(-limits) %>%
  unnest()

helper <- function(vect, value){
  x <- vect
  #x[x == min(x)] <- min(x)-value
  x[x == max(x)] <- max(x)+value
  
  return(x)
}

dummy <- examples2 %>%
  group_by(pair, mutated) %>%
  summarise(log10P = list(range(log10P))) %>%
  ungroup() %>%
  mutate(increment = rep(0.1,8)) %>%
  mutate(log10P = map2(.x = log10P, .y = increment, .f = helper)) %>%
  unnest()

p_values <- examples %>%
  select(pair, mutated, log10P) %>%
  group_by(pair) %>%
  nest() %>%
  ungroup() %>%
  mutate(p_value = map_dbl(.x = data, .f = ~ wilcox.test(log10P ~ mutated, data = .x)$p.value)) %>%
  select(-data) %>%
  inner_join(
    dummy %>% select(-increment) %>% group_by(pair) %>% summarise(log10P = max(log10P)) %>% ungroup(),
    by = "pair") %>%
  rename(y.position = log10P) %>%
  mutate(group1 = "0", group2 = "1") %>%
  mutate(p_value = signif(p_value,2))

N <- examples %>%
  #filter(pair %in% c("STK11--PRKACA", "TP53--CDK1")) %>%
  group_by(pair, mutated) %>%
  tally()

plot1 <- examples2 %>%
  #filter(pair %in% c("STK11--PRKACA", "TP53--CDK1")) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = pair, y = log10P, fill = mutated), notch = FALSE, outlier.shape = NA, alpha = 0.8, size = 1, fatten = 1.5) +
  geom_point(mapping = aes(x = pair, y = log10P, fill = mutated), position = position_jitterdodge(jitter.width = 0.2), size = 1.5, show.legend = F, pch = 21, color = "white", alpha = 0.5) +
  #geom_violin(mapping = aes(x = pair, y = log10P, fill = mutated), alpha = 0.8, size = 1) +
  #geom_jitter(mapping = aes(x = pair, y = log10P, alpha = mutated, group = mutated), width = 0.1, size = 0.5, show.legend = F, color = "black") +
  #geom_blank(data = dummy %>% filter(pair %in% c("STK11--PRKACA", "TP53--CDK1")), mapping = aes(x = mutated, y = log10P)) +
  #facet_wrap(~ pair, scales = "fixed", nrow = 1) +
  #stat_compare_means(mapping = aes(x = mutated, y = log10P), method = "t.test", size = 3.5, label = "p.format") +
  #stat_pvalue_manual(data = p_values %>% filter(pair %in% c("STK11--PRKACA", "TP53--CDK1")), label = "p = {p_value}", xmax = NULL, size = 3.5) +
  geom_text(data = p_values, mapping = aes(x = pair, y = Inf, label = paste0("p = ", p_value)), vjust = 1, size = 4) +
  geom_text(data = N, mapping = aes(x = pair, y = -4.2, label = n, group = mutated), size = 4, position = position_dodge(width = 0.8)) +
  #scale_alpha_manual(values = c(0.5, 0.5)) +
  #scale_fill_manual(values = c("#1b9e77", "#d7191c")) +
  #scale_fill_manual(values = c("#bdbdbd", "#67a9cf"), labels = c("0" = "Non-mutated", "1" = "Mutated")) +
  scale_fill_manual(values = c("#d8b365", "#5ab4ac"), labels = c("0" = "Non-mutated", "1" = "Mutated")) +
  scale_x_discrete(labels = c("STK11--PRKACA" = "STK11\nPRKACA", "TP53--CDK1" = "TP53\nCDK1", "PTEN--AKT3" = "PTEN\nAKT3", "TP53--MAPK13" = "TP53\nMAPK13")) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(colour = "black", size = 18),
    axis.text = element_text(colour = "black", size = 16),
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", size = 16),
    strip.text = element_text(colour = "black", size = 12),
    strip.background = element_blank(),
    legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(size=0.5))) +
  labs(x = "Association", y = "Kinase activity score")

#+ fig.width=5, fig.height=5
plot1

ggsave(filename = "kinase_association_pairs_examples1.png", plot = plot1, path = "./output/plots/genetic_associations/", width = 5, height = 5)
ggsave(filename = "kinase_association_pairs_examples1.pdf", plot = plot1, path = "./output/plots/genetic_associations/", width = 5, height = 5)


#plot <- plotKinase(mut, ka, "PRKACA", "STK11")

#+ fig.width=3, fig.height=3
#plot

