library(tidyverse)
library(ggrepel)
library(ggpubr)
library(viridis)
library(RColorBrewer)

source("./src/utils/getSamples.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load samples annotation
samples_annotation <- read_tsv(file = "./output/files/all_samples_annotation.txt")
samples_annotation <- getSamples(samples_annotation, c("phosphorylation"))


sortPrt <- function(x, y){
  sorted <- sort(c(x,y))
  tibble(a=sorted[1],b=sorted[2])
}


# load string network (vertices names sorted for each edge; duplicates removed)
string_net <- read_tsv("./output/files/string_network_sorted_pairs.txt.gz")


# load kinase-TFs associations
# add string network weight
associations <- read_tsv(file = "./output/files/kinaseNotImputed_TF_activity_associations.txt.gz") %>%
  mutate(x = map2(kinase, tf, sortPrt)) %>%
  unnest() %>%
  select(a, b, everything()) %>%
  left_join(string_net, by = c("a", "b"))


# volcano plot of kinase-TFs associations
vplot <- associations %>%
  mutate(l = pmap_chr(.l = ., .f = ~ if(-log10(..7) > 10){str_c(..3, ..4, sep="--")}else{NA})) %>%
  ggplot(mapping = aes(x = kin_beta, y = -log10(padj), color = weight, label = l)) +
  geom_point(size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, color = "grey60") + 
  geom_text_repel(size = 2, segment.size = 0.2, color = "black") +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 12),
    legend.key.size = unit(0.4, "cm")) +
  labs(x = expression(beta ~ "(kinase activity)"), y = "Adjusted P-value (-log10)", color = "Edge\nweight")

# ggsave(filename = "kinase_TFs_associations_volcano_plot.png", plot = vplot, path = "./output/plots/kinase_TF_pathways_associations/", width = 6, height = 4)
# ggsave(filename = "kinase_TFs_associations_volcano_plot.pdf", plot = vplot, path = "./output/plots/kinase_TF_pathways_associations/", width = 6, height = 4)
# unlink("kinase_TFs_associations_volcano_plot.png")
# unlink("kinase_TFs_associations_volcano_plot.pdf")

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

vplot <- associations %>%
  mutate(pair = pmap_chr(.l = ., .f = ~ if(is.na(..8)){NA}else if(..7 < 1e-10 & ..8 >= 800){str_c(..3, ..4, sep = "--")}else{NA})) %>%
  #mutate(pair = map_chr(.x = pair, .f = ~ if(is.na(.x)){NA}else if(.x %in% c("CDK1--E2F4", "CDK1--TFDP1", "CDK2--MYC", "CDK2--FOXM1")){.x}else{NA})) %>%
  mutate(weight = weight/1000) %>%
  mutate(string = map_chr(.x = weight, .f = helper)) %>%
  mutate(string = fct_relevel(string, "Unknown", "(0,0.4)", "[0.4,0.8)", "[0.8,1)")) %>%
  mutate(padj = -log10(padj)) %>%
  ggplot(mapping = aes(x = kin_beta, y = padj, color = string, label = pair)) +
  geom_point(mapping = aes(alpha = string), size = 1) +
  geom_text_repel(size = 3.5, segment.size = 0.1, color = "black", direction = "y", hjust = -0.3, seed = 123) +
  geom_vline(xintercept = 0, color = "black", linetype = "longdash", size = 0.1) +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "longdash", size = 0.1) +
  theme_classic() +
  theme(
    plot.margin = unit(c(0.1,0.1,1.7,0.1), "cm"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    axis.title = element_text(colour = "black", size = 16),
    axis.text = element_text(colour = "black", size = 14),
    legend.text = element_text(colour = "black", size = 14),
    legend.title = element_text(colour = "black", size = 16),
    legend.position = c(0,0),
    legend.justification = c(0.14,1.8),
    legend.direction = "horizontal") +
  scale_alpha_manual(values = c(0.5, 0.8, 0.9, 1), guide = F) +
  scale_color_manual(values = c("#636363", "#6baed6", "#2171b5", "#08306b")) +
  #scale_x_continuous(limits = c(-2.6,2.6), breaks = seq(-2.5,2.5,1)) +
  labs(x = expression("Kinase activity" ~ beta ~ "(effect size)"), y = "Adjusted P-value (-log10)", color = "String network\n edge weight") +
  guides(color = guide_legend(override.aes = list(alpha=1, size = 2), nrow = 2))

# ggsave(filename = "kinase_TFs_associations_volcano_plot2.png", plot = vplot, path = "./output/plots/kinase_TF_pathways_associations/", width = 4, height = 4)
# ggsave(filename = "kinase_TFs_associations_volcano_plot2.pdf", plot = vplot, path = "./output/plots/kinase_TF_pathways_associations/", width = 4, height = 4)
# unlink("kinase_TFs_associations_volcano_plot2.png")
# unlink("kinase_TFs_associations_volcano_plot2.pdf")


# compare distribution of edge weights between significant and non-significant associations
plot <- associations %>%
  mutate(signf = if_else(padj < 0.05 & abs(kin_beta) > 0.5, "2", "1")) %>%
  ggplot(mapping = aes(x = signf, y = weight, fill = signf)) +
  geom_boxplot(notch = T, outlier.size = 0.1, size = 0.2, outlier.shape = NA) +
  #geom_jitter(width = 0.1, alpha = 0.1, color = "grey") +
  theme_minimal() +
  stat_compare_means(label.y = 800, label.x = 2.5, size = 3) +
  scale_x_discrete(labels = c("Non-significant", "Significant"), name = "Associations") +
  scale_y_continuous(name = "String network edge weight") +
  scale_fill_manual(values = c("#4292c6", "#9ecae1"), guide = F) +
  coord_flip()

# ggsave(filename = "kinase_TFs_associations_weights_distribution1.png", plot = plot, path = "./output/plots/kinase_TF_pathways_associations/", width = 6, height = 2)
# unlink("kinase_TFs_associations_weights_distribution1.png")

plot <- associations %>%
  mutate(signf = if_else(padj < 0.05 & abs(kin_beta) > 0.5, "2", "1")) %>%
  ggplot(mapping = aes(x = weight)) +
  geom_histogram(mapping = aes(y = ..density.., fill = signf), bins = 30, position = "dodge", alpha = 1) +
  geom_density(mapping = aes(color = signf), size = 0.5) +
  theme_minimal() +
  scale_fill_manual(values = c("#4292c6", "#9ecae1"), labels = c("Non-significant", "Significant"), name = "Associations") +
  scale_color_manual(values = c("#4292c6", "#9ecae1"), labels = c("Non-significant", "Significant"), name = "Associations")

# ggsave(filename = "kinase_TFs_associations_weights_distribution2.png", plot = plot, path = "./output/plots/kinase_TF_pathways_associations/", width = 8, height = 4)
# unlink("kinase_TFs_associations_weights_distribution2.png")

plot <- associations %>%
  mutate(signf = if_else(padj < 0.05 & abs(kin_beta) > 0.5, "2", "1")) %>%
  ggplot(mapping = aes(x = signf, y = weight, fill = signf)) +
  geom_boxplot(notch = F, outlier.shape = NA, alpha = 1) +
  #geom_jitter(width = 0.1, alpha = 0.1, color = "grey") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", size = 10),
    #legend.spacing.x = unit(0.05, 'cm'),
    legend.margin = margin(c(0,0,-0.2,0), unit="cm")) +
  stat_compare_means(label.y = 730, label.x = 1.4, size = 3.5, aes(label = paste0("p = ", ..p.format..))) +
  scale_y_continuous(name = "String weight") +
  #scale_fill_manual(values = c("#bdbdbd", "#67a9cf"), labels = c("Others", "Adj.P < 5%, |" ~ beta ~ "| > 0.5")) +
  scale_fill_viridis(discrete = T, labels = c("Others", "Adj.P < 5%, |" ~ beta ~ "| > 0.5")) +
  coord_flip() +
  guides(fill = guide_legend(override.aes = list(alpha=1), ncol = 1, reverse=T))

# ggsave(filename = "kinase_TFs_associations_weights_distribution3.png", plot = plot, path = "./output/plots/kinase_TF_pathways_associations/", width = 2, height = 1.8)
# ggsave(filename = "kinase_TFs_associations_weights_distribution3.pdf", plot = plot, path = "./output/plots/kinase_TF_pathways_associations/", width = 2, height = 1.8)
# unlink("kinase_TFs_associations_weights_distribution3.png")
# unlink("kinase_TFs_associations_weights_distribution3.pdf")

plot <- associations %>%
  mutate(signf = if_else(padj < 0.05 & abs(kin_beta) > 0.5, "2", "1")) %>%
  ggplot(mapping = aes(x = signf, y = weight, fill = signf)) +
  geom_boxplot(notch = F, outlier.shape = NA, alpha = 1) +
  #geom_jitter(width = 0.1, alpha = 0.1, color = "grey") +
  theme_classic() +
  theme(
    axis.title.y = element_text(size = 18, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black"),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(size = 18, colour = "black"),
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 16),
    #legend.spacing.x = unit(0.05, 'cm'),
    legend.margin = margin(c(0,0,-0.2,0), unit="cm")) +
  stat_compare_means(label.y = 950, label.x = 0.8, size = 3.5, aes(label = paste0("p = ", ..p.format..))) +
  scale_y_continuous(name = "String weight") +
  scale_fill_manual(labels = c("1" = "Background","2"="Validation"), values = c("1" = "#fdc086", "2" = "#7fc97f")) +
  guides(fill = guide_legend(override.aes = list(alpha=1), ncol = 1, title.position = "top"))

# ggsave(filename = "kinase_TFs_associations_weights_distribution4.png", plot = plot, path = "./output/plots/kinase_TF_pathways_associations/", width = 2, height = 4, units = "in")
# ggsave(filename = "kinase_TFs_associations_weights_distribution4.pdf", plot = plot, path = "./output/plots/kinase_TF_pathways_associations/", width = 2, height = 4, units = "in")
# unlink("kinase_TFs_associations_weights_distribution4.png")
# unlink("kinase_TFs_associations_weights_distribution4.pdf")


# plot some examples

# load transcription factor activities
tfs <- data.table::fread("./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  as_tibble() %>%
  rename(tf=V1) %>%
  pivot_longer(-tf, names_to = "sample", values_to = "tf_activity")

# load kinase activities
kin_activities <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= 3) %>%
  select(-source_type, -n) %>%
  select(sample, kinase, kin_activity=log10P)

# prepare batch covariate to be regressed-out
batch <- samples_annotation %>%
  mutate(batch = map2_chr(batch, tissue, ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(sample, batch)

plot <- associations %>%
  filter(-log10(padj) > 10) %>%
  inner_join(kin_activities, by = "kinase") %>%
  inner_join(tfs, by = c("tf", "sample")) %>%
  inner_join(batch, by = "sample") %>%
  group_by(kinase, tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(tf_activity_resd = map(.x = data, .f = ~ residuals(lm(tf_activity ~ batch, data = .x)))) %>%
  unnest() %>%
  mutate(pair = str_c(kinase, tf, sep = "--")) %>%
  mutate(pair = fct_reorder(pair, padj, .fun = function(x){unique(x)})) %>%
  ggplot(mapping = aes(x = kin_activity, y = tf_activity_resd)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", lwd = 0.5) +
  stat_cor(size = 2) +
  facet_wrap(~ pair, scales = "free")

# ggsave(filename = "kinase_TFs_associations_examples.png", plot = plot, path = "./output/plots/kinase_TF_pathways_associations/", width = 12, height = 8)
# ggsave(filename = "kinase_TFs_associations_examples.pdf", plot = plot, path = "./output/plots/kinase_TF_pathways_associations/", width = 12, height = 8)
# unlink("kinase_TFs_associations_examples.png")
# unlink("kinase_TFs_associations_examples.pdf")


examples <- associations %>%
  mutate(pair = pmap_chr(.l = ., .f = ~ if(is.na(..8)){NA}else if(..7 < 1e-10 & ..8 >= 800){str_c(..3, ..4, sep = "-")}else{NA})) %>%
  #mutate(pair = map_chr(.x = pair, .f = ~ if(is.na(.x)){NA}else if(.x %in% c("CDK1-E2F4", "CDK1-TFDP1", "CDK2-MYC", "CDK2-FOXM1")){.x}else{NA})) %>%
  filter(!is.na(pair)) %>%
  inner_join(kin_activities, by = "kinase") %>%
  inner_join(tfs, by = c("tf", "sample")) %>%
  inner_join(batch, by = "sample") %>%
  group_by(kinase, tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(tf_activity_resd = map(.x = data, .f = ~ residuals(lm(tf_activity ~ batch, data = .x)))) %>%
  mutate(kin_activity_resd = map(.x = data, .f = ~ residuals(lm(kin_activity ~ batch, data = .x)))) %>%
  unnest() %>%
  mutate(pair = fct_reorder(pair, padj, .fun = function(x){unique(x)})) %>%
  select(a,b,kinase,tf,pair,kin_beta,kin_pval,padj,weight,sample,batch,everything())

helper <- function(vect, value){
  x <- vect
  #x[x == min(x)] <- min(x)-value
  x[x == max(x)] <- max(x)+value
  
  return(x)
}

dummy <- examples %>%
  select(pair, tf_activity_resd, kin_activity_resd, kin_activity, tf_activity) %>%
  group_by(pair) %>%
  summarise_if(is.numeric, .funs = ~ list(range(.x))) %>%
  ungroup() %>%
  mutate(increment = rep(0,10)) %>%
  mutate(tf_activity_resd = map2(.x = tf_activity_resd, .y = increment, .f = helper)) %>%
  mutate(tf_activity = map2(.x = tf_activity, .y = increment, .f = helper)) %>%
  #mutate_if(is.list, .funs = function(x, y){map2(.x = x, .y = y, .f = helper)}, y = rep(1,4)) %>%
  #mutate_if(is.list, .funs = ~ map(.x = ..1, .f = helper, value = 1)) %>%
  unnest()

plot <- examples %>%
  ggplot(mapping = aes(x = kin_activity_resd, y = tf_activity_resd)) +
  #ggplot(mapping = aes(x = kin_activity, y = tf_activity_resd)) +
  #ggplot(mapping = aes(x = kin_activity, y = tf_activity)) +
  geom_point(size = 0.1, alpha = 0.4) +
  geom_smooth(method = "lm", lwd = 0.5, color = "red") +
  geom_blank(data = dummy, mapping = aes(x = kin_activity_resd, y = tf_activity_resd)) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 14),
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 16)) +
  stat_cor(size = 3, label.sep = ", ", label.y.npc = 1, label.x.npc = 0.01) +
  facet_wrap(~ pair, scales = "free") +
  labs(x = "Kinase activity", y = "Transcription factor activity")

# ggsave(filename = "kinase_TFs_associations_examples2.png", plot = plot, path = "./output/plots/kinase_TF_pathways_associations/", width = 4, height = 4)
# ggsave(filename = "kinase_TFs_associations_examples2.pdf", plot = plot, path = "./output/plots/kinase_TF_pathways_associations/", width = 4, height = 4)
# unlink("kinase_TFs_associations_examples2.png")
# unlink("kinase_TFs_associations_examples2.pdf")

