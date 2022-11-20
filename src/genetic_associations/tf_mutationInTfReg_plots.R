#' ---
#' title: "Associations between TF activities and mutations in the same TF"
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


#' load TF activity data
tf <- read_csv(file = "./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  rename(tf = X1) %>%
  pivot_longer(-tf, names_to = "sample", values_to = "activity")


#' load mutation-TF activity associations
associations <- read_tsv(file = "./output/files/tf_mutStatus_inTF.txt.gz") %>%
  mutate(signf = map_dbl(.x = p.adjust, .f = ~ if(.x < 0.05){2}else{1}))


#' volcano plot
volcano <- associations %>%
  mutate(l = pmap_chr(.l = ., .f = ~ if(..4 < 0.05){..1}else{NA})) %>%
  mutate(p.adjust = -log10(p.adjust)) %>%
  ggplot(mapping = aes(x = estimate, y = p.adjust, color = as.character(signf), label = l)) +
  geom_point(size = 1) +
  geom_text_repel(size = 3, segment.size = 0.3, color = "black", point.padding = 0.1, box.padding = 1, seed = 1000, direction = "y") +
  geom_vline(xintercept = 0, color = "black", linetype = "longdash", size = 0.2) +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "longdash", size = 0.2) +
  scale_colour_manual(values = c("#662506", "#006d2c"), guide = F) +
  scale_x_continuous(limits = c(-2,2), breaks = seq(-2,2,1)) +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1.5),
    axis.text = element_text(size = 8, colour = "black"),
    axis.title = element_blank())

#+ fig.width=1.5, fig.height=1.5
volcano

ggsave(filename = "tf_association_mutation_in_tf_volcano_plot.png", plot = volcano, path = "./output/plots/genetic_associations/", width = 1.5, height = 1.5, dpi = 900)
ggsave(filename = "tf_association_mutation_in_tf_volcano_plot.pdf", plot = volcano, path = "./output/plots/genetic_associations/", width = 1.5, height = 1.5)
unlink("tf_association_mutation_in_tf_volcano_plot.png")
unlink("tf_association_mutation_in_tf_volcano_plot.pdf")



#' plot examples of associations
examples <- associations %>%
  mutate(pair = pmap_chr(.l = ., .f = ~ if(..4 < 0.05){..1}else{NA})) %>%
  filter(!is.na(pair)) %>%
  mutate(pair = fct_reorder(pair, p.adjust, median)) %>%
  inner_join(tf, by = "tf") %>%
  inner_join(mut, by = c("tf" = "gene", "sample")) %>%
  mutate(mutated = fct_relevel(as.character(mutated), "0", "1"))

helper <- function(vect, value){
  x <- vect
  #x[x == min(x)] <- min(x)-value
  x[x == max(x)] <- max(x)+value
  
  return(x)
}

dummy <- examples %>%
  group_by(pair, mutated) %>%
  summarise(activity = list(range(activity))) %>%
  ungroup() %>%
  mutate(activity = map(.x = activity, .f = helper, value = 1)) %>%
  unnest()

N <- examples %>%
  group_by(pair, mutated) %>%
  tally() %>%
  ungroup()

plot <- examples %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = pair, y = activity, fill = mutated), notch = FALSE, outlier.shape = NA, alpha = 0.8, size = 1, fatten = 1.5) +
  geom_point(mapping = aes(x = pair, y = activity, fill = mutated), position = position_jitterdodge(jitter.width = 0.2), size = 1.5, show.legend = F, pch = 21, color = "white", alpha = 0.5) +
  geom_text(data = N, mapping = aes(x = pair, y = -9, label = n, group = mutated), size = 4, position = position_dodge(width = 0.8)) +
  #geom_boxplot(mapping = aes(x = mutated, y = activity, fill = mutated), notch = FALSE, outlier.shape = NA, alpha = 0.8, show.legend = F) +
  #geom_jitter(mapping = aes(x = mutated, y = activity, alpha = mutated), width = 0.1, size = 0.5, show.legend = F, color = "black") +
  #geom_blank(data = dummy, mapping = aes(x = mutated, y = activity)) +
  #facet_wrap(~ pair, scales = "free", nrow = 1) +
  stat_compare_means(mapping = aes(x = pair, y = activity, group = mutated), method = "wilcox.test", size = 4, label = "p.format", label.y.npc = 1, vjust = -0.5) +
  #scale_alpha_manual(values = c(0.05, 0.1)) +
  #scale_fill_manual(values = c("#1b9e77", "#d7191c")) +
  scale_fill_manual(values = c("#d8b365", "#5ab4ac"), labels = c("0" = "Non-mutated", "1" = "Mutated")) +
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
  labs(x = "Association", y = "TF activity score")

#+ fig.width=5, fig.height=2.5
plot

ggsave(filename = "tf_association_mutation_in_tf_pairs_examples.png", plot = plot, path = "./output/plots/genetic_associations/", width = 5, height = 5)
ggsave(filename = "tf_association_mutation_in_tf_pairs_examples.pdf", plot = plot, path = "./output/plots/genetic_associations/", width = 5, height = 5)
