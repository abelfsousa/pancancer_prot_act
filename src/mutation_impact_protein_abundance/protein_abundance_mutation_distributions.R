#' ---
#' title: "Protein abundance distribution by mutation type"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=6, fig.height=6)


#' load R packages
library(rstatix)
library(ggpubr)
library(tidyverse)
library(viridis)

set.seed(123)

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load mutation data
cptac_mutations <- read_tsv(file = "./output/files/mutations_protpos.txt.gz")


#' select only the gene, sample and type of mutation\
#' remove duplicates (same gene in the same sample can have multiple mutations of same type)\
mut <- cptac_mutations %>%
  select(sample, gene = gene_symbol, mutation_type = variant_class) %>%
  distinct()


#' load protein abundance
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "log2fc") %>%
  filter(!is.na(log2fc))


#' # -- assess the protein abundance distributions between different mutation types

#' filter the mutations that occurred in proteins\
#' select only the pairs sample-protein that are unique to the same type of mutation\
#' this was done to prevent the assignment of the same protein abundance to different types of mutations
mutations <- mut %>%
  semi_join(protein, by = c("gene", "sample")) %>%
  group_by(mutation_type) %>%
  nest() %>%
  ungroup() %>%
  mutate(data2 = (function(x) map2(.x=1:length(x),.y=x,.f=function(a,b)reduce(.x=c(list(b),x[-a]),.f=setdiff))) (data)) %>%
  select(-data) %>%
  unnest()

#' define a background to compare with (no known mutation for the protein in the sample)
background <- protein %>%
  anti_join(mut, by = c("gene", "sample")) %>%
  mutate(mutation_type = "Background (no mutation)")
  #sample_n(nrow(mutations))

distribution <- protein %>%
  inner_join(mutations, by = c("sample", "gene")) %>%
  bind_rows(background) %>%
  mutate(mutation_type = fct_rev(fct_reorder(mutation_type, log2fc, .fun = median))) %>%
  mutate(mutation_type = fct_rev(fct_relevel(mutation_type, "Background (no mutation)"))) %>%
  ggplot(mapping = aes(x=mutation_type, y = log2fc, fill = mutation_type)) +
  geom_boxplot(notch = F, alpha = 0.8, show.legend = F, color = "black", outlier.size = 1) +
  scale_fill_viridis(discrete = T) +
  coord_flip() +
  theme_classic() +
  labs(x = "", y = "Protein abundance (log2fc)")

#+ fig.width=12, fig.height=4
distribution

ggsave(filename = "prot_abundance_mutation_distribution_box1.png", plot = distribution, path = "./output/plots/mutation_impact_protein_abundance/", height = 4, width = 12)


distribution <- protein %>%
  inner_join(mutations, by = c("sample", "gene")) %>%
  bind_rows(background) %>%
  mutate(mutation_type = fct_rev(fct_reorder(mutation_type, log2fc, .fun = median))) %>%
  mutate(mutation_type = fct_rev(fct_relevel(mutation_type, "Background (no mutation)"))) %>%
  ggplot(mapping = aes(x = log2fc, color = mutation_type)) +
  geom_density(alpha = 0.5) +
  scale_color_viridis(discrete = T) +
  theme_classic() +
  labs(x = "Protein abundance (log2fc)", y = "Density")

#+ fig.width=12, fig.height=4
distribution

ggsave(filename = "prot_abundance_mutation_distribution_density1.png", plot = distribution, path = "./output/plots/mutation_impact_protein_abundance/", height = 4, width = 12)


# plot for paper
protein_mutation <- protein %>%
  inner_join(mutations, by = c("sample", "gene")) %>%
  bind_rows(background)

filter_abundances <- function(vect, lims){
  x <- vect
  x <- x[x >= min(lims) & x <= max(lims)]
  
  return(x)
}

protein_mutation2 <- protein_mutation %>%
  group_by(mutation_type) %>%
  summarise(log2fc = list(log2fc)) %>%
  ungroup() %>%
  mutate(limits = map(.x = log2fc, .f = ~ c(quantile(.x, 0.25)-(1.5*IQR(.x)), quantile(.x, 0.75)+(1.5*IQR(.x))))) %>%
  mutate(log2fc = map2(.x = log2fc, .y = limits, .f = filter_abundances)) %>%
  select(-limits) %>%
  unnest()

pvalues_ <- protein_mutation2 %>%
  t_test(log2fc ~ mutation_type, ref.group = "Background (no mutation)", p.adjust.method = "none") %>%
  add_xy_position(x = "mutation_type")

pvalues <- protein_mutation %>%
  t_test(log2fc ~ mutation_type, ref.group = "Background (no mutation)", p.adjust.method = "none") %>%
  add_xy_position(x = "mutation_type") %>%
  inner_join(pvalues_[, c("group2", "y.position")], by = c("group2")) %>%
  rename(y.position = y.position.y)

N <- protein_mutation %>%
  group_by(mutation_type) %>%
  summarise(n = n()) %>%
  ungroup()

N_ <- protein_mutation2 %>%
  group_by(mutation_type) %>%
  summarise(min = min(log2fc)) %>%
  mutate(min = min(min)+0.1) %>%
  ungroup()

N <- N %>% inner_join(N_)


distribution <- protein_mutation2 %>%
  #mutate(mutation_type = fct_rev(fct_reorder(mutation_type, log2fc, .fun = median))) %>%
  #mutate(mutation_type = fct_rev(fct_relevel(mutation_type, "Background (no mutation)"))) %>%
  ggplot(mapping = aes(x=mutation_type, y = log2fc)) +
  geom_boxplot(mapping = aes(fill = mutation_type), alpha = 0.8, show.legend = F, outlier.shape = NA) +
  stat_pvalue_manual(pvalues, label = "p", tip.length = 0.01, coord.flip = TRUE, size = 5, step.increase = 0.03) +
  geom_text(data = N, mapping = aes(x = seq(1.3, 9.3, by = 1), y = min, label = n), size = 7) +
  coord_flip() +
  scale_fill_viridis(discrete = T) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 24, colour = "black", hjust = 0.5),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 24, colour = "black"),
    axis.text.x = element_text(size = 22, colour = "black"),
    axis.text.y = element_text(size = 18, colour = "black")) +
  labs(y = "Protein abundance (log2fc)", title = "CPTAC tumours")

#+ fig.width=12, fig.height=6
distribution

ggsave(filename = "prot_abundance_mutation_distribution_box2.png", plot = distribution, path = "./output/plots/mutation_impact_protein_abundance/", height = 6, width = 12)
ggsave(filename = "prot_abundance_mutation_distribution_box2.pdf", plot = distribution, path = "./output/plots/mutation_impact_protein_abundance/", height = 6, width = 12)


distribution <- protein %>%
  inner_join(mutations, by = c("sample", "gene")) %>%
  bind_rows(background) %>%
  mutate(mutation_type = fct_rev(fct_reorder(mutation_type, log2fc, .fun = median))) %>%
  mutate(mutation_type = fct_rev(fct_relevel(mutation_type, "Background (no mutation)"))) %>%
  ggplot(mapping = aes(x = log2fc, color = mutation_type)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(-4,4)) +
  scale_color_viridis(discrete = T) +
  theme_classic() +
  labs(x = "Protein abundance (log2fc)", y = "Density")

#+ fig.width=12, fig.height=4
distribution

ggsave(filename = "prot_abundance_mutation_distribution_density2.png", plot = distribution, path = "./output/plots/mutation_impact_protein_abundance/", height = 4, width = 12)


distribution <- protein %>%
  inner_join(mutations, by = c("sample", "gene")) %>%
  bind_rows(background) %>%
  mutate(mutation_type = fct_rev(fct_reorder(mutation_type, log2fc, .fun = median))) %>%
  mutate(mutation_type = fct_rev(fct_relevel(mutation_type, "Background (no mutation)"))) %>%
  ggplot(mapping = aes(x = log2fc, color = mutation_type)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(-4,4)) +
  scale_color_viridis(discrete = T) +
  theme_classic() +
  labs(x = "Protein abundance (log2fc)", y = "Density")

#+ fig.width=12, fig.height=4
distribution

ggsave(filename = "prot_abundance_mutation_distribution_density3.pdf", plot = distribution, path = "./output/plots/mutation_impact_protein_abundance/", height = 4, width = 12)
ggsave(filename = "prot_abundance_mutation_distribution_density3.png", plot = distribution, path = "./output/plots/mutation_impact_protein_abundance/", height = 4, width = 12)


#' # -- nonsense vs missense mutations
nonsense <- mut %>%
  filter(mutation_type == "Nonsense_Mutation") %>%
  select(-mutation_type)

missense <- mut %>%
  filter(mutation_type == "Missense_Mutation") %>%
  select(-mutation_type)

exclude <- intersect(missense, nonsense)

nonsense <- nonsense %>%
  anti_join(exclude, by = c("sample", "gene")) %>%
  inner_join(protein, by = c("sample", "gene")) %>%
  mutate(mutation_type = "Nonsense_Mutation")

missense <- missense %>%
  anti_join(exclude, by = c("sample", "gene")) %>%
  inner_join(protein, by = c("sample", "gene")) %>%
  mutate(mutation_type = "Missense_Mutation")

mutations <- bind_rows(nonsense, missense)

mutations <- mut %>%
  filter(mutation_type %in% c("Nonsense_Mutation", "Missense_Mutation")) %>%
  semi_join(protein, by = c("gene", "sample")) %>%
  group_by(mutation_type) %>%
  nest() %>%
  ungroup() %>%
  mutate(data2 = (function(x) map2(.x=1:length(x),.y=x,.f=function(a,b)reduce(.x=c(list(b),x[-a]),.f=setdiff))) (data)) %>%
  select(-data) %>%
  unnest() %>%
  inner_join(protein, by = c("sample", "gene"))

nonsense_missense <- mutations %>%
  mutate(mutation_type = fct_reorder(mutation_type, log2fc, .fun = median)) %>%
  ggplot(mapping = aes(x = log2fc, fill = mutation_type, color = mutation_type)) +
  geom_density(alpha = 0.5) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "Protein abundance (log2fc)", y = "Density")

#+ fig.width=8, fig.height=2
nonsense_missense

ggsave(filename = "prot_abundance_nonsense_missense_density1.png", plot = nonsense_missense, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)


nonsense_missense <- mutations %>%
  mutate(mutation_type = fct_reorder(mutation_type, log2fc, .fun = median)) %>%
  ggplot(mapping = aes(x = mutation_type, y = log2fc, fill = mutation_type)) +
  geom_boxplot(notch = T, alpha = 0.8, color = "grey60") +
  scale_fill_viridis(discrete = T, guide=F) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()) +
  coord_flip() +
  stat_compare_means(method = "t.test", label.x = 2.4, label.y = 10, size = 3) +
  labs(x = "", y = "Protein abundance (log2fc)")

#+ fig.width=8, fig.height=2
nonsense_missense

ggsave(filename = "prot_abundance_nonsense_missense_box1.png", plot = nonsense_missense, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)


nonsense_missense <- mutations %>%
  mutate(mutation_type = fct_reorder(mutation_type, log2fc, .fun = median)) %>%
  ggplot(mapping = aes(x = log2fc, fill = mutation_type, color = mutation_type)) +
  geom_density(alpha = 0.5) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "Protein abundance (log2fc)", y = "Density") +
  scale_x_continuous(limits = c(-2,2))

#+ fig.width=8, fig.height=2
nonsense_missense

ggsave(filename = "prot_abundance_nonsense_missense_density2.png", plot = nonsense_missense, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)


nonsense_missense <- mutations %>%
  mutate(mutation_type = fct_reorder(mutation_type, log2fc, .fun = median)) %>%
  ggplot(mapping = aes(x = mutation_type, y = log2fc, fill = mutation_type)) +
  geom_boxplot(notch = T, alpha = 0.8, color = "grey60") +
  scale_fill_viridis(discrete = T, guide = F) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()) +
  coord_flip() +
  scale_y_continuous(limits = c(-2,2)) +
  stat_compare_means(method = "t.test", label.x = 2.4, label.y = 1, size = 3) +
  labs(x = "", y = "Protein abundance (log2fc)")

#+ fig.width=8, fig.height=2
nonsense_missense

ggsave(filename = "prot_abundance_nonsense_missense_box2.png", plot = nonsense_missense, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)



#' # -- nonsense mutations vs background
nonsense <- mut %>%
  filter(mutation_type == "Nonsense_Mutation") %>%
  inner_join(protein, by = c("sample", "gene"))

background <- protein %>%
  anti_join(mut, by = c("gene", "sample")) %>%
  mutate(mutation_type = "Background (no mutation)")
  #sample_n(nrow(nonsense))

nonsense_backg <- nonsense %>%
  bind_rows(background) %>%
  ggplot(mapping = aes(x = log2fc, fill = mutation_type, color = mutation_type)) +
  geom_density(alpha = 0.5) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "Protein abundance (log2fc)", y = "Density")

#+ fig.width=8, fig.height=2
nonsense_backg

ggsave(filename = "prot_abundance_nonsense_backg_density1.png", plot = nonsense_backg, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)


nonsense_backg <- nonsense %>%
  bind_rows(background) %>%
  ggplot(mapping = aes(x = log2fc, fill = mutation_type, color = mutation_type)) +
  geom_density(alpha = 0.5) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  theme_classic() +
  scale_x_continuous(limits = c(-4,4)) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "Protein abundance (log2fc)", y = "Density")

#+ fig.width=8, fig.height=2
nonsense_backg

ggsave(filename = "prot_abundance_nonsense_backg_density2.png", plot = nonsense_backg, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)


nonsense_backg <- nonsense %>%
  bind_rows(background) %>%
  ggplot(mapping = aes(x = mutation_type, y = log2fc, fill = mutation_type)) +
  geom_boxplot(notch = T, color = "grey60", alpha = 0.8) +
  scale_fill_viridis(discrete = T) +
  stat_compare_means(method = "t.test", label.x = 2.4, label.y = 5, size = 3) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()) +
  labs(x = "", y = "Protein abundance (log2fc)")

#+ fig.width=8, fig.height=2
nonsense_backg

ggsave(filename = "prot_abundance_nonsense_backg_box1.png", plot = nonsense_backg, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)


nonsense_backg <- nonsense %>%
  bind_rows(background) %>%
  ggplot(mapping = aes(x = mutation_type, y = log2fc, fill = mutation_type)) +
  geom_boxplot(notch = T, color = "grey60", alpha = 0.8) +
  scale_fill_viridis(discrete = T) +
  stat_compare_means(method = "t.test", label.x = 2.4, label.y = 2, size = 3) +
  coord_flip() +
  theme_classic() +
  scale_y_continuous(limits = c(-4,4)) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()) +
  labs(x = "", y = "Protein abundance (log2fc)")

#+ fig.width=8, fig.height=2
nonsense_backg

ggsave(filename = "prot_abundance_nonsense_backg_box2.png", plot = nonsense_backg, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)



#' # -- frameshift mutations vs background
frameshift <- mut %>%
  filter(str_detect(mutation_type, "Frame_Shift")) %>%
  mutate(mutation_type = "Frame_Shift") %>%
  inner_join(protein, by = c("sample", "gene"))

background <- protein %>%
  anti_join(mut, by = c("gene", "sample")) %>%
  mutate(mutation_type = "Background (no mutation)")
  #sample_n(nrow(frameshift))

frameshfit_backg <- frameshift %>%
  bind_rows(background) %>%
  ggplot(mapping = aes(x = log2fc, fill = mutation_type, color = mutation_type)) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "Protein abundance (log2fc)", y = "Density")

#+ fig.width=8, fig.height=2
frameshfit_backg

ggsave(filename = "prot_abundance_frameshift_backg_density1.png", plot = frameshfit_backg, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)


frameshfit_backg <- frameshift %>%
  bind_rows(background) %>%
  ggplot(mapping = aes(x = log2fc, fill = mutation_type, color = mutation_type)) +
  geom_density(alpha = 0.5) +
  scale_color_viridis(discrete = T) +
  scale_fill_viridis(discrete = T) +
  theme_classic() +
  scale_x_continuous(limits = c(-4,4)) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "Protein abundance (log2fc)", y = "Density")

#+ fig.width=8, fig.height=2
frameshfit_backg

ggsave(filename = "prot_abundance_frameshift_backg_density2.png", plot = frameshfit_backg, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)


frameshfit_backg <- frameshift %>%
  bind_rows(background) %>%
  ggplot(mapping = aes(x = mutation_type, y = log2fc, fill = mutation_type)) +
  geom_boxplot(notch = T, color = "grey60", alpha = 0.8) +
  scale_fill_viridis(discrete = T) +
  coord_flip() +
  theme_classic() +
  stat_compare_means(method = "t.test", label.x = 2.4, label.y = 5, size = 3) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()) +
  labs(x = "", y = "Protein abundance (log2fc)")

#+ fig.width=8, fig.height=2
frameshfit_backg

ggsave(filename = "prot_abundance_frameshift_backg_box1.png", plot = frameshfit_backg, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)


frameshfit_backg <- frameshift %>%
  bind_rows(background) %>%
  ggplot(mapping = aes(x = mutation_type, y = log2fc, fill = mutation_type)) +
  geom_boxplot(notch = T, color = "grey60", alpha = 0.8) +
  scale_fill_viridis(discrete = T) +
  coord_flip() +
  theme_classic() +
  scale_y_continuous(limits = c(-4,4)) +
  stat_compare_means(method = "t.test", label.x = 2.4, label.y = 2, size = 3) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()) +
  labs(x = "", y = "Protein abundance (log2fc)")

#+ fig.width=8, fig.height=2
frameshfit_backg

ggsave(filename = "prot_abundance_frameshift_backg_box2.png", plot = frameshfit_backg, path = "./output/plots/mutation_impact_protein_abundance/", height = 2, width = 8)

