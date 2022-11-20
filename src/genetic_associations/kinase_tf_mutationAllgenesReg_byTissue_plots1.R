#

# load R packages
library(tidyverse)
library(UpSetR)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load mutation to kinase/TF activity associations
kinase_associations <- map_dfr(
  .x = c("5" = 1, "4" = 2, "3" = 3, "10" = 4),
  .f = ~ read_tsv(file = str_c("./output/files/ka_mutStatus_allGenes_byTissue", .x, ".txt.gz"), progress = F, col_types = cols()),
  .id = "mutations") %>%
  filter(p.adjust < 0.05)

tf_associations <- map_dfr(
  .x = c("5" = 1, "4" = 2, "3" = 3, "10" = 4),
  .f = ~ read_tsv(file = str_c("./output/files/tf_mutStatus_allGenes_byTissue", .x, ".txt.gz"), progress = F, col_types = cols()),
  .id = "mutations") %>%
  filter(p.adjust < 0.05)


# plot number of associations by tissues and dataset (batch)
kinase_associations %>%
  mutate(tissue = fct_rev(fct_infreq(tissue)), batch = fct_rev(fct_infreq(batch)), mutations = fct_inseq(mutations)) %>%
  ggplot(mapping = aes(x = tissue, fill = batch)) +
  geom_bar(position = "dodge") +
  facet_wrap(facets = vars(mutations), scales = "fixed") +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE))
ggsave("by_tissue_kinase_mutation_associations_barplot.png", path = "./output/plots/genetic_associations/", width = 10, height = 8)
ggsave("by_tissue_kinase_mutation_associations_barplot.pdf", path = "./output/plots/genetic_associations/", width = 10, height = 8)

tf_associations %>%
  mutate(tissue = fct_rev(fct_infreq(tissue)), batch = fct_rev(fct_infreq(batch)), mutations = fct_inseq(mutations)) %>%
  ggplot(mapping = aes(x = tissue, fill = batch)) +
  geom_bar(position = "dodge") +
  facet_wrap(facets = vars(mutations), scales = "fixed") +
  coord_flip() +
  guides(fill = guide_legend(reverse = TRUE))
ggsave("by_tissue_tf_mutation_associations_barplot.png", path = "./output/plots/genetic_associations/", width = 10, height = 8)
ggsave("by_tissue_tf_mutation_associations_barplot.pdf", path = "./output/plots/genetic_associations/", width = 10, height = 8)


# plot the highest number of associations obtained
kinase_associations %>%
  filter(mutations == 4) %>%
  mutate(tissue = fct_rev(fct_infreq(tissue)), batch = fct_rev(fct_infreq(batch)), mutations = fct_inseq(mutations)) %>%
  ggplot(mapping = aes(x = tissue, fill = batch)) +
  geom_bar(position = "dodge") +
  facet_wrap(facets = vars(mutations), scales = "fixed", labeller = labeller(mutations = c("4" = ">= 4 samples (mutation)"))) +
  coord_flip() +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "right") +
  guides(fill = guide_legend(reverse = TRUE))
ggsave("by_tissue_kinase_mutation_4_associations_barplot.png", path = "./output/plots/genetic_associations/", width = 6, height = 4)
ggsave("by_tissue_kinase_mutation_4_associations_barplot.pdf", path = "./output/plots/genetic_associations/", width = 6, height = 4)

tf_associations %>%
  filter(mutations == 10) %>%
  mutate(tissue = fct_rev(fct_infreq(tissue)), batch = fct_rev(fct_infreq(batch)), mutations = fct_inseq(mutations)) %>%
  ggplot(mapping = aes(x = tissue, fill = batch)) +
  geom_bar(position = "dodge") +
  facet_wrap(facets = vars(mutations), scales = "fixed", labeller = labeller(mutations = c("10" = ">= 10 samples (mutation)"))) +
  coord_flip() +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    legend.position = "right") +
  guides(fill = guide_legend(reverse = TRUE))
ggsave("by_tissue_tf_mutation_10_associations_barplot.png", path = "./output/plots/genetic_associations/", width = 6, height = 4)
ggsave("by_tissue_tf_mutation_10_associations_barplot.pdf", path = "./output/plots/genetic_associations/", width = 6, height = 4)


# plot the intersection of the associations between the different tissues/datasets
kinases <- kinase_associations %>%
  filter(mutations == 4) %>%
  select(tissue, batch, kinase, gene) %>%
  unite(col = "kinase_gene", kinase, gene, sep = "_") %>%
  group_by(tissue, batch) %>%
  summarise(data = list(kinase_gene)) %>%
  ungroup() %>%
  mutate(data = set_names(data, nm = str_c(tissue, batch, sep = "_")))

association_sets <- upset(fromList(kinases$data), order.by = "freq")
ggsave(filename = "by_tissue_kinase_mutation_4_associations_upsetplot.png", plot = print(association_sets), path = "./output/plots/genetic_associations/", height = 4, width = 4)
ggsave(filename = "by_tissue_kinase_mutation_4_associations_upsetplot.pdf", plot = print(association_sets), path = "./output/plots/genetic_associations/", height = 4, width = 4)

kinases <- kinase_associations %>%
  filter(mutations == 10) %>%
  select(tissue, batch, kinase, gene) %>%
  unite(col = "kinase_gene", kinase, gene, sep = "_") %>%
  group_by(tissue, batch) %>%
  summarise(data = list(kinase_gene)) %>%
  ungroup() %>%
  mutate(data = set_names(data, nm = str_c(tissue, batch, sep = "_")))

association_sets <- upset(fromList(kinases$data), order.by = "freq")
ggsave(filename = "by_tissue_kinase_mutation_10_associations_upsetplot.png", plot = print(association_sets), path = "./output/plots/genetic_associations/", height = 4, width = 4)
ggsave(filename = "by_tissue_kinase_mutation_10_associations_upsetplot.pdf", plot = print(association_sets), path = "./output/plots/genetic_associations/", height = 4, width = 4)

tfs <- tf_associations %>%
  filter(mutations == 10) %>%
  select(tissue, batch, tf, gene) %>%
  unite(col = "tf_gene", tf, gene, sep = "_") %>%
  group_by(tissue, batch) %>%
  summarise(data = list(tf_gene)) %>%
  ungroup() %>%
  mutate(data = set_names(data, nm = str_c(tissue, batch, sep = "_")))

association_sets <- upset(fromList(tfs$data), order.by = "freq", nsets = 9)
ggsave(filename = "by_tissue_tf_mutation_10_associations_upsetplot.png", plot = print(association_sets), path = "./output/plots/genetic_associations/", height = 4, width = 4)
ggsave(filename = "by_tissue_tf_mutation_10_associations_upsetplot.pdf", plot = print(association_sets), path = "./output/plots/genetic_associations/", height = 4, width = 4)


# get the kinase-mutation associations that are specific of each tissue
specific_associations <- function(x){
  associations <- x
  res <- map2(.x=seq_along(associations), .y=associations, .f=function(a,b){reduce(.x=c(list(b),associations[-a]),.f=setdiff)})
  return(res)
}

kin_spec_by_tissue1 <- kinase_associations %>%
  filter(mutations == 4) %>%
  select(tissue, batch, kinase, gene) %>%
  group_by(tissue, batch) %>%
  nest() %>%
  #summarise(data = list(tibble(kinase, gene))) %>%
  #summarise(data = (function(x, y){list(tibble(x, y))}) (kinase, gene)) %>%
  #summarise(data = list((function(x, y){tibble(x, y)}) (kinase, gene))) %>%
  ungroup() %>%
  #mutate(data2 = (function(x) map2(.x=seq_along(x),.y=x,.f=function(a,b)reduce(.x=c(list(b),x[-a]),.f=setdiff))) (data))
  mutate(data2 = specific_associations(data)) %>%
  select(-data) %>%
  unnest(data2)

kin_spec_by_tissue2 <- kinase_associations %>%
  filter(mutations == 10) %>%
  select(tissue, batch, kinase, gene) %>%
  group_by(tissue, batch) %>%
  nest() %>%
  #summarise(data = list(tibble(kinase, gene))) %>%
  #summarise(data = (function(x, y){list(tibble(x, y))}) (kinase, gene)) %>%
  #summarise(data = list((function(x, y){tibble(x, y)}) (kinase, gene))) %>%
  ungroup() %>%
  #mutate(data2 = (function(x) map2(.x=seq_along(x),.y=x,.f=function(a,b)reduce(.x=c(list(b),x[-a]),.f=setdiff))) (data))
  mutate(data2 = specific_associations(data)) %>%
  select(-data) %>%
  unnest(data2)

tf_spec_by_tissue <- tf_associations %>%
  filter(mutations == 10) %>%
  select(tissue, batch, tf, gene) %>%
  group_by(tissue, batch) %>%
  nest() %>%
  #summarise(data = list(tibble(tf, gene))) %>%
  #summarise(data = (function(x, y){list(tibble(x, y))}) (tf, gene)) %>%
  #summarise(data = list((function(x, y){tibble(x, y)}) (tf, gene))) %>%
  ungroup() %>%
  #mutate(data2 = (function(x) map2(.x=seq_along(x),.y=x,.f=function(a,b)reduce(.x=c(list(b),x[-a]),.f=setdiff))) (data))
  mutate(data2 = specific_associations(data)) %>%
  select(-data) %>%
  unnest(data2)
