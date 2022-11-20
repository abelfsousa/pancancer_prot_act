# assess the kinase/TF activity distributions between different mutation types

# load R packages
library(ggpubr)
library(rstatix)
library(viridis)
library(tidyverse)

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# -- for CPTAC tumours


# load CPTAC samples
source("./src/utils/getSamples.R")
cptac_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
  getSamples("protein") %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue)


# load mutation data
# select only the gene, sample and type of mutation
# remove duplicates (same gene in the same sample can have multiple mutations of same type)
cptac_mutations <- read_tsv(file = "./output/files/mutations_protpos.txt.gz") %>%
  select(sample, gene = gene_symbol, mutation_type = variant_class) %>%
  distinct()


# load CPTAC TF activities
tf_activity <- read_csv(file = "./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "activity") %>%
  rename(tf=X1) %>%
  select(sample, gene = tf, activity)


# load CPTAC kinase activities (with quantile-normalized protein regressed-out phosphorylation data)
k_subN <- 3
kin_activity <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  select(sample, gene = kinase, activity = log10P)


# filter the mutations that occurred in kinases
# select only the pairs sample-kinase that are unique to the same type of mutation
# this was done to prevent the assignment of the same kinase activity to different types of mutations
kin_mutations <- cptac_mutations %>%
  semi_join(kin_activity, by = c("gene", "sample")) %>%
  group_by(mutation_type) %>%
  nest() %>%
  ungroup() %>%
  mutate(data2 = (function(x) map2(.x=1:length(x),.y=x,.f=function(a,b)reduce(.x=c(list(b),x[-a]),.f=setdiff))) (data)) %>%
  select(-data) %>%
  filter(map_lgl(data2, ~ nrow(.x) > 10 )) %>%
  unnest() %>%
  mutate(protein = "Kinases")

# define a background to compare with (no known mutation for the kinase in the sample)
kin_background <- kin_activity %>%
  anti_join(cptac_mutations, by = c("gene", "sample")) %>%
  mutate(mutation_type = "Background (no mutation)") %>%
  mutate(protein = "Kinases")

kinase_mutation_data <- kin_activity %>%
  inner_join(kin_mutations, by = c("sample", "gene")) %>%
  bind_rows(kin_background)

distribution <- kinase_mutation_data %>%
  mutate(mutation_type = fct_rev(fct_reorder(mutation_type, activity, .fun = median))) %>%
  mutate(mutation_type = fct_rev(fct_relevel(mutation_type, "Background (no mutation)"))) %>%
  ggplot(mapping = aes(x=mutation_type, y = activity, fill = mutation_type)) +
  geom_boxplot(notch = F, alpha = 0.8, show.legend = F, color = "black", outlier.size = 1) +
  scale_fill_viridis(discrete = T) +
  coord_flip() +
  theme_classic()
distribution


# filter the mutations that occurred in TFs
# select only the pairs sample-TF that are unique to the same type of mutation
# this was done to prevent the assignment of the same TF activity to different types of mutations
tf_mutations <- cptac_mutations %>%
  semi_join(tf_activity, by = c("gene", "sample")) %>%
  group_by(mutation_type) %>%
  nest() %>%
  ungroup() %>%
  mutate(data2 = (function(x) map2(.x=1:length(x),.y=x,.f=function(a,b)reduce(.x=c(list(b),x[-a]),.f=setdiff))) (data)) %>%
  select(-data) %>%
  filter(map_lgl(data2, ~ nrow(.x) > 10 )) %>%
  unnest() %>%
  mutate(protein = "TFs")

# define a background to compare with (no known mutation for the TF in the sample)
tf_background <- tf_activity %>%
  anti_join(cptac_mutations, by = c("gene", "sample")) %>%
  mutate(mutation_type = "Background (no mutation)") %>%
  mutate(protein = "TFs")

tf_mutation_data <- tf_activity %>%
  inner_join(tf_mutations, by = c("sample", "gene")) %>%
  bind_rows(tf_background)

distribution <- tf_mutation_data %>%
  mutate(mutation_type = fct_rev(fct_reorder(mutation_type, activity, .fun = median))) %>%
  mutate(mutation_type = fct_rev(fct_relevel(mutation_type, "Background (no mutation)"))) %>%
  ggplot(mapping = aes(x=mutation_type, y = activity, fill = mutation_type)) +
  geom_boxplot(notch = F, alpha = 0.8, show.legend = F, color = "black", outlier.size = 1) +
  scale_fill_viridis(discrete = T) +
  coord_flip() +
  theme_classic()
distribution


# join kinase and TF activities and plot distributions
kinase_tf_mutation_data <- bind_rows(kinase_mutation_data, tf_mutation_data)

filter_activities <- function(vect, lims){
  x <- vect
  x <- x[x >= min(lims) & x <= max(lims)]
  
  return(x)
}

kinase_tf_mutation_data2 <- kinase_tf_mutation_data %>%
  group_by(protein, mutation_type) %>%
  summarise(activity = list(activity)) %>%
  ungroup() %>%
  mutate(limits = map(.x = activity, .f = ~ c(quantile(.x, 0.25)-(1.5*IQR(.x)), quantile(.x, 0.75)+(1.5*IQR(.x))))) %>%
  mutate(activity = map2(.x = activity, .y = limits, .f = filter_activities)) %>%
  select(-limits) %>%
  unnest()

pvalues_ <- kinase_tf_mutation_data2 %>%
  group_by(protein) %>%
  t_test(activity ~ mutation_type, ref.group = "Background (no mutation)", p.adjust.method = "none") %>%
  add_xy_position(x = "mutation_type")

pvalues <- kinase_tf_mutation_data %>%
  group_by(protein) %>%
  t_test(activity ~ mutation_type, ref.group = "Background (no mutation)", p.adjust.method = "none") %>%
  add_xy_position(x = "mutation_type") %>%
  inner_join(pvalues_[, c("protein", "group2", "y.position")], by = c("protein", "group2")) %>%
  rename(y.position = y.position.y)

N <- kinase_tf_mutation_data %>%
  group_by(protein, mutation_type) %>%
  summarise(n = n()) %>%
  ungroup()

N_ <- kinase_tf_mutation_data2 %>%
  group_by(protein, mutation_type) %>%
  summarise(min = min(activity)) %>%
  mutate(min = min(min)+0.5) %>%
  ungroup()

N <- N %>% inner_join(N_)

distribution <- kinase_tf_mutation_data2 %>%
  ggplot(mapping = aes(x=mutation_type, y = activity)) +
  geom_boxplot(mapping = aes(fill = mutation_type), alpha = 0.8, show.legend = F, color = "black", outlier.shape = NA) +
  stat_pvalue_manual(pvalues, label = "p", tip.length = 0.01, coord.flip = TRUE, size = 6, step.increase = 0.03) +
  geom_text(data = N, mapping = aes(x = rep(seq(1.3, 7.3, by = 1),2), y = min, label = n), size = 5) +
  facet_wrap(~ protein, scales = "free_x") +
  scale_fill_viridis(discrete = T) +
  coord_flip() +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, colour = "black", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 22, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 22, colour = "black"),
    axis.text.x = element_text(size = 20, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black")) +
  labs(y = "Activity score", title = "CPTAC tumours")
distribution

ggsave(filename = "cptac_kin_tf_act_mutation_distribution_box.png", plot = distribution, path = "./output/plots/mutation_impact_kin_tf_activity/", height = 5, width = 12)
ggsave(filename = "cptac_kin_tf_act_mutation_distribution_box.pdf", plot = distribution, path = "./output/plots/mutation_impact_kin_tf_activity/", height = 5, width = 12)


# load CPTAC tumours purity scores
purity_scores <- read_tsv(file = "./output/files/samples_purity.txt")
pure_samples <- purity_scores %>%
  filter(score > 0.8)

kinase_tf_mutation_data_pure <- kinase_tf_mutation_data %>%
  filter(sample %in% pure_samples$sample)

kinase_tf_mutation_data_pure2 <- kinase_tf_mutation_data_pure %>%
  group_by(protein, mutation_type) %>%
  summarise(activity = list(activity)) %>%
  ungroup() %>%
  mutate(limits = map(.x = activity, .f = ~ c(quantile(.x, 0.25)-(1.5*IQR(.x)), quantile(.x, 0.75)+(1.5*IQR(.x))))) %>%
  mutate(activity = map2(.x = activity, .y = limits, .f = filter_activities)) %>%
  select(-limits) %>%
  unnest()

pvalues_ <- kinase_tf_mutation_data_pure2 %>%
  group_by(protein) %>%
  t_test(activity ~ mutation_type, ref.group = "Background (no mutation)", p.adjust.method = "none") %>%
  add_xy_position(x = "mutation_type")

pvalues <- kinase_tf_mutation_data_pure %>%
  group_by(protein) %>%
  t_test(activity ~ mutation_type, ref.group = "Background (no mutation)", p.adjust.method = "none") %>%
  add_xy_position(x = "mutation_type") %>%
  inner_join(pvalues_[, c("protein", "group2", "y.position")], by = c("protein", "group2")) %>%
  rename(y.position = y.position.y)

N <- kinase_tf_mutation_data_pure %>%
  group_by(protein, mutation_type) %>%
  summarise(n = n()) %>%
  ungroup()

N_ <- kinase_tf_mutation_data_pure2 %>%
  group_by(protein, mutation_type) %>%
  summarise(min = min(activity)) %>%
  mutate(min = min(min)+0.5) %>%
  ungroup()

N <- N %>% inner_join(N_)

distribution <- kinase_tf_mutation_data_pure2 %>%
  ggplot(mapping = aes(x=mutation_type, y = activity)) +
  geom_boxplot(mapping = aes(fill = mutation_type), alpha = 0.8, show.legend = F, color = "black", outlier.shape = NA) +
  stat_pvalue_manual(pvalues, label = "p", tip.length = 0.01, coord.flip = TRUE, size = 6, step.increase = 0.03) +
  geom_text(data = N, mapping = aes(x = rep(seq(1.3, 7.3, by = 1),2), y = min, label = n), size = 5) +
  facet_wrap(~ protein, scales = "free_x") +
  scale_fill_viridis(discrete = T) +
  coord_flip() +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, colour = "black", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 18, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 16, colour = "black")) +
  labs(y = "Activity", title = "CPTAC tumours")
distribution

ggsave(filename = "cptac_kin_tf_act_mutation_distribution_box_pure_samples.png", plot = distribution, path = "./output/plots/mutation_impact_kin_tf_activity/", height = 5, width = 12)
ggsave(filename = "cptac_kin_tf_act_mutation_distribution_box_pure_samples.pdf", plot = distribution, path = "./output/plots/mutation_impact_kin_tf_activity/", height = 5, width = 12)



# -- for NCI60/CRC65 cell lines


# load NCI60/CRC65 cell lines
nci60_crc65 <- read_tsv(file = "./output/files/nci60_crc65_cell_lines.txt")

shared_cell_lines <- nci60_crc65 %>%
  group_by(batch) %>%
  summarise(cell_line = list(cell_line)) %>%
  pull(cell_line) %>%
  reduce(intersect)

alternative_names <- nci60_crc65 %>%
  filter(!is.na(alternative_names)) %>%
  mutate(alternative_names = str_split(alternative_names, ";")) %>%
  unnest() %>%
  distinct()


# load NCI60/CRC65 mutations
# select only the gene, sample and type of mutation
# remove duplicates (same gene in the same sample can have multiple mutations of same type)
cells_mutations <- read_tsv(file = "./output/files/nci60_crc65_mutations.txt.gz") %>%
  select(sample, gene, mutation_type) %>%
  distinct()


# load TF activity inference data from NCI60/CRC65 dataset
tf_activity_crc65 <- read_csv(file = "./data/dorothea/nci60_crc65/TF_activity_log2FC_allsamples_CRC65.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "activity") %>%
  rename(gene=X1) %>%
  filter(!sample %in% shared_cell_lines) %>%
  select(sample, gene, activity)

tf_activity_nci60 <- read_csv(file = "./data/dorothea/nci60_crc65/TF_activity_log2FC_allsamples_NCI60.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "activity") %>%
  rename(gene=X1) %>%
  select(sample, gene, activity)

tf_activity_cells <- bind_rows(tf_activity_crc65, tf_activity_nci60)


# load kinase activity inference data from NCI60/CRC65 dataset
# select quantifications with more than 3 substrates
kin_activity_cells <- read_tsv(file = "./output/files/nci60_crc65_kinase_activities.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= 3) %>%
  select(-source_type, -n) %>%
  filter(!(batch == "CRC65" & sample %in% shared_cell_lines)) %>%
  select(sample, gene = kinase, activity = log10P)


# filter the mutations that occurred in kinases
# select only the pairs sample-kinase that are unique to the same type of mutation
# this was done to prevent the assignment of the same kinase activity to different types of mutations
kin_mutations_cells <- cells_mutations %>%
  semi_join(kin_activity_cells, by = c("gene", "sample")) %>%
  group_by(mutation_type) %>%
  nest() %>%
  ungroup() %>%
  mutate(data2 = (function(x) map2(.x=1:length(x),.y=x,.f=function(a,b)reduce(.x=c(list(b),x[-a]),.f=setdiff))) (data)) %>%
  select(-data) %>%
  filter(map_lgl(data2, ~ nrow(.x) > 5 )) %>%
  unnest() %>%
  mutate(protein = "Kinases")

# define a background to compare with (no known mutation for the kinase in the sample)
kin_background_cells <- kin_activity_cells %>%
  anti_join(cells_mutations, by = c("gene", "sample")) %>%
  mutate(mutation_type = "Background (no mutation)") %>%
  mutate(protein = "Kinases")

kinase_mutation_data_cells <- kin_activity_cells %>%
  inner_join(kin_mutations_cells, by = c("sample", "gene")) %>%
  bind_rows(kin_background_cells)

distribution <- kinase_mutation_data_cells %>%
  mutate(mutation_type = fct_rev(fct_reorder(mutation_type, activity, .fun = median))) %>%
  mutate(mutation_type = fct_rev(fct_relevel(mutation_type, "Background (no mutation)"))) %>%
  ggplot(mapping = aes(x=mutation_type, y = activity, fill = mutation_type)) +
  geom_boxplot(notch = F, alpha = 0.8, show.legend = F, color = "black", outlier.size = 1) +
  scale_fill_viridis(discrete = T) +
  coord_flip() +
  theme_classic()
distribution


# filter the mutations that occurred in TFs
# select only the pairs sample-TF that are unique to the same type of mutation
# this was done to prevent the assignment of the same TF activity to different types of mutations
tf_mutations_cells <- cells_mutations %>%
  semi_join(tf_activity_cells, by = c("gene", "sample")) %>%
  group_by(mutation_type) %>%
  nest() %>%
  ungroup() %>%
  mutate(data2 = (function(x) map2(.x=1:length(x),.y=x,.f=function(a,b)reduce(.x=c(list(b),x[-a]),.f=setdiff))) (data)) %>%
  select(-data) %>%
  filter(map_lgl(data2, ~ nrow(.x) > 5 )) %>%
  unnest() %>%
  mutate(protein = "TFs")

# define a background to compare with (no known mutation for the TF in the sample)
tf_background_cells <- tf_activity_cells %>%
  anti_join(cells_mutations, by = c("gene", "sample")) %>%
  mutate(mutation_type = "Background (no mutation)") %>%
  mutate(protein = "TFs")

tf_mutation_data_cells <- tf_activity_cells %>%
  inner_join(tf_mutations_cells, by = c("sample", "gene")) %>%
  bind_rows(tf_background_cells)

distribution <- tf_mutation_data_cells %>%
  mutate(mutation_type = fct_rev(fct_reorder(mutation_type, activity, .fun = median))) %>%
  mutate(mutation_type = fct_rev(fct_relevel(mutation_type, "Background (no mutation)"))) %>%
  ggplot(mapping = aes(x=mutation_type, y = activity, fill = mutation_type)) +
  geom_boxplot(notch = F, alpha = 0.8, show.legend = F, color = "black", outlier.size = 1) +
  scale_fill_viridis(discrete = T) +
  coord_flip() +
  theme_classic()
distribution



# join kinase and TF activities and plot distributions
kinase_tf_mutation_data_cells <- bind_rows(kinase_mutation_data_cells, tf_mutation_data_cells)

filter_activities <- function(vect, lims){
  x <- vect
  x <- x[x >= min(lims) & x <= max(lims)]
  
  return(x)
}

kinase_tf_mutation_data_cells2 <- kinase_tf_mutation_data_cells %>%
  group_by(protein, mutation_type) %>%
  summarise(activity = list(activity)) %>%
  ungroup() %>%
  mutate(limits = map(.x = activity, .f = ~ c(quantile(.x, 0.25)-(1.5*IQR(.x)), quantile(.x, 0.75)+(1.5*IQR(.x))))) %>%
  mutate(activity = map2(.x = activity, .y = limits, .f = filter_activities)) %>%
  select(-limits) %>%
  unnest()

pvalues_ <- kinase_tf_mutation_data_cells2 %>%
  group_by(protein) %>%
  t_test(activity ~ mutation_type, ref.group = "Background (no mutation)", p.adjust.method = "none") %>%
  add_xy_position(x = "mutation_type")

pvalues <- kinase_tf_mutation_data_cells %>%
  group_by(protein) %>%
  t_test(activity ~ mutation_type, ref.group = "Background (no mutation)", p.adjust.method = "none") %>%
  add_xy_position(x = "mutation_type") %>%
  inner_join(pvalues_[, c("protein", "group2", "y.position")], by = c("protein", "group2")) %>%
  rename(y.position = y.position.y)

pvalues <- pvalues %>%
  mutate(xmax = c(2:5,2:7)) %>%
  mutate(y.position = map2_dbl(protein, y.position, ~ if(.x == "Kinases"){.y}else{.y-1}))

N <- kinase_tf_mutation_data_cells %>%
  group_by(protein, mutation_type) %>%
  summarise(n = n()) %>%
  ungroup()

N_ <- kinase_tf_mutation_data_cells2 %>%
  group_by(protein, mutation_type) %>%
  summarise(min = min(activity)) %>%
  mutate(min = min(min)) %>%
  ungroup() %>%
  mutate(min = map2_dbl(protein, min, ~ if(.x == "Kinases"){.y+0.5}else{.y+1}))

N <- N %>% inner_join(N_)

distribution <- kinase_tf_mutation_data_cells2 %>%
  ggplot(mapping = aes(x=mutation_type, y = activity)) +
  geom_boxplot(mapping = aes(fill = mutation_type), alpha = 0.8, show.legend = F, color = "black", outlier.shape = NA) +
  stat_pvalue_manual(pvalues, label = "p", tip.length = 0.01, coord.flip = TRUE, size = 5, step.increase = c(rep(0, 4), rep(0.04, 6))) +
  geom_text(data = N, mapping = aes(x = c(seq(1.3, 5.3, by = 1),seq(1.3, 7.3, by = 1)), y = min, label = n), size = 5) +
  facet_wrap(~ protein, scales = "free") +
  scale_fill_viridis(discrete = T) +
  coord_flip() +
  theme_classic() +
  theme(
    plot.title = element_text(size = 22, colour = "black", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 22, colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 22, colour = "black"),
    axis.text.x = element_text(size = 20, colour = "black"),
    axis.text.y = element_text(size = 16, colour = "black")) +
  labs(y = "Activity score", title = "NCI60/CRC65 cancer cell lines")
distribution

ggsave(filename = "nci60_crc65_kin_tf_act_mutation_distribution_box.png", plot = distribution, path = "./output/plots/mutation_impact_kin_tf_activity/", height = 5, width = 12)
ggsave(filename = "nci60_crc65_kin_tf_act_mutation_distribution_box.pdf", plot = distribution, path = "./output/plots/mutation_impact_kin_tf_activity/", height = 5, width = 12)
