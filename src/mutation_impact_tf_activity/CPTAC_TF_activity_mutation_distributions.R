#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)

#' load R packages
library(tidyverse)

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load CPTAC samples
source("./src/utils/getSamples.R")
cptac_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")
cptac_samples <- getSamples(cptac_samples, c("protein","mRNA")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = str_c(batch, tissue, sep = "-")) %>%
  select(-tissue)


#' load CPTAC TF activities
tf_activity <- read_csv(file = "./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "tf_activity") %>%
  rename(tf=X1) %>%
  select(sample, tf, tf_activity)


#' load CPTAC mutations\
#' select only the gene, sample and type of mutation\
#' remove duplicates (same gene in the same sample can have multiple mutations of same type)
cptac_mutations <- read_tsv(file = "./output/files/mutations_protpos.txt.gz") %>%
  select(sample, gene = gene_symbol, mutation_type = variant_class) %>%
  distinct()


#' filter the mutations that occurred in TFs\
#' select only the pairs sample-TF that are unique to the same type of mutation\
#' this was done to prevent the assignment of the same TF activity to different types of mutations
mutations <- cptac_mutations %>%
  semi_join(tf_activity, by = c("gene" = "tf", "sample")) %>%
  group_by(mutation_type) %>%
  nest() %>%
  ungroup() %>%
  mutate(data2 = (function(x) map2(.x=1:length(x),.y=x,.f=function(a,b)reduce(.x=c(list(b),x[-a]),.f=setdiff))) (data)) %>%
  select(-data) %>%
  unnest()


#' define a background to compare with (no known mutation for the gene in the sample)
background <- tf_activity %>%
  anti_join(cptac_mutations, by = c("tf" = "gene", "sample")) %>%
  mutate(mutation_type = "background (no mutation)")


#' # -- assess the TF activity distributions between different mutation types
distribution <- tf_activity %>%
  inner_join(mutations, by = c("sample", "tf" = "gene")) %>%
  ggplot(mapping = aes(x=mutation_type, y = tf_activity, fill = mutation_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, show.legend = F) +
  geom_jitter(size = 1, width = 0.1, alpha = 0.1, show.legend = F) +
  coord_flip()

#+ fig.width=10, fig.height=4
distribution


distribution <- tf_activity %>%
  inner_join(mutations, by = c("sample", "tf" = "gene")) %>%
  ggplot(mapping = aes(x=tf_activity, fill = mutation_type, color = mutation_type)) +
  geom_density(alpha = 0.5)

#+ fig.width=8, fig.height=3
distribution


#' # -- add the background
distribution <- tf_activity %>%
  inner_join(mutations, by = c("sample", "tf" = "gene")) %>%
  bind_rows(background) %>%
  ggplot(mapping = aes(x=mutation_type, y = tf_activity, fill = mutation_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, show.legend = F) +
  geom_jitter(size = 1, width = 0.1, alpha = 0.1, show.legend = F) +
  coord_flip()

#+ fig.width=10, fig.height=4
distribution


distribution <- tf_activity %>%
  inner_join(mutations, by = c("sample", "tf" = "gene")) %>%
  bind_rows(background) %>%
  ggplot(mapping = aes(x=tf_activity, fill = mutation_type, color = mutation_type)) +
  geom_density(alpha = 0.5)

#+ fig.width=8, fig.height=3
distribution
