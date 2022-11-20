# load R packages
library(tidyverse)
library(lmtest)

source("./src/utils/getSamples.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load samples annotation
samples_annotation <- read_tsv(file = "./output/files/all_samples_annotation.txt")
samples_annotation <- getSamples(samples_annotation, c("phosphorylation", "mRNA"))


# load pathway activities
#pathways <- data.table::fread("./data/progeny/pathway_activities_log2fpkm.csv") %>%
pathways <- data.table::fread("./data/progeny/pathway_prog_log2FC.csv") %>%
  as_tibble() %>%
  rename(sample=V1) %>%
  pivot_longer(-sample, names_to = "pathway", values_to = "pathway_activity") %>%
  mutate(sample = str_replace(sample, "^X{1}", "")) %>%
  mutate(pathway = str_replace(pathway, "JAK.STAT", "JAK-STAT"))


# load kinase activities and imputed values
# matrix used for PCA analysis and UMAP
kin_activities <- data.table::fread("./data/Danish/kinaseActMatImputed.tsv") %>%
  as_tibble() %>%
  rename(kinase = V1) %>%
  pivot_longer(-kinase, names_to = "sample", values_to = "kin_activity")


# function to fit a multiple linear model
reg <- function(df, covars){
  dat <- df %>%
    inner_join(covars, by = "sample")
  
  # fit the models
  m1 <- lm(pathway_activity ~ batch + kin_activity, data = dat)
  
  # compute and extract multiple statistical values
  m1 <- broom::tidy(m1)
  
  beta <- m1[m1$term == "kin_activity", "estimate", drop=T]
  pval <- m1[m1$term == "kin_activity", "p.value", drop=T]
  
  res <- tibble(kin_beta=beta, kin_pval=pval)
  
  res
}


# associate kinase and pathway activities using linear models
# use the kinase activities with imputations (matrix used for PCA and UMAP)
# use experimental batch as covariate
associations <- kin_activities %>%
  inner_join(pathways, by = "sample") %>%
  group_by(kinase, pathway) %>%
  nest() %>%
  ungroup() %>%
  mutate(models = map(.x=data, .f = reg, covars = samples_annotation[, c("sample", "batch")])) %>%
  select(-data) %>%
  unnest() %>%
  mutate(padj = p.adjust(kin_pval, method = "BH"))

write_tsv(associations, "./output/files/kinaseImputed_pathway_activity_associations.txt.gz")
