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
samples_annotation <- getSamples(samples_annotation, c("phosphorylation"))


# load transcription factor activities
#tfs <- data.table::fread("./data/progeny/TF_activity_log2FC.csv") %>%
tfs <- data.table::fread("./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  as_tibble() %>%
  rename(tf=V1) %>%
  pivot_longer(-tf, names_to = "sample", values_to = "tf_activity")


# load kinase-activity inference data without imputations
# (quantile-normalized protein regressed-out phosphorylation data)
# select quantifications with more than 3 substrates
k_subN <- 3
kin_activities <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  select(sample, kinase, kin_activity=log10P)


# load kinase activities with imputed values
# matrix used for PCA analysis and UMAP
kin_activities_imputed <- data.table::fread("./data/Danish/kinaseActMatImputed.tsv") %>%
  as_tibble() %>%
  rename(kinase = V1) %>%
  pivot_longer(-kinase, names_to = "sample", values_to = "kin_activity")


# function to fit a multiple linear model
reg <- function(df, y, x){
  dat <- df
  
  # fit the models
  if(length(unique(dat$batch)) == 1){
    m1 <- lm(as.formula(str_c(y, "~", x)), data = dat)
  } else {
    m1 <- lm(as.formula(str_c(y, "~batch+", x)), data = dat)
  }
  
  # compute and extract multiple statistical values
  m1 <- broom::tidy(m1)
  
  beta <- m1[m1$term == x, "estimate", drop=T]
  pval <- m1[m1$term == x, "p.value", drop=T]
  
  
  res <- tibble(beta=beta, pval=pval)
  
  res
}


# function to compare two nested linear models using a LRT
lr_test <- function(df, covars){
  dat <- df %>%
    inner_join(covars, by = "sample")
  
  # fit the models
  m1 <- lm(tf_activity ~ batch, data = dat)
  m2 <- lm(tf_activity ~ batch + kin_activity, data = dat)
  
  # compare the two models using a likelihood ratio test
  lrt <- lrtest(m2, m1)
  lrt_pval <- lrt[2,5]
  
  # compute and extract multiple statistical values
  m2 <- broom::tidy(m2)
  
  beta <- m2[m2$term == "kin_activity", "estimate", drop=T]
  pval <- m2[m2$term == "kin_activity", "p.value", drop=T]
  
  
  res <- tibble(kin_beta=beta, kin_pval=pval, lrt_pval = lrt_pval)
  
  res
}


# associate kinase and TF activities using linear models
# use the kinase activities without imputations
# select kinase-TF pairs with more than 10 samples
# use experimental batch as covariate
# TF activity as dependent and kinase activity as independent variable
associations <- kin_activities %>%
  inner_join(tfs, by = "sample") %>%
  inner_join(samples_annotation[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase, tf) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, ~ nrow(.x) > 10)) %>%
  mutate(models = map(.x=data, .f = reg, y = "tf_activity", x = "kin_activity")) %>%
  select(-data) %>%
  unnest() %>%
  rename(kin_beta=beta, kin_pval=pval) %>%
  mutate(padj = p.adjust(kin_pval, method = "BH"))
 write_tsv(associations, "./output/files/kinaseNotImputed_TF_activity_associations.txt.gz")


# associate kinase and TF activities using linear models
# use the kinase activities with imputations (matrix used for PCA and UMAP)
# use experimental batch as covariate
# TF activity as dependent and kinase activity as independent variable
associations2 <- kin_activities_imputed %>%
  inner_join(tfs, by = "sample") %>%
  inner_join(samples_annotation[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase, tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(models = map(.x=data, .f = reg, y = "tf_activity", x = "kin_activity")) %>%
  select(-data) %>%
  unnest() %>%
  rename(kin_beta=beta, kin_pval=pval) %>%
  mutate(padj = p.adjust(kin_pval, method = "BH"))
write_tsv(associations2, "./output/files/kinaseImputed_TF_activity_associations.txt.gz")

# kinase activity as dependent and tf activity as independent variable
associations3 <- kin_activities_imputed %>%
  inner_join(tfs, by = "sample") %>%
  inner_join(samples_annotation[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase, tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(models = map(.x=data, .f = reg, y = "kin_activity", x = "tf_activity")) %>%
  select(-data) %>%
  unnest() %>%
  rename(tf_beta=beta, tf_pval=pval) %>%
  mutate(padj = p.adjust(tf_pval, method = "BH"))
write_tsv(associations3, "./output/files/kinaseImputed_TF_activity_associations2.txt.gz")


# associate kinase and TF activities using two nested linear models
# compare them using a likelihood ratio test
# use the kinase activities with imputations (matrix used for PCA and UMAP)
# use experimental batch as covariate
associations4 <- kin_activities_imputed %>%
  inner_join(tfs, by = "sample") %>%
  group_by(kinase, tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(models = map(.x=data, .f = lr_test, covars = samples_annotation[, c("sample", "batch")])) %>%
  select(-data) %>%
  unnest() %>%
  mutate(padj = p.adjust(lrt_pval, method = "BH"))
write_tsv(associations4, "./output/files/kinaseImputed_TF_activity_associations_lrtest.txt.gz")
