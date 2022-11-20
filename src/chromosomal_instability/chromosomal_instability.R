#' ---
#' title: "Chromosomal instability"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=12, fig.height=6)


#' load R packages
library(ComplexHeatmap)
library(tidyverse)
library(lmtest)


source("./src/utils/getSamples.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load samples annotation
samples_annotation <- read_tsv(file = "./output/files/all_samples_annotation.txt")
metadata <- getSamples(samples_annotation, c("phosphorylation"))


#' load chromosomal instability scores
cin <- read_tsv(file = "./output/files/cin.txt.gz")


#' load purity scores
purity <- read_tsv(file = "./output/files/samples_purity.txt")


#' load kinase-activity inference data without imputations\
#' (quantile-normalized protein regressed-out phosphorylation data)\
#' select quantifications with more than 3 substrates
k_subN <- 3
ka <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  select(sample, kinase, activity=log10P)


#' load kinase activities and imputed values\
#' matrix used for PCA analysis and UMAP
kin_activities <- data.table::fread("./data/Danish/kinaseActMatImputed.tsv") %>%
  as_tibble() %>%
  rename(kinase = V1) %>%
  pivot_longer(-kinase, names_to = "sample", values_to = "activity")


#' load kinase-activity principal components\
#' select 1:10 principal components
pcs <- data.table::fread(file = "./data/Danish/kinasePCAMat.tsv") %>%
  as_tibble() %>%
  rename(sample = V1) %>%
  select(1:11) %>%
  pivot_longer(-sample, names_to = "pc", values_to = "value")


#' load UMAP projections
umap <- read_tsv("./output/files/kinActivities_umap_projections.txt") %>%
  pivot_longer(-sample, names_to = "umap", values_to = "value")


#' function to fit the linear associations and compare two nested models
reg1 <- function(df){
  dat <- df
  
  # fit the models
  if(length(unique(dat$batch)) == 1){
    m1 <- lm(activity ~ purity, data = dat)
    m2 <- lm(activity ~ purity + cin, data = dat)
  } else {
    m1 <- lm(activity ~ purity + batch, data = dat)
    m2 <- lm(activity ~ purity + batch + cin, data = dat)
  }
  
  # compare the two models using an ANOVA test
  anv <- anova(m1, m2)
  anv_pval <- anv[[6]][2]
  
  # compare the two models using a likelihood ratio test
  lrt <- lrtest(m2, m1)
  lrt_pval <- lrt[2,5]
  
  # compute and extract multiple statistical values
  m2 <- summary(m2)
  Ftest <- pf(m2$fstatistic[1], m2$fstatistic[2], m2$fstatistic[3], lower.tail = F)
  
  m2 <- broom::tidy(m2)
  purity_beta <- m2[m2$term == "purity", "estimate", drop=T]
  purity_pval <- m2[m2$term == "purity", "p.value", drop=T]
  
  cin_beta <- m2[m2$term == "cin", "estimate", drop=T]
  cin_pval <- m2[m2$term == "cin", "p.value", drop=T]
  
  
  res <- tibble(anv_pval=anv_pval, lrt_pval=lrt_pval, Ftest=Ftest, purity_beta=purity_beta, purity_pval=purity_pval, cin_beta=cin_beta, cin_pval=cin_pval)
  res
}


#' function to fit a simple linear model
reg2 <- function(df, covar){
  dat <- df
  
  # fit the models
  if(length(unique(dat$batch)) == 1){
    m1 <- lm(as.formula(str_c("activity ~ ", covar, sep = "")), data = dat)
  } else {
    m1 <- lm(as.formula(str_c("activity ~ batch + ", covar, sep = "")), data = dat)
  }
  
  # compute and extract multiple statistical values
  m1 <- summary(m1)
  Ftest <- pf(m1$fstatistic[1], m1$fstatistic[2], m1$fstatistic[3], lower.tail = F)
  
  m1 <- broom::tidy(m1)

  beta <- m1[m1$term == covar, "estimate", drop=T]
  pval <- m1[m1$term == covar, "p.value", drop=T]
  
  
  res <- tibble(Ftest=Ftest, beta=beta, pval=pval) %>%
    set_names(nm = c("Ftest", str_c(covar, "beta", sep = "_"), str_c(covar, "pval", sep = "_")))
  res
}


#' # correlate UMAP projections with chromosomal instability scores
correlations <- umap %>%
  inner_join(cin[, c("sample", "scaled")], by = "sample") %>%
  group_by(umap) %>%
  summarise(cor = cor(value, scaled, use = "p")) %>%
  ungroup()

correlations


#' # correlate principal components with chromosomal instability scores
correlations <- pcs %>%
  inner_join(cin[, c("sample", "scaled")], by = "sample") %>%
  group_by(pc) %>%
  summarise(cor = cor(value, scaled, use = "p")) %>%
  ungroup()

correlations


#' # correlate UMAP projections with purity scores
correlations <- umap %>%
  inner_join(purity[, c("sample", "score")], by = "sample") %>%
  group_by(umap) %>%
  summarise(cor = cor(value, score, use = "p")) %>%
  ungroup()

correlations


#' # correlate principal components with purity scores
correlations <- pcs %>%
  inner_join(purity[, c("sample", "score")], by = "sample") %>%
  group_by(pc) %>%
  summarise(cor = cor(value, score, use = "p")) %>%
  ungroup()

correlations


#' # estimate the association between kinase activity and CIN\
#' use a linear regression with purity and experimental batch as covariates\
#' use the kinase activities without imputations\
#' select kinases with at least 10 samples
associations <- ka %>%
  inner_join(cin[, c("sample", "scaled")], by = "sample") %>%
  rename(cin=scaled) %>%
  inner_join(purity[, c("sample", "score")], by = "sample") %>%
  rename(purity=score) %>%
  inner_join(metadata[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, ~ nrow(.x) > 10)) %>%
  mutate(res = map(data, reg1)) %>%
  select(-data) %>%
  unnest(res) %>%
  mutate(padj = p.adjust(anv_pval, method = "BH"))


#' significant associations (FDR < 0.05)
associations <- associations %>%
  filter(padj < 0.05)

associations


#' # estimate the association between kinase activity and CIN (without purity as a covariate)\
#' use a linear regression with the experimental batch as covariate\
#' use the kinase activities without imputations\
#' select kinases with at least 10 samples
associations <- ka %>%
  inner_join(cin[, c("sample", "scaled")], by = "sample") %>%
  rename(cin=scaled) %>%
  inner_join(purity[, c("sample", "score")], by = "sample") %>%
  rename(purity=score) %>%
  inner_join(metadata[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, ~ nrow(.x) > 10)) %>%
  mutate(res = map(data, reg2, covar = "cin")) %>%
  select(-data) %>%
  unnest(res) %>%
  mutate(padj = p.adjust(cin_pval, method = "BH"))


#' significant associations (FDR < 0.05)
associations <- associations %>%
  filter(padj < 0.05)

associations


#' # estimate the association between kinase activity and purity\
#' use a linear regression with the experimental batch as covariate\
#' use the kinase activities without imputations\
#' select kinases with at least 10 samples
associations <- ka %>%
  inner_join(cin[, c("sample", "scaled")], by = "sample") %>%
  rename(cin=scaled) %>%
  inner_join(purity[, c("sample", "score")], by = "sample") %>%
  rename(purity=score) %>%
  inner_join(metadata[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, ~ nrow(.x) > 10)) %>%
  mutate(res = map(data, reg2, covar = "purity")) %>%
  select(-data) %>%
  unnest(res) %>%
  mutate(padj = p.adjust(purity_pval, method = "BH"))


#' significant associations (FDR < 0.05)
associations <- associations %>%
  filter(padj < 0.05) %>%
  arrange(padj)

associations

kinases1 <- associations %>%
  pull(kinase)


#' # estimate the association between kinase activity and CIN\
#' use a linear regression with purity and experimental batch as covariates\
#' use the kinase activities with imputations\
associations <- kin_activities %>%
  inner_join(cin[, c("sample", "scaled")], by = "sample") %>%
  rename(cin=scaled) %>%
  inner_join(purity[, c("sample", "score")], by = "sample") %>%
  rename(purity=score) %>%
  inner_join(metadata[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(res = map(data, reg1)) %>%
  select(-data) %>%
  unnest(res) %>%
  mutate(padj = p.adjust(cin_pval, method = "BH"))


#' significant associations (FDR < 0.05)
#associations <- associations %>%
#  filter(padj < 0.05)

#associations


plot <- associations %>%
  filter(cin_pval < 0.05) %>%
  select(kinase) %>% 
  inner_join(kin_activities, by = "kinase") %>%
  inner_join(cin[, c("sample", "scaled")], by = "sample") %>%
  rename(cin=scaled) %>%
  inner_join(purity[, c("sample", "score")], by = "sample") %>%
  rename(purity=score) %>%
  inner_join(metadata[, c("sample", "batch")], by = "sample") %>%
  filter(!is.na(purity)) %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(activity_resd = map(.x = data, .f = ~ residuals(lm(activity ~ batch + purity, data = .x)))) %>%
  unnest() %>%
  ggplot(mapping = aes(x = cin, y = activity_resd)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor() +
  facet_wrap(~ kinase)

#+ fig.width=8, fig.height=6
plot


#' # estimate the association between kinase activity and CIN (without purity as a covariate)\
#' use a linear regression with the experimental batch as covariate\
#' use the kinase activities with imputations\
associations <- kin_activities %>%
  inner_join(cin[, c("sample", "scaled")], by = "sample") %>%
  rename(cin=scaled) %>%
  inner_join(purity[, c("sample", "score")], by = "sample") %>%
  rename(purity=score) %>%
  inner_join(metadata[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(res = map(data, reg2, covar = "cin")) %>%
  select(-data) %>%
  unnest(res) %>%
  mutate(padj = p.adjust(cin_pval, method = "BH"))


#' significant associations (FDR < 0.05)
associations <- associations %>%
  filter(padj < 0.05)

associations


#' # estimate the association between kinase activity and purity\
#' use a linear regression with the experimental batch as covariate\
#' use the kinase activities with imputations\
associations <- kin_activities %>%
  inner_join(cin[, c("sample", "scaled")], by = "sample") %>%
  rename(cin=scaled) %>%
  inner_join(purity[, c("sample", "score")], by = "sample") %>%
  rename(purity=score) %>%
  inner_join(metadata[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(res = map(data, reg2, covar = "purity")) %>%
  select(-data) %>%
  unnest(res) %>%
  mutate(padj = p.adjust(purity_pval, method = "BH"))


#' significant associations (FDR < 0.05)
associations <- associations %>%
  filter(padj < 0.05) %>%
  arrange(padj)

associations


kinases2 <- associations %>%
  pull(kinase)

length(intersect(kinases1, kinases2))


#' # plot some examples
plot <- associations %>%
  inner_join(kin_activities, by = "kinase") %>%
  inner_join(cin[, c("sample", "scaled")], by = "sample") %>%
  rename(cin=scaled) %>%
  inner_join(purity[, c("sample", "score")], by = "sample") %>%
  rename(purity=score) %>%
  inner_join(metadata[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(activity_resd = map(.x = data, .f = ~ residuals(lm(activity ~ batch, data = .x)))) %>%
  unnest() %>%
  ggplot(mapping = aes(x = purity, y = activity_resd)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", lwd = 0.5) +
  ggpubr::stat_cor(size = 2) +
  facet_wrap(~ kinase, scales = "free")

#+ fig.width=8, fig.height=6
plot
