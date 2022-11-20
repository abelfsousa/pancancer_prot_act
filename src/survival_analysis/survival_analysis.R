#' ---
#' title: "Survival analysis"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=4, fig.height=6)


#' load R packages
library(tidyverse)
library(survminer)
library(survival)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
#nest <- nest_legacy
#unnest <- unnest_legacy


#' load clinical data
clinical_data <- read_tsv(file = "./output/files/clinical_data.txt")


#' select tumour samples to analyze
samples <- clinical_data %>%
  filter(!batch %in% c("ccle", "colon-oppt")) %>%
  select(sample, cancer = tissue, OS_time, OS) %>%
  filter(!(is.na(OS) | is.na(OS_time)))


#' load kinase-activity inference data\
#' (quantile-normalized protein regressed-out phosphorylation data)\
#' select quantifications with more than 3 substrates
k_subN <- 3
ka <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  select(kinase, sample, log10P)


#' load kinase-activity principal components
pcs <- data.table::fread(file = "./data/Danish/kinasePCAMat.tsv") %>%
  as_tibble() %>%
  rename(sample = V1) %>%
  select(1:11) %>%
  pivot_longer(-sample, names_to = "pc", values_to = "value")


#' load UMAP projections
umap <- read_tsv("./output/files/kinActivities_umap_projections.txt")


#' function to classify the activities into neutral, activated and deactivated
act <- function(p, cutoff){
  if(p < -cutoff){activity <- "deactivated"}
  else if(p > cutoff){activity <- "activated"}
  else{activity <- "neutral"}
  activity
}


#' function to select kinases to test
selKinases <- function(df, ct_act, ct_os = 0, logic_act = "or"){
  tab <- table(df$activity)
  act = tab["activated"]
  deact = tab["deactivated"]
  if(is.na(act)) act = 0
  if(is.na(deact)) deact = 0
  
  tab <- table(df$OS)
  alive = tab["0"]
  dead = tab["1"]
  if(is.na(alive)) alive = 0
  if(is.na(dead)) dead = 0
  
  if(logic_act == "or") res <- (act > ct_act | deact > ct_act) & (dead > ct_os)
  if(logic_act == "and") res <- (act > ct_act & deact > ct_act) & (dead > ct_os)
  
  res
}


#' function to compare Kaplan-Meier survival curves
compare_KMCurves <- function(df){
  dat <- df

  km_test <- survdiff(Surv(OS_time, OS) ~ activity, data=dat)
  chq <- km_test$chisq
  df <- length(km_test$n)-1
  pval <- pchisq(chq, df, lower.tail = F)
  
  pval
}


#' function to plot Kaplan-Meier curves
plot_KMCurves <- function(kinase, cancer, df){
  dat <- df
  
  km_fit <- survfit(Surv(OS_time, OS) ~ activity, data=dat)
  plot <- ggsurvplot(km_fit, data = dat, pval = T, risk.table = T) +
    labs(title = str_c(kinase, cancer, sep = " - "))
  plot
}


#' function to fit the cox proportional hazards regression
cox_reg <- function(df, adjust = FALSE){
  dat <- df
  
  # compute cox hazard regression
  if(adjust){
    cox_test <- coxph(Surv(OS_time, OS) ~ cancer + covar, data=dat)
  } else {
    cox_test <- coxph(Surv(OS_time, OS) ~ covar, data=dat)
  }
  
  # statistics summary
  sm <- summary(cox_test)
  
  act_pval <- sm$coefficients["covar", "Pr(>|z|)"]
  act_beta <- sm$coefficients["covar", "coef"]
  act_hr <- sm$coefficients["covar", "exp(coef)"]
  hrCI95_lower <- sm$conf.int["covar", "lower .95"]
  hrCI95_upper <- sm$conf.int["covar", "upper .95"]
  lrt_pval <- sm$logtest["pvalue"]
  
  # test the proportional hazards assumption
  # independence of Schoenfeld residuals over time
  test_ph <- cox.zph(cox_test)
  ph_act_pval <- test_ph$table["covar", "p"]
  ph_global_pval <- test_ph$table["GLOBAL", "p"]
  
  
  tibble(ph_act_pval=ph_act_pval, ph_global_pval=ph_global_pval, act_pval=act_pval, act_beta=act_beta, act_hr=act_hr, hrCI95_lower=hrCI95_lower, hrCI95_upper=hrCI95_upper,  lrt_pval=lrt_pval)
}




#' ## (A) survival analysis using Kaplan-Meier survival estimates - by kinase and cancer type\
#' ### classify the activities into neutral, activated and deactivated (cutoff = 1.75)\
#' select kinases with at least 5 events of activation or deactivation (214/1115 tests)\
#' adjust P-values by cancer type
surv <- ka %>%
  mutate(activity = map_chr(log10P, act, cutoff = 1.75)) %>%
  inner_join(samples, by = "sample") %>%
  group_by(cancer, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, selKinases, ct_act = 5)) %>%
  mutate(pval = map_dbl(data, compare_KMCurves)) %>%
  group_by(cancer) %>%
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  arrange(padj)

#' associations with FDR < 0.1 (log-rank test P-value)
surv <- surv %>%
  filter(padj < 0.1)
surv

#' plot Kaplan-Meier curves for the associations with FDR < 10%
#pwalk(.l = surv, .f = ~ print(plot_KMCurves(..1, ..2, ..3)))
surv <- surv %>%
  mutate(KM_plots = pmap(.l = ., .f = ~ plot_KMCurves(..1, ..2, ..3)))

#+ fig.width=4, fig.height=6
surv$KM_plots


#' ## (B) survival analysis using Kaplan-Meier survival estimates - by kinase and cancer type\
#' ### classify the activities into neutral, activated and deactivated (cutoff = 1.5)\
#' select kinases with at least 5 events of activation or deactivation (308/1115 tests)\
#' adjust P-values by cancer type
surv <- ka %>%
  mutate(activity = map_chr(log10P, act, cutoff = 1.5)) %>%
  inner_join(samples, by = "sample") %>%
  group_by(cancer, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, selKinases, ct_act = 5)) %>%
  mutate(pval = map_dbl(data, compare_KMCurves)) %>%
  group_by(cancer) %>%
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  arrange(padj)

#' associations with FDR < 0.1 (log-rank test P-value)
surv <- surv %>%
  filter(padj < 0.1)
surv

#' plot Kaplan-Meier curves for the associations with FDR < 10%
#pwalk(.l = surv, .f = ~ print(plot_KMCurves(..1, ..2, ..3)))
surv <- surv %>%
  mutate(KM_plots = pmap(.l = ., .f = ~ plot_KMCurves(..1, ..2, ..3)))

#+ fig.width=4, fig.height=6
surv$KM_plots


#' ## (C) survival analysis using Kaplan-Meier survival estimates - by kinase and cancer type\
#' ### classify the activities into neutral, activated and deactivated (cutoff = 1)\
#' select kinases with at least 5 events of activation or deactivation (623/1115 tests)\
#' adjust P-values by cancer type
surv <- ka %>%
  mutate(activity = map_chr(log10P, act, cutoff = 1)) %>%
  inner_join(samples, by = "sample") %>%
  group_by(cancer, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, selKinases, ct_act = 5)) %>%
  mutate(pval = map_dbl(data, compare_KMCurves)) %>%
  group_by(cancer) %>%
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  ungroup() %>%
  arrange(padj)

#' associations with FDR < 0.1 (log-rank test P-value)
surv <- surv %>%
  filter(padj < 0.1)
surv

#' plot Kaplan-Meier curves for the associations with FDR < 0.1
#pwalk(.l = surv, .f = ~ print(plot_KMCurves(..1, ..2, ..3)))
surv <- surv %>%
  mutate(KM_plots = pmap(.l = ., .f = ~ plot_KMCurves(..1, ..2, ..3)))

#+ fig.width=4, fig.height=6
surv$KM_plots


#' ## (D) survival analysis by *cox proportional hazards* regression analysis - by kinase and cancer type\
#' ### classify the activities into neutral, activated and deactivated (cutoff = 1.75)\
#' select kinases with at least 5 events of activation or deactivation (214/1115 tests)\
#' adjust P-values by cancer type
surv <- ka %>%
  mutate(activity = map_chr(log10P, act, cutoff = 1.75)) %>%
  inner_join(samples, by = "sample") %>%
  rename(covar = log10P) %>%
  group_by(cancer, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, selKinases, ct_act = 5, logic = "or")) %>%
  mutate(res = map(data, cox_reg)) %>%
  unnest(res) %>%
  group_by(cancer) %>%
  mutate(padj = p.adjust(lrt_pval, method = "BH")) %>%
  ungroup() %>%
  arrange(padj)

#' associations with FDR < 0.1 (global statistical significance [likelihood-ratio test])
surv <- surv %>%
  filter(padj < 0.1)
as.data.frame(surv[, -c(3)])

#' plot Kaplan-Meier curves for the associations with FDR < 0.1
#pwalk(.l = surv, .f = ~ print(plot_KMCurves(..1, ..2, ..3)))
surv <- surv %>%
  mutate(KM_plots = pmap(.l = ., .f = ~ plot_KMCurves(..1, ..2, ..3)))

#+ fig.width=4, fig.height=6
surv$KM_plots


#' ## (E) survival analysis by *cox proportional hazards* regression analysis - by kinase and cancer type\
#' ### classify the activities into neutral, activated and deactivated (cutoff = 1.5)\
#' select kinases with at least 5 events of activation or deactivation (308/1115 tests)\
#' adjust P-values by cancer type
surv <- ka %>%
  mutate(activity = map_chr(log10P, act, cutoff = 1.5)) %>%
  inner_join(samples, by = "sample") %>%
  rename(covar = log10P) %>%
  group_by(cancer, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, selKinases, ct_act = 5, logic = "or")) %>%
  mutate(res = map(data, cox_reg)) %>%
  unnest(res) %>%
  group_by(cancer) %>%
  mutate(padj = p.adjust(lrt_pval, method = "BH")) %>%
  ungroup() %>%
  arrange(padj)

#' associations with FDR < 0.1 (global statistical significance [likelihood-ratio test])
surv <- surv %>%
  filter(padj < 0.1)
as.data.frame(surv[, -c(3)])

#' plot Kaplan-Meier curves for the associations with FDR < 0.1
#pwalk(.l = surv, .f = ~ print(plot_KMCurves(..1, ..2, ..3)))
surv <- surv %>%
  mutate(KM_plots = pmap(.l = ., .f = ~ plot_KMCurves(..1, ..2, ..3)))

#+ fig.width=4, fig.height=6
surv$KM_plots


#' ## (F) survival analysis by *cox proportional hazards* regression analysis - by kinase and cancer type\
#' ### classify the activities into neutral, activated and deactivated (cutoff = 1)\
#' select kinases with at least 5 events of activation or deactivation (623/1115 tests)\
#' adjust P-values by cancer type
surv <- ka %>%
  mutate(activity = map_chr(log10P, act, cutoff = 1)) %>%
  inner_join(samples, by = "sample") %>%
  rename(covar = log10P) %>%
  group_by(cancer, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, selKinases, ct_act = 5, logic = "or")) %>%
  mutate(res = map(data, cox_reg)) %>%
  unnest(res) %>%
  group_by(cancer) %>%
  mutate(padj = p.adjust(lrt_pval, method = "BH")) %>%
  ungroup() %>%
  arrange(padj)

#' associations with FDR < 0.1 (global statistical significance [likelihood-ratio test])
surv <- surv %>%
  filter(padj < 0.1)
as.data.frame(surv[, -c(3)])

#' plot Kaplan-Meier curves for the associations with FDR < 0.1
#pwalk(.l = surv, .f = ~ print(plot_KMCurves(..1, ..2, ..3)))
surv <- surv %>%
  mutate(KM_plots = pmap(.l = ., .f = ~ plot_KMCurves(..1, ..2, ..3)))

#+ fig.width=4, fig.height=6
surv$KM_plots


#' ## (G) survival analysis by *cox proportional hazards* regression analysis - by kinase\
#' ### classify the activities into neutral, activated and deactivated (cutoff = 1.75)\
#' select kinases with at least 10 events of activation or deactivation (94/208 tests)\
#' cancer type as covariate\
surv <- ka %>%
  mutate(activity = map_chr(log10P, act, cutoff = 1.75)) %>%
  inner_join(samples, by = "sample") %>%
  rename(covar = log10P) %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, selKinases, ct_act = 10)) %>%
  mutate(res = map(data, cox_reg, adjust = TRUE)) %>%
  unnest(res) %>%
  mutate(padj = p.adjust(lrt_pval, method = "BH")) %>%
  arrange(padj)

#' associations with FDR < 0.05 (global statistical significance [likelihood-ratio test]) and activity beta P-value < 0.05
surv <- surv %>%
  filter(padj < 0.05 & act_pval < 0.05)
as.data.frame(surv[, -c(2)])

#' plot Kaplan-Meier curves for the associations with FDR < 0.05 and activity beta < 0.05
#pwalk(.l = surv, .f = ~ print(plot_KMCurves(..1, ..2, ..3)))
surv <- surv %>%
  mutate(KM_plots = pmap(.l = ., .f = ~ plot_KMCurves(..1, "", ..2)))

#+ fig.width=6, fig.height=6
surv$KM_plots


#' ## (H) survival analysis by *cox proportional hazards* regression analysis - by kinase\
#' ### classify the activities into neutral, activated and deactivated (cutoff = 1.5)\
#' select kinases with at least 10 events of activation or deactivation (114/208 tests)\
#' cancer type as covariate\
surv <- ka %>%
  mutate(activity = map_chr(log10P, act, cutoff = 1.5)) %>%
  inner_join(samples, by = "sample") %>%
  rename(covar = log10P) %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, selKinases, ct_act = 10)) %>%
  mutate(res = map(data, cox_reg, adjust = TRUE)) %>%
  unnest(res) %>%
  mutate(padj = p.adjust(lrt_pval, method = "BH")) %>%
  arrange(padj)

#' associations with FDR < 0.05 (global statistical significance [likelihood-ratio test]) and activity beta P-value < 0.05
surv <- surv %>%
  filter(padj < 0.05 & act_pval < 0.05)
as.data.frame(surv[, -c(2)])

#' plot Kaplan-Meier curves for the associations with FDR < 0.05 and activity beta < 0.05
#pwalk(.l = surv, .f = ~ print(plot_KMCurves(..1, ..2, ..3)))
surv <- surv %>%
  mutate(KM_plots = pmap(.l = ., .f = ~ plot_KMCurves(..1, "", ..2)))

#+ fig.width=6, fig.height=6
surv$KM_plots


#' ## (I) survival analysis by *cox proportional hazards* regression analysis - by kinase\
#' ### classify the activities into neutral, activated and deactivated (cutoff = 1)\
#' select kinases with at least 10 events of activation or deactivation (143/208 tests)\
#' cancer type as covariate\
surv <- ka %>%
  mutate(activity = map_chr(log10P, act, cutoff = 1)) %>%
  inner_join(samples, by = "sample") %>%
  rename(covar = log10P) %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, selKinases, ct_act = 10)) %>%
  slice(-126) %>%
  mutate(res = map(data, cox_reg, adjust = TRUE)) %>%
  unnest(res) %>%
  mutate(padj = p.adjust(lrt_pval, method = "BH")) %>%
  arrange(padj)

#' associations with FDR < 0.05 (global statistical significance [likelihood-ratio test]) and activity beta P-value < 0.05
surv <- surv %>%
  filter(padj < 0.05 & act_pval < 0.05)
as.data.frame(surv[, -c(2)])

#' plot Kaplan-Meier curves for the associations with FDR < 0.05 and activity beta < 0.05
#pwalk(.l = surv, .f = ~ print(plot_KMCurves(..1, ..2, ..3)))
surv <- surv %>%
  mutate(KM_plots = pmap(.l = ., .f = ~ plot_KMCurves(..1, "", ..2)))

#+ fig.width=6, fig.height=6
surv$KM_plots


#' ## (J) survival analysis by *cox proportional hazards* regression analysis - by PCs
#' cancer type as covariate\
surv <- pcs %>%
  inner_join(samples, by = "sample") %>%
  rename(covar = value) %>%
  group_by(pc) %>%
  nest() %>%
  ungroup() %>%
  mutate(res = map(data, cox_reg, adjust = TRUE)) %>%
  unnest(res) %>%
  mutate(padj = p.adjust(lrt_pval, method = "BH")) %>%
  arrange(padj)


#' associations with global statistical significance (likelihood-ratio test) < 0.05 and PC beta P-value < 0.05
surv <- surv %>%
  filter(lrt_pval < 0.05 & act_pval < 0.05)
as.data.frame(surv[, -c(2)])


#cox_test <- coxph(Surv(OS_time, OS) ~ 1, data=surv$data[[1]])
#cox_test <- surv_fit(cox_test, surv$data[[1]])
#ggsurvplot(fit=cox_test, data=surv$data[[1]], ggtheme = theme_minimal())

#cox_test <- coxph(Surv(OS_time, OS) ~ cancer + covar, data=surv$data[[1]])
#cox_test <- surv_fit(cox_test, surv$data[[1]])
#ggsurvplot(fit=cox_test, data=surv$data[[1]], ggtheme = theme_minimal())


#' ## (K) survival analysis by *cox proportional hazards* regression analysis - by UMAP projections\
#' cancer type as covariate\
surv <- umap %>%
  pivot_longer(-sample, names_to = "projection", values_to = "value") %>%
  inner_join(samples, by = "sample") %>%
  rename(covar = value) %>%
  group_by(projection) %>%
  nest() %>%
  ungroup() %>%
  mutate(res = map(data, cox_reg, adjust = TRUE)) %>%
  unnest(res) %>%
  arrange(act_pval)

#' associations with global statistical significance (likelihood-ratio test) < 0.05 and UMAP beta P-value < 0.05
surv <- surv %>%
  filter(lrt_pval < 0.05 & act_pval < 0.05)
as.data.frame(surv[, -c(2)])



#' ## (L) survival analysis using Kaplan-Meier survival estimates - by kinase and cancer type
#' ### classify the activities into neutral, activated and deactivated (cutoff = 1.75)\
#' select kinases with at least 5 events of activation or deactivation (125/208 tests)\
compare_KMCurves <- function(df){
  dat <- df
  
  km_test <- survdiff(Surv(OS_time, OS) ~ cancer + activity, data=dat)
  chq <- km_test$chisq
  df <- length(km_test$n)-1
  pval <- pchisq(chq, df, lower.tail = F)
  
  pval
}

surv <- ka %>%
  mutate(activity = map_chr(log10P, act, cutoff = 1.75)) %>%
  inner_join(samples, by = "sample") %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, selKinases, ct_act = 5)) %>%
  mutate(pval = map_dbl(data, compare_KMCurves)) %>%
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  arrange(padj)

#' associations with FDR < 0.01 (log-rank test P-value)
surv <- surv %>%
  filter(padj < 0.01)
surv

km_fit <- survfit(Surv(OS_time, OS) ~ cancer + activity, data=surv$data[[1]])
plot <- ggsurvplot(km_fit, data = surv$data[[1]], pval = T, risk.table = T)

#+ fig.width=8, fig.height=8
plot$plot + theme(legend.position = "right")

