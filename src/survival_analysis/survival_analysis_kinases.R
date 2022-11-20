# Survival analysis using kinases

# load R packages
library(tidyverse)
library(survminer)
library(survival)
library(ComplexHeatmap)
library(RColorBrewer)



# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
#nest <- nest_legacy
#unnest <- unnest_legacy


# load clinical data
clinical_data <- read_tsv(file = "./output/files/clinical_data.txt")


# select tumour samples to analyze
samples <- clinical_data %>%
  filter(!batch %in% c("ccle", "colon-oppt")) %>%
  select(sample, cancer = tissue, OS_time, OS) %>%
  filter(!(is.na(OS) | is.na(OS_time)))


# select tissues with at least 10 deaths?
# samples <- samples %>%
#   group_by(cancer) %>%
#   nest() %>%
#   ungroup() %>%
#   mutate(deaths = map_dbl(data, ~ .x %>% filter(OS == 1) %>% nrow())) %>%
#   filter(deaths > 10) %>%
#   select(-deaths) %>%
#   unnest(cols = data)


# load kinase-activity inference data
# (quantile-normalized protein regressed-out phosphorylation data)
# select quantifications with more than 3 substrates
k_subN <- 3
kins <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  select(kinase, sample, kin_activity=log10P)


# function to classify the activities into neutral, activated and deactivated
act <- function(p, cutoff){
  if(p < -cutoff){activity <- "deactivated"}
  else if(p > cutoff){activity <- "activated"}
  else{activity <- "neutral"}
  activity
}


# function to select kinases to test
selKins <- function(df, ct_act = 0, ct_death = 0, logic_act = "or"){
  tab <- table(df$activity_class)
  act = tab["activated"]
  deact = tab["deactivated"]
  if(is.na(act)) act = 0
  if(is.na(deact)) deact = 0
  
  tab <- table(df$OS)
  alive = tab["0"]
  dead = tab["1"]
  if(is.na(alive)) alive = 0
  if(is.na(dead)) dead = 0
  
  if(logic_act == "or") res <- (act > ct_act | deact > ct_act) & (dead > ct_death)
  if(logic_act == "and") res <- (act > ct_act & deact > ct_act) & (dead > ct_death)
  
  res
}


# function to compare Kaplan-Meier survival curves
compare_KMCurves <- function(df){
  dat <- df
  
  km_test <- survdiff(Surv(OS_time, OS) ~ activity_class, data=dat)
  chq <- km_test$chisq
  df <- length(km_test$n)-1
  pval <- pchisq(chq, df, lower.tail = F)
  
  pval
}


# function to plot Kaplan-Meier curves
plot_KMCurves <- function(kin, cancer, df){
  dat <- df
  
  km_fit <- survfit(Surv(OS_time, OS) ~ activity_class, data=dat)
  plot <- ggsurvplot(km_fit, data = dat, pval = T, risk.table = F, censor = F) +
    labs(title = str_c(kin, cancer, sep = " - "))
  plot
}


# survival analysis using Kaplan-Meier survival estimates - by kinase and cancer type\
# classify the activities into neutral, activated and deactivated (cutoff = 1.75)\
# adjust P-values
km_surv <- kins %>%
  rename(activity_score = kin_activity) %>%
  mutate(activity_class = map_chr(activity_score, act, cutoff = 1.75)) %>%
  inner_join(samples, by = "sample") %>%
  group_by(cancer, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, selKins, ct_death = 5, ct_act = 5, logic_act = "or")) %>%
  mutate(pval = map_dbl(data, compare_KMCurves)) %>%
  #group_by(cancer) %>%
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  #ungroup() %>%
  arrange(padj)

write_tsv(km_surv %>% select(-data), "./output/files/kinases_survival_curves_logrank_test.txt")

# associations with FDR < 0.2 (adjusted log-rank test P-value)
km_surv <- km_surv %>%
  filter(padj < 0.2)


# plot heatmap
mat <- km_surv %>%
  select(kinase, cancer, padj) %>%
  mutate(padj = -log10(padj)) %>%
  pivot_wider(names_from = "kinase", values_from = "padj") %>%
  as.data.frame() %>%
  column_to_rownames(var = "cancer") %>%
  as.matrix()
rownames(mat) <- c("Lung", "Kidney", "Brain", "Breast")

heatmap <- Heatmap(
  col = c("#fdbb84", "#d7301f"),
  na_col = "#fff7ec",
  column_title = "Kinases",
  column_title_gp = gpar(fontsize = 16),
  show_heatmap_legend = T,
  border = T,
  rect_gp = gpar(col = "black", lwd = 1.5),
  column_names_centered = T,
  matrix = mat,
  name = "-log10 adjusted P-value",
  heatmap_legend_param = list(legend_width = unit(6, "cm"), direction = "horizontal", title_position = "topcenter", labels_gp = gpar(fontsize = 15), title_gp = gpar(fontsize = 17)),
  cluster_columns = F,
  cluster_rows = F,
  row_names_gp = gpar(fontsize = 16),
  column_names_gp = gpar(fontsize = 16),
  column_names_rot = 45)

draw(heatmap, heatmap_legend_side = "bottom")

pdf(file = "./output/plots/survival_analysis/kinases_KM_analysis_by_kinase_tissue_logRanK_test.pdf", height = 3, width = 5)
draw(heatmap, heatmap_legend_side = "bottom")
dev.off()

png(file="./output/plots/survival_analysis/kinases_KM_analysis_by_kinase_tissue_logRanK_test.png", units = "in", res = 300, height = 3, width = 5)
draw(heatmap, heatmap_legend_side = "bottom")
dev.off()


# plot Kaplan-Meier curves for the previous associations
#pwalk(.l = km_surv, .f = ~ print(plot_KMCurves(..1, ..2, ..3)))
#km_surv <- km_surv %>%
#  mutate(KM_plots = pmap(.l = ., .f = ~ plot_KMCurves(..1, ..2, ..3)))


# plot Kaplan-Meier curves for PRKCA associations

dat <- km_surv %>%
  filter(kinase == "PRKCA", cancer == "lung") %>%
  unnest(cols = c(data))

km_fit <- survfit(Surv(OS_time, OS) ~ activity_class, data=dat)
km_plot <- ggsurvplot(km_fit, 
                      data = dat,
                      pval = F,
                      risk.table = F,
                      censor = F, 
                      conf.int = F,
                      palette = c("red", "blue", "black"),
                      legend.labs = c("Activated", "Deactivated", "Neutral"),
                      tables.y.text.col = F,
                      risk.table.y.text.col = TRUE,
                      risk.table.title = "",
                      title = "PRKCA - Lung",
                      legend = c(0.2,0.2),
                      legend.title = "Activity",
                      xlab = "Time (days)")
km_plot

pdf(file = "./output/plots/survival_analysis/kinases_KM_analysis_PRKCA_lung.pdf", height = 4, width = 4)
print(km_plot, newpage = F)
dev.off()

png(file = "./output/plots/survival_analysis/kinases_KM_analysis_PRKCA_lung.png", height = 4, width = 4, units = "in", res = 300)
km_plot
dev.off()


dat <- km_surv %>%
  filter(kinase == "PRKCA", cancer == "kidney") %>%
  unnest(cols = c(data))

km_fit <- survfit(Surv(OS_time, OS) ~ activity_class, data=dat)
km_plot <- ggsurvplot(km_fit, 
                      data = dat,
                      pval = F,
                      risk.table = F, 
                      censor = F, 
                      conf.int = F,
                      palette = c("red", "blue", "black"),
                      legend.labs = c("Activated", "Deactivated", "Neutral"),
                      tables.y.text.col = F,
                      risk.table.y.text.col = TRUE,
                      risk.table.title = "",
                      title = "PRKCA - Kidney",
                      legend = c(0.2,0.2),
                      legend.title = "Activity",
                      xlab = "Time (days)")
km_plot

pdf(file = "./output/plots/survival_analysis/kinases_KM_analysis_PRKCA_kidney.pdf", height = 4, width = 4)
print(km_plot, newpage = F)
dev.off()

png(file = "./output/plots/survival_analysis/kinases_KM_analysis_PRKCA_kidney.png", height = 4, width = 4, units = "in", res = 300)
km_plot
dev.off()


# plot Kaplan-Meier curves for MAP3K8 and AURKB associations

dat <- km_surv %>%
  filter(kinase == "MAP3K8", cancer == "kidney") %>%
  unnest(cols = c(data))

km_fit <- survfit(Surv(OS_time, OS) ~ activity_class, data=dat)
km_plot <- ggsurvplot(km_fit, 
                      data = dat,
                      pval = F,
                      risk.table = F,
                      censor = F, 
                      conf.int = F,
                      palette = c("red", "black"),
                      legend.labs = c("Activated", "Neutral"),
                      tables.y.text.col = F,
                      risk.table.y.text.col = TRUE,
                      risk.table.title = "",
                      title = "MAP3K8 - Kidney",
                      legend = c(0.2,0.2),
                      legend.title = "Activity",
                      xlab = "Time (days)")
km_plot

pdf(file = "./output/plots/survival_analysis/kinases_KM_analysis_MAP3K8_kidney.pdf", height = 4, width = 4)
print(km_plot, newpage = F)
dev.off()

png(file = "./output/plots/survival_analysis/kinases_KM_analysis_MAP3K8_kidney.png", height = 4, width = 4, units = "in", res = 300)
km_plot
dev.off()


dat <- km_surv %>%
  filter(kinase == "AURKB", cancer == "breast") %>%
  unnest(cols = c(data))

km_fit <- survfit(Surv(OS_time, OS) ~ activity_class, data=dat)
km_plot <- ggsurvplot(km_fit, 
                      data = dat,
                      pval = F,
                      risk.table = T, 
                      censor = F, 
                      conf.int = F,
                      palette = c("red", "blue", "black"),
                      legend.labs = c("Activated", "Deactivated", "Neutral"),
                      tables.y.text.col = F,
                      risk.table.y.text.col = TRUE,
                      risk.table.title = "",
                      title = "AURKB - Breast",
                      legend = c(0.2,0.2),
                      legend.title = "Activity",
                      xlab = "Time (days)")
km_plot

pdf(file = "./output/plots/survival_analysis/kinases_KM_analysis_AURKB_breast.pdf", height = 6, width = 6)
print(km_plot, newpage = F)
dev.off()

png(file = "./output/plots/survival_analysis/kinases_KM_analysis_AURKB_breast.png", height = 6, width = 6, units = "in", res = 300)
km_plot
dev.off()


# survival analysis using cox proportional hazards regression models - by kinase and cancer type adjusting for additional covariates

# function to fit the cox proportional hazards regression
cox_reg <- function(df, covariates, adjust = F){
  
  # remove factors with less then 2 levels
  # remove numeric variables with null variance
  covariates <- covariates %>%
    filter(sample %in% df$sample) %>%
    mutate_if(.predicate = is.character, .funs = as.factor) %>%
    mutate_if(.predicate = is.factor, .funs = fct_drop) %>%
    select_if(.predicate = ~ if(!is.factor(.x)){TRUE}else{if(nlevels(.x) == 1){FALSE}else{TRUE}}) %>%
    select_if(.predicate = ~ if(!is.numeric(.x)){TRUE}else{if(var(.x, na.rm=T) > 0){TRUE}else{FALSE}})
  
  dat <- df %>%
    inner_join(covariates, by = "sample")
  
  # compute cox hazard regression
  if(adjust){
    f <- as.formula(str_c("Surv(OS_time, OS) ~ ",  str_c(c("activity_score", colnames(dat)[6:length(colnames(dat))]), collapse = "+")))
    cox_test <- tryCatch(coxph(formula = f, data=dat), error=function(e){NA})
    if(is.na(cox_test)){
      f <- as.formula(str_c("Surv(OS_time, OS) ~ ",  str_c(c("activity_score", c("gender", "age")), collapse = "+")))
      cox_test <- tryCatch(coxph(formula = f, data=dat), error=function(e){NA})
    }
  } else {
    cox_test <- tryCatch(coxph(Surv(OS_time, OS) ~ activity_score, data=dat), error=function(e){NA})
  }
  
  # statistics summary
  if(is.na(cox_test)){
    act_pval <- NA
    act_beta <- NA
    act_hr <- NA
    hrCI95_lower <- NA
    hrCI95_upper <- NA
    lrt_pval <- NA
  } else {
    sm <- summary(cox_test)
    act_pval <- sm$coefficients["activity_score", "Pr(>|z|)"]
    act_beta <- sm$coefficients["activity_score", "coef"]
    act_hr <- sm$coefficients["activity_score", "exp(coef)"]
    hrCI95_lower <- sm$conf.int["activity_score", "lower .95"]
    hrCI95_upper <- sm$conf.int["activity_score", "upper .95"]
    lrt_pval <- sm$logtest["pvalue"]
  }
  
  # test the proportional hazards assumption
  # independence of Schoenfeld residuals over time
  test_ph <- tryCatch(cox.zph(cox_test), error=function(e){NA})
  
  if(is.na(test_ph)){
    ph_act_pval <- NA
    ph_global_pval <- NA
  } else {
    ph_act_pval <- test_ph$table["activity_score", "p"]
    ph_global_pval <- test_ph$table["GLOBAL", "p"]
  }
  
  res <- tibble(act_pval=act_pval, act_beta=act_beta, act_hr=act_hr, hrCI95_lower=hrCI95_lower, hrCI95_upper=hrCI95_upper, lrt_pval=lrt_pval, ph_act_pval=ph_act_pval, ph_global_pval=ph_global_pval)
  
  return(res)
}

covars <- clinical_data %>%
  filter(sample %in% samples$sample) %>%
  select(sample, age, gender)

samp_mut = 100
mutated_genes <- read_tsv(file = "./output/files/mutations_matrix.txt.gz") %>%
  filter(n >= samp_mut)

genotypes <- mutated_genes %>%
  select(-n) %>%
  as.data.frame() %>%
  column_to_rownames(var = "gene") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  as_tibble()

covars <- covars %>%
  inner_join(genotypes, by = "sample")

cox_surv <- kins %>%
  rename(activity_score = kin_activity) %>%
  mutate(activity_class = map_chr(activity_score, act, cutoff = 1.75)) %>%
  inner_join(samples, by = "sample") %>%
  group_by(cancer, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, selKins, ct_death = 5, ct_act = 5, logic_act = "or")) %>%
  mutate(res = map(.x=data, .f=cox_reg, covariates = covars, adjust=T)) %>%
  unnest(res) %>%
  filter(!is.na(act_pval)) %>%
  #filter(!is.na(ph_act_pval)) %>%
  #group_by(cancer) %>%
  mutate(padj1 = p.adjust(act_pval, method = "BH"), padj2 = p.adjust(lrt_pval, method = "BH")) %>%
  #ungroup() %>%
  arrange(padj1)

write_tsv(cox_surv %>% select(-data), "./output/files/kinases_cox_regression.txt")

cox_surv <- cox_surv %>%
  #filter(padj1 < 0.05, padj2 < 0.05) %>%
  filter(padj1 < 0.05)
  #group_by(cancer) %>%
  #slice_min(order_by = padj1, n=10) %>%
  #ungroup()
cox_surv

dat <- cox_surv %>%
  filter(kinase == "KSR1", cancer == "liver") %>%
  unnest(cols = c(data))

km_fit <- survfit(Surv(OS_time, OS) ~ activity_class, data=dat)
km_plot <- ggsurvplot(km_fit, 
                      data = dat,
                      pval = F,
                      risk.table = F, 
                      censor = F, 
                      conf.int = F,
                      palette = c("red", "blue", "black"),
                      legend.labs = c("Activated", "Deactivated", "Neutral"),
                      tables.y.text.col = F,
                      risk.table.y.text.col = TRUE,
                      risk.table.title = "",
                      title = "KSR1 - Liver",
                      legend = c(0.2,0.2),
                      legend.title = "Activity",
                      xlab = "Time (days)")
km_plot


dat <- kins %>%
  rename(activity_score = kin_activity) %>%
  mutate(activity_class = map_chr(activity_score, act, cutoff = 1.75)) %>%
  inner_join(samples, by = "sample") %>%
  group_by(cancer, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, selKins, ct_death = 5, ct_act = 5, logic_act = "or")) %>%
  filter(kinase == "CAMK2A", cancer == "lung") %>%
  select(-kinase,-cancer) %>%
  unnest(cols = c(data))

covariates <- covars %>%
  filter(sample %in% dat$sample) %>%
  mutate_if(.predicate = is.character, .funs = as.factor) %>%
  mutate_if(.predicate = is.factor, .funs = fct_drop) %>%
  select_if(.predicate = ~ if(!is.factor(.x)){TRUE}else{if(nlevels(.x) == 1){FALSE}else{TRUE}}) %>%
  select_if(.predicate = ~ if(!is.numeric(.x)){TRUE}else{if(var(.x, na.rm=T) > 0){TRUE}else{FALSE}})

dat <- dat %>%
  inner_join(covariates, by = "sample")

f <- as.formula(str_c("Surv(OS_time, OS) ~ ",  str_c(c("activity_score", colnames(dat)[6:length(colnames(dat))]), collapse = "+")))
cox_test <- coxph(formula = f, data=dat)
