# Survival analysis using RPPA phosphosites from TCGA

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


# load TCGA clinical data
clinical_data <- read_tsv(file = "./data/clinical/tcga_publications/survival_data/tcga_clinical_data.txt")
clinical_data <- clinical_data %>%
  select(sample = bcr_patient_barcode, cancer_type = type, OS, OS_time=OS.time) %>%
  filter(if_all(.cols = starts_with("OS"), .fns = ~ !is.na(.x))) %>%
  filter(if_all(.cols = starts_with("OS"), .fns = ~ !str_detect(.x, "[A-Z]+"))) %>%
  mutate(across(.cols = starts_with("OS"), .fns = as.numeric)) %>%
  mutate(sample = str_replace_all(sample, "-", "."))


# load RPPA phosphorylation data
load("./data/protein/tcga/rppa/pan_can_tcga_rppa_mat_updated.Rdata")
pan_can_tcga_rppa_mat <- pan_can_tcga_rppa_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  as_tibble() %>%
  #select_at(.vars = c(1, which(str_sub(colnames(.),14,15) == "01"))) %>%
  select(feature, ends_with("-01")) %>%
  pivot_longer(-feature, names_to = "sample", values_to = "expr") %>%
  mutate(sample = str_replace_all(str_sub(sample, 1, 12), "-", "."))

load("./data/protein/tcga/rppa/matched_prot_phosphoprot_list_updated.Rdata")
matched_prot_phosphoprot_list <- matched_prot_phosphoprot_list %>%
  as_tibble()

load("./data/protein/tcga/rppa/kin_ps_relation.Rdata")
kin_ps_relation <- kin_ps_relation %>%
  as_tibble() %>%
  select(kinase=Protein, psite=Ab, activity=Activity)

rppa_pho <- pan_can_tcga_rppa_mat %>%
  inner_join(matched_prot_phosphoprot_list, by = c("feature" = "fil_rppa_phospho_ab")) %>%
  select(psite = feature, protein = fil_rppa_prot_ab, gene = fil_gene_prot_list, sample, psite_expr = expr) %>%
  mutate(across(.cols = where(is.factor), .fns = as.character))


# function to classify the activities into neutral, activated and deactivated
act <- function(z, z_cutoff){
  if(z < -z_cutoff){activity <- "deactivated"}
  else if(z > z_cutoff){activity <- "activated"}
  else{activity <- "neutral"}
  activity
}


# function to select phosphosites to test
selPsites <- function(df, ct_act = 0, ct_death = 0, logic_act = "or"){
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


# survival analysis using Kaplan-Meier survival estimates - for all possible phosphosites and cancer types
# classify the phosphosites as neutral, activated and deactivated across samples (cutoff = zscore =/>/<2 )
# adjust P-values
km_surv <- rppa_pho %>%
  inner_join(clinical_data, by = "sample") %>%
  group_by(cancer_type, psite) %>%
  mutate(z = scale(psite_expr)[,1]) %>%
  ungroup() %>%
  mutate(activity_class = map_chr(.x = z, .f = act, z_cutoff = 2)) %>%
  group_by(cancer_type, protein, gene, psite) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, selPsites, ct_death = 5, ct_act = 5, logic_act = "or")) %>%
  mutate(pval = map_dbl(data, compare_KMCurves)) %>%
  mutate(padj = p.adjust(pval, method = "BH")) %>%
  arrange(padj)


# survival analysis using Kaplan-Meier survival estimates - for the phosphosites/cancer type pairs reported in the article
# classify the phosphosites as neutral, activated and deactivated across samples (cutoff = zscore =/>/<2 )
km_surv <- rppa_pho %>%
  inner_join(clinical_data, by = "sample") %>%
  filter(cancer_type == "KIRC" & gene %in% c("PRKCA", "RAF1")) %>%
  group_by(cancer_type, psite) %>%
  mutate(z = scale(psite_expr)[,1]) %>%
  ungroup() %>%
  mutate(activity_class = map_chr(.x = z, .f = act, z_cutoff = 2)) %>%
  group_by(cancer_type, protein, gene, psite) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, selPsites, ct_death = 5, ct_act = 5, logic_act = "or")) %>%
  mutate(pval = map_dbl(data, compare_KMCurves))


# plot Kaplan-Meier curves for PRKCA association in kidney cancer
dat <- km_surv %>%
  filter(gene == "PRKCA", cancer_type == "KIRC") %>%
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
                      title = "PRKCA pS657 - Kidney",
                      legend = c(0.2,0.2),
                      legend.title = "Activity",
                      xlab = "Time (days)")
km_plot

pdf(file = "./output/plots/survival_analysis/kinases_KM_analysis_PRKCA_kidney_rppa.pdf", height = 6, width = 6)
print(km_plot, newpage = F)
dev.off()

png(file = "./output/plots/survival_analysis/kinases_KM_analysis_PRKCA_kidney_rppa.png", height = 6, width = 6, units = "in", res = 300)
km_plot
dev.off()


# plot Kaplan-Meier curves for RAF1 association in kidney cancer
dat <- km_surv %>%
  filter(gene == "RAF1", cancer_type == "KIRC") %>%
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
                      title = "RAF1 - Kidney",
                      legend = c(0.2,0.2),
                      legend.title = "Activity",
                      xlab = "Time (days)")
km_plot


# rppa_pho %>%
#   group_by(psite) %>%
#   mutate(z = scale(psite_expr)[,1]) %>%
#   ungroup() %>%
#   mutate(activity_class = if_else(z < -2, "deactivated", if_else(z > 2, "activated", "neutral"))) 
# 
# rppa_pho %>%
#   group_by(psite) %>%
#   mutate(z = scale(psite_expr)[,1]) %>%
#   #mutate(z = (function(x) {scale(x)[,1]}) (psite_expr)) %>%
#   ungroup() %>%
#   mutate(activity_class = if_else(z < -2, "deactivated", if_else(z > 2, "activated", "neutral"))) %>%
#   #mutate(activity_class = if_else(z > 2, "activated", if_else(z < -2, "deactivated", "neutral"))) %>%
#   group_by(psite) %>%
#   mutate(MEAN = mean(psite_expr), SD = sd(psite_expr)) %>%
#   ungroup() %>%
#   mutate(activity_class2 = if_else(psite_expr < (MEAN-(2*SD)), "deactivated", if_else(psite_expr > (MEAN+(2*SD)), "activated", "neutral"))) %>%
#   summarise(equal = (sum(activity_class == activity_class2) == nrow(.)))
# 
# # function to classify the activities into neutral, activated and deactivated
# act <- function(df, z_cutoff){
#   
#   df <- df %>%
#     mutate(z = scale(psite_expr)[,1]) %>%
#     mutate(act = if_else(z < -z_cutoff, "deactivated", if_else(z > z_cutoff, "activated", "neutral"))) %>%
#     select(z, act)
#   
#   df
# }
# 
# rppa_pho %>%
#   group_by(psite) %>%
#   nest() %>%
#   ungroup() %>%
#   mutate(activity_class = map(.x = data, .f = act, z_cutoff = 2)) %>%
#   unnest(cols = c(data, act))
# 
# 
# # function to classify the activities into neutral, activated and deactivated
# act <- function(z, z_cutoff){
#   if(z < -z_cutoff){activity <- "deactivated"}
#   else if(z > z_cutoff){activity <- "activated"}
#   else{activity <- "neutral"}
#   activity
# }
# 
# rppa_pho %>%
#   group_by(psite) %>%
#   mutate(z = scale(psite_expr)[,1]) %>%
#   ungroup() %>%
#   mutate(activity_class = map_chr(.x = z, .f = act, z_cutoff = 2))
# 

