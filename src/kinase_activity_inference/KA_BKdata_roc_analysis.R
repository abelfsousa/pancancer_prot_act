#' ---
#' title: "ROC curves from kinase-activity inference using benchmark phosphorylation data"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)


#' load R packages
library(RColorBrewer)
library(tidyverse)
library(viridis)
library(PRROC)
library(ROCR)


#' source R kinase-activity inference pipeline
source("./src/utils/infer_kinase_activity.R")



#' load kinase activation scores
# non-weighted z-test
nwz_KAinferences <- read_tsv(file = "./output/files/BKphospho_kinase_activity_ztestNW.txt.gz") %>%
  mutate(test = "z-test")

# weighted z-test
wz_KAinferences <- read_tsv(file = "./output/files/BKphospho_kinase_activity_ztestW.txt.gz") %>%
  mutate(test = "z-test-weighted")


#' load known condition-specific kinase regulation
kreg <- read_tsv(file = "./data/kinase_substrate/benchmark_dataset/kinase_condition_pairs.txt") %>%
  select(sample=Condition,kinase=Kinase) %>%
  mutate(regulated = 1)



#' # roc plot by (non-)weighted z-test and comparing KS list - Randomization
roc <- nwz_KAinferences %>%
  bind_rows(wz_KAinferences) %>%
  filter(source_type %in% c("in_vitro", "in_vivo", "database", "text-mining")) %>%
  mutate(source_type = str_replace(source_type, "in_vitro", "in vitro")) %>%
  mutate(source_type = str_replace(source_type, "in_vivo", "in vivo")) %>%
  nest(data = c(sample,kinase,n,log10P)) %>%
  mutate(roc = map(.x=data, .f=compute_prroc, gold_std=kreg, method="roc2", randomizeNeg = T, randomizations = 100, remove_KS_pairs=F, SubN=NULL)) %>%
  select(-data) %>%
  unnest(cols = roc) %>%
  group_by(test, source_type, sampling) %>%
  nest() %>%
  ungroup() %>%
  mutate(group = 1:n()) %>%
  unnest(cols = data) %>%
  mutate(group = as.factor(group)) %>%
  mutate(source_type = fct_rev(fct_reorder(source_type, auc, mean))) %>%
  mutate(test = fct_rev(fct_reorder(test, auc, mean)))

roc_auc <- roc %>%
  group_by(source_type, test, sampling) %>%
  summarise(auc = unique(auc)) %>%
  summarise(auc = mean(auc)) %>%
  mutate(auc = format(round(auc, 2), nsmall = 2)) %>%
  ungroup() %>%
  arrange(test, source_type) %>%
  #mutate(x = 0.1, y = rep(c(1, 0.95, 0.9, 0.85), 2))
  mutate(x = 0.7,
         y = (function(x, y, z) {
           a <- x;
           for(i in 1:(z-1)){
             b <- a[i] + y
             a <- c(a, b)}
           rep(rev(a),2)}) (0.18, 0.08, 4))

roc_curve <- roc %>%
  ggplot() +
  geom_line(mapping=aes(x = fpr, y = tpr, color = source_type, group = group, alpha = source_type)) +
  geom_label(data=roc_auc, mapping=aes(x = x, y = y, label = auc, color = source_type), size = 6, label.padding = unit(0.15, "lines"), show.legend = FALSE, fontface = "bold") +
  annotate("text", x = 0.81, y = 0.5, label = "Mean AUC:", size = 6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  theme_classic() +
  facet_wrap(~ test, labeller = labeller(test = c("z-test" = "Z-test", "z-test-weighted" = "Weighted Z-test"))) +
  theme(
    plot.margin = unit(c(0, 0.5, 0, 0), "cm"),
    axis.title=element_text(colour="black", size=20),
    axis.text=element_text(colour="black", size=18),
    strip.text=element_text(colour="black", size=20),
    strip.background=element_blank(),
    legend.text=element_text(size=18),
    legend.title=element_text(size=20),
    panel.spacing = unit(2, "lines"),
    legend.position = "bottom") +
  scale_alpha_manual(values = c(0.2,0.8,0.2,0.8), guide = F) +
  #scale_color_brewer(type = "div", palette = "Spectral", labels = c("Text-mining", "Databases", "In vivo", "In vitro")) +
  #scale_color_viridis(discrete = T, labels = c("Text-mining", "Databases", "In vivo", "In vitro")) +
  scale_color_manual(values = c("#67a9cf", "#d8b365", "#de2d26", "#78c679"), labels = c("Text-mining", "Databases", "In vivo", "In vitro")) +
  labs(x = "False positive rate", y = "True positive rate", color = "Source of kinase targets") +
  guides(color = guide_legend(override.aes = list(alpha=1, size = 2), nrow = 2))

#+ fig.width=8, fig.height=4
roc_curve

ggsave("bk_phospho_kinase_lists_roc_randomized.png", plot=roc_curve, path="./output/plots/kinase_activity_inference/", height=6, width=8)
ggsave("bk_phospho_kinase_lists_roc_randomized.pdf", plot=roc_curve, path="./output/plots/kinase_activity_inference/", height=6, width=8)



#' # roc plot by (non-)weighted z-test and comparing KS list - Randomization\
#' # another version with some controls
getCommonPairs <- function(x){
  l <- x
  y <- map(l, ~ select(.x, sample, kinase))
  common_pairs <- reduce(y, intersect)
  l <- map(l, ~ semi_join(.x, common_pairs, by = c("sample", "kinase")) %>% arrange(sample, kinase))
  
  l
}

allEqual <- function(x){
  l <- c()
  for(i in 1:(length(x)-1)){
    for(j in (i+1):length(x)){
      l <- c(l, identical(x[[i]], x[[j]]))
    }
  }
  res <- all(l)
  #res <- l
  res
}

sample_pairs <- function(KA_inference, gold_std, randomizations = 100, seed=123){
  
  KA_inference <- KA_inference %>%
    left_join(gold_std, by = c("sample", "kinase")) %>%
    mutate(regulated = replace_na(regulated, 0))
  
  
  positives <- KA_inference %>%
    filter(regulated == 1)
  negatives <- KA_inference %>%
    filter(regulated == 0)
  
  set.seed(seed)
  
  neg_sets <- tibble(sampling = 1:randomizations) %>%
    mutate(data = map(
      .x = sampling,
      .f = ~ negatives %>% sample_n(size = nrow(positives)) %>% bind_rows(positives))) %>%
  unnest(cols = data)
  
  return(neg_sets)
}

calculate_roc <- function(KA_inference, gold_std, randomizations = 100, seed=123){
  
  KA_inference <- KA_inference %>%
    left_join(gold_std, by = c("sample", "kinase")) %>%
    mutate(regulated = replace_na(regulated, 0))
  
  
  positives <- KA_inference %>%
    filter(regulated == 1)
  negatives <- KA_inference %>%
    filter(regulated == 0)
  
  set.seed(seed)
  
  neg_sets <- tibble(sampling = 1:randomizations) %>%
    mutate(data = map(
      .x = sampling,
      .f = ~ negatives %>% sample_n(size = nrow(positives))))
    #unnest(cols = data)
  
  res <- neg_sets %>%
    mutate(curve = map(
      .x = data,
      .f = function(x){
        roccurve <- roc.curve(positives[, "log10P", drop=T], x[, "log10P", drop=T], curve = T)
        auc <- roccurve$auc
        res <- tibble(fpr=roccurve$curve[,1], tpr=roccurve$curve[,2], cutoff=roccurve$curve[,3], auc=auc)})) %>%
    select(-data) %>%
    unnest(cols = curve)
  
  return(res)
}


roc <- nwz_KAinferences %>%
  bind_rows(wz_KAinferences) %>%
  filter(source_type %in% c("in_vitro", "in_vivo", "database", "text-mining")) %>%
  mutate(source_type = str_replace(source_type, "in_vitro", "in vitro")) %>%
  mutate(source_type = str_replace(source_type, "in_vivo", "in vivo")) %>%
  filter(n > 0) %>%
  mutate(log10P = abs(log10P)) %>%
  nest(data = c(sample,kinase,n,log10P)) %>%
  mutate(data = getCommonPairs(data))

allEqual(map(roc$data, ~ select(.x,1:2)))

all_pairs <- roc %>%
  mutate(pairs = map(.x=data, .f=sample_pairs, gold_std=kreg, randomizations = 100))

allEqual(map(all_pairs$pairs, ~ select(.x,-4:-5)))

roc <- roc %>%
  mutate(roc = map(.x=data, .f=calculate_roc, gold_std=kreg, randomizations = 100)) %>%
  select(-data) %>%
  unnest(cols = roc) %>%
  group_by(test, source_type, sampling) %>%
  nest() %>%
  ungroup() %>%
  mutate(group = 1:n()) %>%
  unnest(cols = data) %>%
  mutate(group = as.factor(group)) %>%
  mutate(source_type = fct_rev(fct_reorder(source_type, auc, mean))) %>%
  mutate(test = fct_rev(fct_reorder(test, auc, mean)))


roc_auc <- roc %>%
  group_by(source_type, test, sampling) %>%
  summarise(auc = unique(auc)) %>%
  summarise(auc = mean(auc)) %>%
  mutate(auc = format(round(auc, 2), nsmall = 2)) %>%
  ungroup() %>%
  arrange(test, source_type) %>%
  #mutate(x = 0.1, y = rep(c(1, 0.95, 0.9, 0.85), 2))
  mutate(x = 0.7,
         y = (function(x, y, z) {
           a <- x;
           for(i in 1:(z-1)){
             b <- a[i] + y
             a <- c(a, b)}
           rep(rev(a),2)}) (0.18, 0.08, 4))

roc_curve <- roc %>%
  ggplot() +
  geom_line(mapping=aes(x = fpr, y = tpr, color = source_type, group = group, alpha = source_type)) +
  geom_label(data=roc_auc, mapping=aes(x = x, y = y, label = auc, color = source_type), size = 6, label.padding = unit(0.15, "lines"), show.legend = FALSE, fontface = "bold") +
  annotate("text", x = 0.81, y = 0.5, label = "Mean AUC:", size = 6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  theme_classic() +
  facet_wrap(~ test, labeller = labeller(test = c("z-test" = "Z-test", "z-test-weighted" = "Weighted Z-test"))) +
  theme(
    plot.margin = unit(c(0, 0.5, 0, 0), "cm"),
    axis.title=element_text(colour="black", size=20),
    axis.text=element_text(colour="black", size=18),
    strip.text=element_text(colour="black", size=20),
    strip.background=element_blank(),
    legend.text=element_text(size=18),
    legend.title=element_text(size=20),
    panel.spacing = unit(2, "lines"),
    legend.position = "bottom") +
  scale_alpha_manual(values = c(0.2,0.8,0.2,0.8), guide = F) +
  #scale_color_brewer(type = "div", palette = "Spectral", labels = c("Text-mining", "Databases", "In vivo", "In vitro")) +
  #scale_color_viridis(discrete = T, labels = c("Text-mining", "Databases", "In vivo", "In vitro")) +
  scale_color_manual(values = c("#67a9cf", "#d8b365", "#de2d26", "#78c679"), labels = c("Text-mining", "Databases", "In vivo", "In vitro")) +
  labs(x = "False positive rate", y = "True positive rate", color = "Source of kinase targets") +
  guides(color = guide_legend(override.aes = list(alpha=1, size = 2), nrow = 2))

#+ fig.width=8, fig.height=4
roc_curve

ggsave("bk_phospho_kinase_lists_roc_randomized2.png", plot=roc_curve, path="./output/plots/kinase_activity_inference/", height=6, width=8)
ggsave("bk_phospho_kinase_lists_roc_randomized2.pdf", plot=roc_curve, path="./output/plots/kinase_activity_inference/", height=6, width=8)




#' # roc plot for non-weighted z-test and comparing KS list - randomization and average ROC curve
roc <- nwz_KAinferences %>%
  filter(source_type %in% c("in_vitro", "in_vivo", "database", "text-mining")) %>%
  mutate(source_type = str_replace(source_type, "in_vitro", "in vitro")) %>%
  mutate(source_type = str_replace(source_type, "in_vivo", "in vivo")) %>%
  nest(data = c(sample,kinase,n,log10P)) %>%
  mutate(roc = map(.x=data, .f=compute_rocr, gold_std=kreg, randomizations = 100, remove_KS_pairs=F, SubN=NULL))

pdf(file = "./output/plots/kinase_activity_inference/bk_phospho_kinase_lists_roc_randomized_averaged.pdf", height=6, width=6)
plot(roc$roc[[1]][[1]], avg='vertical', spread.estimate='stddev', lwd=3, main='Kinase Targets Performance', col='orange', xaxis.cex.axis = 1.3, yaxis.cex.axis = 1.3, cex.lab = 1.5, cex.main = 1.5)
plot(roc$roc[[2]][[1]], avg='vertical', spread.estimate='stddev', lwd=3, main='', col='red', add=T)
plot(roc$roc[[3]][[1]], avg='vertical', spread.estimate='stddev', lwd=3, main='', col='#0099FF', add=T)
plot(roc$roc[[4]][[1]], avg='vertical', spread.estimate='stddev', lwd=3, main='', col='#009933', add=T)
abline(coef = c(0,1), lty = "dashed")
legend(0.65, 0.4, title = "Mean AUC", legend=c(round(roc$roc[[4]][[2]],2), round(roc$roc[[3]][[2]],2), round(roc$roc[[2]][[2]],2), round(roc$roc[[1]][[2]],2)), col=c("#009933", "#0099FF", "red", "orange"), lty = 1, lwd = 3, cex=1.3)
dev.off()


#' # another version with some controls (using the same kinase-condition pairs across lists)
compute_rocr_control <- function(kinase_list){
  
  sets <- kinase_list %>%
    select(sampling, log10P, regulated) %>%
    chop(c(log10P, regulated))
  
  predictions <- as.list(sets$log10P)
  labels <- as.list(sets$regulated)
  
  pred <- prediction(predictions, labels)
  perf <- performance(pred, measure = "tpr", x.measure = "fpr")
  perf2 <- performance(pred, measure = "auc")
  
  auc_avg <- sum(unlist(perf2@y.values))/length(perf2@y.values)
  
  return(list(perf, auc_avg))
}

kinase_sets <- nwz_KAinferences %>%
  filter(source_type %in% c("in_vitro", "in_vivo", "database", "text-mining")) %>%
  mutate(source_type = str_replace(source_type, "in_vitro", "in vitro")) %>%
  mutate(source_type = str_replace(source_type, "in_vivo", "in vivo")) %>%
  filter(n > 0) %>%
  mutate(log10P = abs(log10P)) %>%
  nest(data = c(sample,kinase,n,log10P)) %>%
  mutate(data = getCommonPairs(data))

allEqual(map(kinase_sets$data, ~ select(.x,1:2)))

all_pairs <- kinase_sets %>%
  mutate(pairs = map(.x=data, .f=sample_pairs, gold_std=kreg, randomizations = 100))

allEqual(map(all_pairs$pairs, ~ select(.x,-4:-5)))

all_pairs <- all_pairs %>%
  mutate(roc = map(.x = pairs, .f = compute_rocr_control))

pdf(file = "./output/plots/kinase_activity_inference/bk_phospho_kinase_lists_roc_randomized_averaged2.pdf", height=6, width=6)
plot(all_pairs$roc[[1]][[1]], avg='vertical', spread.estimate='stddev', lwd=3, main='Kinase Targets Performance using\nsame kinase-condition pairs', col='orange', xaxis.cex.axis = 1.3, yaxis.cex.axis = 1.3, cex.lab = 1.5, cex.main = 1.5)
plot(all_pairs$roc[[2]][[1]], avg='vertical', spread.estimate='stddev', lwd=3, main='', col='red', add=T)
plot(all_pairs$roc[[3]][[1]], avg='vertical', spread.estimate='stddev', lwd=3, main='', col='#0099FF', add=T)
plot(all_pairs$roc[[4]][[1]], avg='vertical', spread.estimate='stddev', lwd=3, main='', col='#009933', add=T)
abline(coef = c(0,1), lty = "dashed")
legend(0.65, 0.4, title = "Mean AUC", legend=c(round(all_pairs$roc[[4]][[2]],2), round(all_pairs$roc[[3]][[2]],2), round(all_pairs$roc[[2]][[2]],2), round(all_pairs$roc[[1]][[2]],2)), col=c("#009933", "#0099FF", "red", "orange"), lty = 1, lwd = 3, cex=1.3)
dev.off()





#' # roc plot by (non-)weighted z-test and comparing KS list
roc <- nwz_KAinferences %>%
  bind_rows(wz_KAinferences) %>%
  filter(source_type %in% c("in_vitro", "in_vivo", "database", "text-mining")) %>%
  nest(data = c(sample,kinase,n,log10P)) %>%
  mutate(roc = map(.x=data, .f=compute_prroc, gold_std=kreg, method="roc1")) %>%
  select(-data) %>%
  unnest(cols = roc) %>%
  mutate(auc = round(auc, 2)) %>%
  rename(list=source_type) %>%
  mutate(list = fct_rev(fct_reorder(list, auc, mean))) %>%
  mutate(test = fct_rev(fct_reorder(test, auc, mean)))

roc_auc <- roc %>%
  group_by(list, test) %>%
  summarise(auc = unique(auc)) %>%
  ungroup() %>%
  arrange(test, list) %>%
  mutate(x = 0.1, y = rep(c(1, 0.90, 0.8, 0.7), 2))

roc_curve <- roc %>%
  ggplot(mapping=aes(x = fpr, y = tpr, color = list)) +
  geom_line() +
  geom_text(data=roc_auc, mapping=aes(x = x, y = y, label = auc, color = list), show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~ test, labeller = labeller(test = c("z-test" = "z-test", "z-test-weighted" = "weighted z-test"))) +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    strip.text=element_text(colour="black", size=13),
    strip.background=element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=13),
    panel.spacing = unit(1, "lines")) +
  scale_color_viridis(discrete = T) +
  labs(x = "FPR", y = "TPR", title = "")

#+ fig.width=7, fig.height=4
roc_curve

ggsave("bk_phospho_kinase_lists_roc1.png", plot=roc_curve, path="./output/plots/kinase_activity_inference/", height=4, width=7)
unlink("bk_phospho_kinase_lists_roc1.png")



#' # roc plot by KS list and comparing (non-)weighted z-test
roc <- nwz_KAinferences %>%
  bind_rows(wz_KAinferences) %>%
  filter(source_type %in% c("in_vitro", "in_vivo", "database", "text-mining")) %>%
  nest(data = c(sample,kinase,n,log10P)) %>%
  mutate(roc = map(.x=data, .f=compute_prroc, gold_std=kreg, method="roc1")) %>%
  select(-data) %>%
  unnest(cols = roc) %>%
  mutate(auc = round(auc, 2)) %>%
  mutate(test = fct_rev(fct_reorder(test, auc, mean))) %>%
  mutate(source_type = fct_rev(fct_reorder(source_type, auc, mean)))

roc_auc <- roc %>%
  group_by(source_type, test) %>%
  summarise(auc = unique(auc)) %>%
  ungroup() %>%
  mutate(x = 0.1, y = rep(c(1, 0.9), 4))

roc_curve <- roc %>%
  ggplot(mapping=aes(x = fpr, y = tpr, color = test)) +
  geom_line() +
  geom_text(data=roc_auc, mapping=aes(x = x, y = y, label = auc, color = test), show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~ source_type, labeller = labeller(source_type = c("text-mining" = "text-mining", "database" = "database", "in_vivo" = "in vivo", "in_vitro" = "in vitro"))) +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    strip.text=element_text(colour="black", size=13),
    strip.background=element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=13),
    panel.spacing = unit(1, "lines")) +
  scale_color_viridis(discrete = T) +
  labs(x = "FPR", y = "TPR", title = "")

#+ fig.width=7, fig.height=6
roc_curve

ggsave("bk_phospho_kinase_lists_roc2.png", plot=roc_curve, path="./output/plots/kinase_activity_inference/", height=6, width=7)
unlink("bk_phospho_kinase_lists_roc2.png")


#' # prc plot by (non-)weighted z-test and comparing KS list
prc <- nwz_KAinferences %>%
  bind_rows(wz_KAinferences) %>%
  filter(source_type %in% c("in_vitro", "in_vivo", "database", "text-mining")) %>%
  nest(data = c(sample,kinase,n,log10P)) %>%
  mutate(prc = map(.x=data, .f=compute_prroc, gold_std=kreg, method = "prc")) %>%
  select(-data) %>%
  unnest(cols = prc) %>%
  mutate(auc = round(auc, 2)) %>%
  rename(list=source_type) %>%
  mutate(list = fct_rev(fct_reorder(list, auc, mean))) %>%
  mutate(test = fct_rev(fct_reorder(test, auc, mean)))

prc_auc <- prc %>%
  group_by(list, test) %>%
  summarise(auc = unique(auc)) %>%
  ungroup() %>%
  arrange(test, list) %>%
  mutate(x = 0.1, y = rep(c(1, 0.90, 0.8, 0.7), 2))

prc_curve <- prc %>%
  ggplot(mapping=aes(x = recall, y = precision, color = list)) +
  geom_line() +
  geom_text(data=prc_auc, mapping=aes(x = x, y = y, label = auc, color = list), show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~ test, labeller = labeller(test = c("z-test" = "z-test", "z-test-weighted" = "weighted z-test"))) +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    strip.text=element_text(colour="black", size=13),
    strip.background=element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=13),
    panel.spacing = unit(1, "lines")) +
  scale_color_viridis(discrete = T) +
  labs(x = "Recall", y = "Precision", title = "")

#+ fig.width=7, fig.height=4
prc_curve

ggsave("bk_phospho_kinase_lists_prc1.png", plot=prc_curve, path="./output/plots/kinase_activity_inference/", height=4, width=7)
unlink("bk_phospho_kinase_lists_prc1.png")


#' # prc plot by KS list and comparing (non-)weighted z-test
prc <- nwz_KAinferences %>%
  bind_rows(wz_KAinferences) %>%
  filter(source_type %in% c("in_vitro", "in_vivo", "database", "text-mining")) %>%
  nest(data = c(sample,kinase,n,log10P)) %>%
  mutate(prc = map(.x=data, .f=compute_prroc, gold_std=kreg, method = "prc")) %>%
  select(-data) %>%
  unnest(cols = prc) %>%
  mutate(auc = round(auc, 2)) %>%
  mutate(test = fct_rev(fct_reorder(test, auc, mean))) %>%
  mutate(source_type = fct_rev(fct_reorder(source_type, auc, mean)))

prc_auc <- prc %>%
  group_by(source_type, test) %>%
  summarise(auc = unique(auc)) %>%
  ungroup() %>%
  mutate(x = 0.1, y = rep(c(1, 0.9), 4))

prc_curve <- prc %>%
  ggplot(mapping=aes(x = recall, y = precision, color = test)) +
  geom_line() +
  geom_text(data=prc_auc, mapping=aes(x = x, y = y, label = auc, color = test), show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~ source_type, labeller = labeller(source_type = c("text-mining" = "text-mining", "database" = "database", "in_vivo" = "in vivo", "in_vitro" = "in vitro"))) +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12),
    strip.text=element_text(colour="black", size=13),
    strip.background=element_blank(),
    legend.text=element_text(size=12),
    legend.title=element_text(size=13),
    panel.spacing = unit(1, "lines")) +
  scale_color_viridis(discrete = T) +
  labs(x = "Recall", y = "Precision", title = "")

#+ fig.width=7, fig.height=6
prc_curve

ggsave("bk_phospho_kinase_lists_prc2.png", plot=prc_curve, path="./output/plots/kinase_activity_inference/", height=6, width=7)
unlink("bk_phospho_kinase_lists_prc2.png")
