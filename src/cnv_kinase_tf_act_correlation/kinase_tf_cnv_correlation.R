#' ---
#' title: "CNV to kinase/TF activity correlation"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=6, fig.height=6)


#' load R packages
library(tidyverse)
library(ggpubr)
library(viridis)

source(file = "./src/utils/getSamples.R")

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy



#' load CNV data
cnv <- read_tsv(file = "./output/files/cnv.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "cnv")


#' load kinase-activity inference data\
#' (quantile-normalized protein regressed-out phosphorylation data)
ka <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz")

k_subN <- 3
k_sampN <- 0
kin_activity <- ka %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  group_by(kinase) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > k_sampN) %>%
  select(-n) %>%
  rename(activity = log10P)


#' load CPTAC TF activities
tf_activity <- read_csv(file = "./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "activity") %>%
  rename(tf=X1) %>%
  select(sample, tf, activity)


#' load metadata
cptac_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")
cptac_samples <- getSamples(cptac_samples, c("cnv", "phosphorylation", "mRNA")) %>%
  select(sample, batch, tissue) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(-tissue)


#' kinase to CNV correlation
kin_cnv <- kin_activity %>%
  inner_join(cnv, by = c("sample", "kinase" = "gene")) %>%
  inner_join(cptac_samples, by = "sample")

#' all samples
kin_cnv_corr_all <- kin_cnv %>%
  group_by(kinase) %>%
  nest() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$cnv, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(id = kinase, cnv_cor_value = estimate) %>%
  mutate(protein_type = "kinase", cor_type = "cnv_kinase")

#' by batch
kin_cnv_corr_by_batch <- kin_cnv %>%
  group_by(batch, kinase) %>%
  nest() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$cnv, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(batch, id = kinase, cnv_cor_value = estimate) %>%
  mutate(protein_type = "kinase", cor_type = "cnv_kinase")


#' TF to CNV correlation
tf_cnv <- tf_activity %>%
  inner_join(cnv, by = c("sample", "tf" = "gene")) %>%
  inner_join(cptac_samples, by = "sample")

#' all samples
tf_cnv_corr_all <- tf_cnv %>%
  group_by(tf) %>%
  nest() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$cnv, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(id = tf, cnv_cor_value = estimate) %>%
  mutate(protein_type = "TF", cor_type = "cnv_tf")

#' by batch
tf_cnv_corr_by_batch <- tf_cnv %>%
  group_by(batch, tf) %>%
  nest() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$activity, .x$cnv, method = "pearson")))) %>%
  select(-data) %>%
  unnest(cor) %>%
  select(batch, id = tf, cnv_cor_value = estimate) %>%
  mutate(protein_type = "TF", cor_type = "cnv_tf")


#' join kinase and TF correlations (all samples)
kin_tf_cnv_cor_all <- kin_cnv_corr_all %>%
  bind_rows(tf_cnv_corr_all)

#' plot correlation distributions
kin_tf_cnv_cor_all_plot <- kin_tf_cnv_cor_all %>%
  ggplot(mapping = aes(x = protein_type, y = cnv_cor_value, fill = protein_type)) +
  geom_boxplot(size = 2, outlier.size = 5, show.legend = F, alpha = 0.7) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 1),
    plot.title = element_text(size = 28, colour = "black", hjust = 0.5),
    axis.title.x = element_text(size = 28, colour = "black"),
    axis.text.x = element_text(size = 26, colour = "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    legend.position = "none") +
  scale_x_discrete(labels = c(kinase = "Kinases", TF = "TFs")) +
  scale_fill_viridis(discrete = T) +
  labs(y = "Pearson's r", title =  "CNV vs protein activity")

#+ fig.width=14, fig.height=5
kin_tf_cnv_cor_all_plot

ggsave(filename = "kinase_tf_cnv_correlation_boxplot.png", plot = kin_tf_cnv_cor_all_plot, path = "./output/plots/cnv_kinase_tf_act_correlation/", height = 5, width = 14)
ggsave(filename = "kinase_tf_cnv_correlation_boxplot.pdf", plot = kin_tf_cnv_cor_all_plot, path = "./output/plots/cnv_kinase_tf_act_correlation/", height = 5, width = 14)


#' join kinase and TF correlations (by batch)
kin_tf_cnv_cor_by_batch <- kin_cnv_corr_by_batch %>%
  bind_rows(tf_cnv_corr_by_batch)

#' plot correlation distributions
N <- kin_tf_cnv_cor_by_batch %>%
  group_by(batch, protein_type) %>%
  tally() %>%
  ungroup() %>%
  mutate(batch = fct_recode(batch, Brain = "cbttc", `Colorectal (ccle)` = "ccle-colorectal", `Colon`= "colon-opportunities", `Kidney` = "discovery-ccrcc", `Uterus` = "discovery-ucec", `Liver` = "hcc-proteogenomics", `Breast (tcga)` = "tcga-brca", `Ovary` = "tcga-ov")) %>%
  mutate(batch = factor(x = as.character(batch), levels = sort(unique(as.character(batch)))))

kin_tf_cnv_cor_by_batch_plot <- kin_tf_cnv_cor_by_batch %>%
  mutate(batch = fct_recode(batch, Brain = "cbttc", `Colorectal (ccle)` = "ccle-colorectal", `Colon`= "colon-opportunities", `Kidney` = "discovery-ccrcc", `Uterus` = "discovery-ucec", `Liver` = "hcc-proteogenomics", `Breast (tcga)` = "tcga-brca", `Ovary` = "tcga-ov")) %>%
  mutate(batch = factor(x = as.character(batch), levels = sort(unique(as.character(batch))))) %>%
  ggplot(mapping = aes(x = protein_type, y = cnv_cor_value, fill = protein_type)) +
  geom_boxplot(size = 1, outlier.size = 2, show.legend = F, alpha = 0.7) +
  geom_text(data = N, mapping = aes(x = rep(c(1.3, 2.3),8), y = -0.6, label = n), size = 6) +
  theme_classic() +
  coord_flip() +
  facet_wrap(facets = vars(batch), scales = "fixed", nrow = 2) +
  theme(
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 1),
    plot.title = element_text(size = 28, colour = "black", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(size = 24, colour = "black"),
    axis.title.x = element_text(size = 28, colour = "black"),
    axis.text.x = element_text(size = 26, colour = "black"),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 24, colour = "black"),
    legend.position = "none") +
  scale_x_discrete(labels = c(kinase = "Kinases", TF = "TFs")) +
  scale_fill_viridis(discrete = T) +
  labs(y = "Pearson's r", title =  "CNV vs protein activity")

#+ fig.width=14, fig.height=5
kin_tf_cnv_cor_by_batch_plot

ggsave(filename = "kinase_tf_cnv_correlation_by_batch_boxplot.png", plot = kin_tf_cnv_cor_by_batch_plot, path = "./output/plots/cnv_kinase_tf_act_correlation/", height = 5, width = 14)
ggsave(filename = "kinase_tf_cnv_correlation_by_batch_boxplot.pdf", plot = kin_tf_cnv_cor_by_batch_plot, path = "./output/plots/cnv_kinase_tf_act_correlation/", height = 5, width = 14)
