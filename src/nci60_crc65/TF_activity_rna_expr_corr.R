#' ---
#' title: "Correlation of TF activity to rna expression - NCI60/CRC65 data"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, tidy.opts=list(width.cutoff=40), tidy=TRUE)


#' load R packages
library(tidyverse)
library(ggpubr)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load NCI60/CRC65 cell lines
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


#' load TF activity inference data from NCI60/CRC65 dataset
crc65_TF <- read_csv(file = "./data/dorothea/nci60_crc65/TF_activity_log2FC_allsamples_CRC65.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "tf_activity") %>%
  rename(tf=X1) %>%
  filter(!sample %in% shared_cell_lines) %>%
  select(sample, tf, tf_activity) %>%
  mutate(batch = "CRC65")

nci60_TF <- read_csv(file = "./data/dorothea/nci60_crc65/TF_activity_log2FC_allsamples_NCI60.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "tf_activity") %>%
  rename(tf=X1) %>%
  select(sample, tf, tf_activity) %>%
  mutate(batch = "NCI60")

nci60_crc65_TF <- bind_rows(crc65_TF, nci60_TF)


#' load rna abundance
rna <- read_tsv(file = "./output/files/nci60_crc65_rna_log2fc.txt.gz") %>%
  filter(!(batch == "CRC65" & sample %in% shared_cell_lines))


#' # scatterplot - correlation of kinase activity to protein expression by sample
tf_rna_expr <- nci60_crc65_TF %>%
  inner_join(rna, by = c("batch", "tf" = "gene", "sample")) %>%
  filter(!is.na(log2fc)) %>%
  ggplot(mapping = aes(x = tf_activity, y = log2fc)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  theme_classic() +
  stat_cor(size = 3) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10),
    aspect.ratio = 1) +
  labs(x = "TF activity estimate", y = "mRNA expression (log2fc)")


#+ fig.width=4, fig.height=4
tf_rna_expr

ggsave(filename = "TFact_rna_expr_scatterplot.png", plot = tf_rna_expr, path = "./output/plots/nci60_crc65/", height = 4, width = 4)
unlink("TFact_rna_expr_scatterplot.png")



#' # boxplot - distribution of TF activity to RNA expression correlation across samples by TF
tf_rna_expr <- nci60_crc65_TF %>%
  inner_join(rna, by = c("batch", "tf" = "gene", "sample")) %>%
  filter(!is.na(log2fc)) %>%
  group_by(tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$tf_activity, .x$log2fc)))) %>%
  select(-data) %>%
  unnest(cor) %>%
  ggplot(mapping = aes(x = factor(0), y = estimate)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "", y = "Pearson's r (TF activity vs RNA expression)")


#+ fig.width=4, fig.height=4
tf_rna_expr

ggsave(filename = "TFact_rna_expr_cor_boxpl.png", plot = tf_rna_expr, path = "./output/plots/nci60_crc65/", height = 4, width = 4)
unlink("TFact_rna_expr_cor_boxpl.png")

