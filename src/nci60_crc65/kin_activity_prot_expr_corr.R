#' ---
#' title: "Correlation of kinase activity to protein expression - NCI60/CRC65 data"
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


#' load kinase activity inference data from NCI60/CRC65 dataset\
#' select quantifications with more than 3 substrates
nci60_crc65_KA <- read_tsv(file = "./output/files/nci60_crc65_kinase_activities.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= 3) %>%
  select(-source_type, -n) %>%
  filter(!(batch == "CRC65" & sample %in% shared_cell_lines)) %>%
  select(batch, sample, kinase, kin_activity=log10P)


#' load protein abundance
protein <- read_tsv(file = "./output/files/nci60_crc65_protein_log2fc.txt.gz") %>%
  filter(!(batch == "CRC65" & sample %in% shared_cell_lines))


#' # scatterplot - correlation of kinase activity to protein expression
ka_prot_expr <- nci60_crc65_KA %>%
  inner_join(protein, by = c("batch", "kinase" = "gene", "sample")) %>%
  filter(!is.na(log2fc)) %>%
  ggplot(mapping = aes(x = kin_activity, y = log2fc)) +
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
  labs(x = "Kinase activity score (log10)", y = "Protein expression (log2fc)")

#+ fig.width=4, fig.height=4
ka_prot_expr

ggsave(filename = "KA_prot_expr_cor_scatter.png", plot = ka_prot_expr, path = "./output/plots/nci60_crc65/", height = 4, width = 4)
unlink("KA_prot_expr_cor_scatter.png")


#' # scatterplot - correlation of kinase activity to protein expression by batch
ka_prot_expr <- nci60_crc65_KA %>%
  inner_join(protein, by = c("batch", "kinase" = "gene", "sample")) %>%
  filter(!is.na(log2fc)) %>%
  ggplot(mapping = aes(x = kin_activity, y = log2fc)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(~ batch, scales = "free") +
  theme_classic() +
  stat_cor(size = 3) +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 10),
    aspect.ratio = 1) +
  labs(x = "Kinase activity score (log10)", y = "Protein expression (log2fc)")

#+ fig.width=6, fig.height=6
ka_prot_expr

ggsave(filename = "KA_prot_expr_cor_scatter_batch.png", plot = ka_prot_expr, path = "./output/plots/nci60_crc65/", height = 6, width = 6)
unlink("KA_prot_expr_cor_scatter_batch")


#' # boxplot - distribution of kinase activity to protein expression correlation across samples by kinase
ka_prot_expr <- nci60_crc65_KA %>%
  inner_join(protein, by = c("batch", "kinase" = "gene", "sample")) %>%
  filter(!is.na(log2fc)) %>%
  group_by(kinase) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > 10) %>%
  select(-n) %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$kin_activity, .x$log2fc)))) %>%
  select(-data) %>%
  unnest(cor) %>%
  ggplot(mapping = aes(x = factor(0), y = estimate)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "KA vs protein expression", y = "Pearson's r")

#+ fig.width=4, fig.height=4
ka_prot_expr

ggsave(filename = "KA_prot_expr_cor_boxpl.png", plot = ka_prot_expr, path = "./output/plots/nci60_crc65/", height = 4, width = 4)
unlink("KA_prot_expr_cor_boxpl.png")



#' # boxplot - distribution of kinase activity to protein expression correlation (by batch) across samples by kinase
ka_prot_expr <- nci60_crc65_KA %>%
  inner_join(protein, by = c("batch", "kinase" = "gene", "sample")) %>%
  filter(!is.na(log2fc)) %>%
  group_by(batch, kinase) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > 10) %>%
  select(-n) %>%
  group_by(batch, kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$kin_activity, .x$log2fc)))) %>%
  select(-data) %>%
  unnest(cor) %>%
  ggplot(mapping = aes(x = batch, y = estimate)) +
  geom_boxplot() +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8)) +
  labs(x = "KA vs protein expression", y = "Pearson's r")

#+ fig.width=4, fig.height=4
ka_prot_expr

ggsave(filename = "KA_prot_expr_cor_boxpl_batch.png", plot = ka_prot_expr, path = "./output/plots/nci60_crc65/", height = 4, width = 4)
unlink("KA_prot_expr_cor_boxpl_batch")

