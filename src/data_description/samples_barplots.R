#' ---
#' title: "Samples barplots"
#' author: "Abel Sousa"
#' ---

#+ setup, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)


#' load R packages
library(RColorBrewer)
library(tidyverse)
library(cowplot)
library(viridis)


#' # load samples annotation
all_samples_annotation <- read_tsv(file = "./output/files/all_samples_annotation.txt")


#' # mRNA and protein data
omics1 <- all_samples_annotation %>%
  filter(data %in% c("protein", "mRNA", "clinical")) %>%
  rename(info = data) %>%
  group_by(info, tissue) %>%
  summarise(sample = list(sample)) %>%
  ungroup() %>%
  spread(key = "info", value = "sample") %>%
  mutate(overlap = pmap(.l = ., .f = ~ reduce(.x = compact(list(..2, ..3, ..4)), .f = ~ intersect(.x, .y)))) %>%
  #group_by(tissue) %>%
  #mutate(overlap = list(reduce(compact(list(mRNA[[1]], protein[[1]], clinical[[1]])), intersect))) %>%
  #ungroup() %>%
  gather(key = "data", value = "samples", -c(tissue)) %>%
  mutate(n = map_dbl(samples, length))

omics1_bp1 <- omics1 %>%
  mutate(tissue = fct_reorder(tissue, n, mean)) %>%
  mutate(data = fct_rev(fct_relevel(data, c("mRNA", "protein", "clinical", "overlap")))) %>%
  ggplot() +
  geom_bar(mapping = aes(x = tissue, y = n, fill = data), stat = "identity", position = "dodge", color = "black") +
  theme_classic() +
  coord_flip() +
  scale_fill_viridis(discrete = "T", option = "B") +
  theme(
    axis.title=element_text(colour="black", size=16),
    axis.text=element_text(colour="black", size=14),
    legend.title=element_text(colour="black", size=16),
    legend.text=element_text(colour="black", size=14)) +
  labs(x = "Tissue", y = "Number of samples", fill = "Data type") +
  guides(fill = guide_legend(reverse = TRUE))

omics1_bp2 <- omics1 %>%
  filter(data == "overlap") %>%
  mutate(tissue = fct_reorder(tissue, n, mean)) %>%
  ggplot() +
  geom_bar(mapping = aes(x = tissue, y = n), stat = "identity", position = "dodge", color = "grey", fill = "black") +
  theme_classic() +
  coord_flip() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12)) +
  labs(x = "", y = "", title = "Overlap")

omics1_bp12 <- ggdraw(omics1_bp1) +
  draw_plot(omics1_bp2, .42, .12, .4, .5)

#+ fig.width=9, fig.height=5
omics1_bp12

ggsave(filename="ProtRNA_number_samples_study.png", plot = omics1_bp12, path = "./output/plots/data_description/", width=9, height=5)
unlink("ProtRNA_number_samples_study.png")




#' # all omics
omics2 <- all_samples_annotation %>%
  rename(info = data) %>%
  group_by(info, tissue) %>%
  summarise(sample = list(sample)) %>%
  ungroup() %>%
  spread(key = "info", value = "sample") %>%
  mutate(overlap = pmap(.l = ., .f = ~ reduce(.x = compact(list(..2, ..3, ..4, ..5, ..6, ..7)), .f = ~ intersect(.x, .y)))) %>%
  #group_by(tissue) %>%
  #mutate(overlap = list(reduce(compact(list(clinical[[1]], cnv[[1]], mRNA[[1]], mutation[[1]], phosphorylation[[1]], protein[[1]])), intersect))) %>%
  #ungroup() %>%
  gather(key = "data", value = "samples", -c(tissue)) %>%
  mutate(n = map_dbl(samples, length))

omics2_bp1 <- omics2 %>%
  mutate(tissue = fct_reorder(tissue, n, mean)) %>%
  mutate(data = fct_rev(fct_relevel(data, c("mutation", "cnv", "mRNA", "protein", "phosphorylation", "clinical", "overlap")))) %>%
  ggplot() +
  geom_bar(mapping = aes(x = tissue, y = n, fill = data), stat = "identity", position = "dodge", color = "black") +
  theme_classic() +
  coord_flip() +
  scale_fill_viridis(discrete = "T", option = "B") +
  theme(
    axis.title=element_text(colour="black", size=16),
    axis.text=element_text(colour="black", size=14),
    legend.title=element_text(colour="black", size=16),
    legend.text=element_text(colour="black", size=14)) +
  labs(x = "Tissue", y = "Number of samples", fill = "Data type") +
  guides(fill = guide_legend(reverse = TRUE))

omics2_bp2 <- omics2 %>%
  filter(data == "overlap") %>%
  mutate(tissue = fct_reorder(tissue, n, mean)) %>%
  ggplot() +
  geom_bar(mapping = aes(x = tissue, y = n), stat = "identity", position = "dodge", color = "grey", fill = "black") +
  theme_classic() +
  coord_flip() +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12)) +
  labs(x = "", y = "", title = "Overlap")

omics2_bp12 <- ggdraw(omics2_bp1) +
  draw_plot(omics2_bp2, .36, .12, .4, .5)

#+ fig.width=9, fig.height=5
omics2_bp12

ggsave(filename="AllOmics_number_samples_study.png", plot = omics2_bp12, path = "./output/plots/data_description/", width=9, height=5)
unlink("AllOmics_number_samples_study.png")





#' # number of cancer studies per tissue
tissues_studies <- all_samples_annotation %>%
  filter(data == "protein") %>%
  select(tissue, batch) %>%
  distinct() %>%
  mutate(tissue = fct_rev(fct_infreq(tissue))) %>%
  ggplot() +
  geom_bar(mapping = aes(x = tissue, y = ..count.., fill = batch), position = "stack") +
  theme_classic() +
  coord_flip() +
  scale_fill_viridis(discrete = T) +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12)) +
  labs(x = "", y = "", title = "")

#+ fig.width=5, fig.height=5
tissues_studies

ggsave(filename="tissues_studies.png", plot = tissues_studies, path = "./output/plots/data_description/", width=5, height=5)
unlink("tissues_studies.png")




#' # number of cancer subtypes per tissue
tissues_cancers <- all_samples_annotation %>%
  filter(data == "protein") %>%
  select(tissue, cancer) %>%
  distinct() %>%
  mutate(tissue = fct_rev(fct_relevel(tissue, "colorectal", "breast", "brain", "kidney", "liver", "lung", "ovary", "stomach", "uterus"))) %>%
  ggplot() +
  geom_bar(mapping = aes(x = tissue, y = ..count.., fill = cancer), position = "stack") +
  theme_classic() +
  coord_flip() +
  scale_fill_viridis(discrete = T) +
  theme(
    axis.title=element_text(colour="black", size=15),
    axis.text=element_text(colour="black", size=12)) +
  labs(x = "", y = "", title = "")

#+ fig.width=5, fig.height=5
tissues_cancers

ggsave(filename="tissues_cancers.png", plot = tissues_cancers, path = "./output/plots/data_description/", width=5, height=5)
unlink("tissues_cancers.png")



