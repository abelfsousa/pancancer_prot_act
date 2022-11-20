#' ---
#' title: "Kinase-substrate lists"
#' author: "Abel Sousa"
#' ---


#+ setup, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)


#' load R packages
library(tidyverse)
library(VennDiagram)
library(RColorBrewer)
library(viridis)
library(cowplot)
source("./src/utils/psite_prot_seq_match.R")


#' load kinase substrate lists
in_vv_pmapper <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz")


#' number of phosphosites per data source
sources_bp1 <- in_vv_pmapper %>%
  mutate(source_type = fct_rev(fct_infreq(f=source_type))) %>%
  mutate(source = fct_rev(fct_infreq(f=source))) %>%
  ggplot() +
  geom_bar(mapping = aes(x = source_type, fill = source), position = "dodge") +
  theme_classic() +
  coord_flip() +
  theme(
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14, colour = "black"),
    #legend.title = element_text(size = 8, colour = "black"),
    #legend.text = element_text(size = 6, colour = "black"),
    legend.position = "bottom") +
  scale_fill_viridis(discrete = T, 
                     name = "Source",
                     labels = c("in_vitro" = "in vitro", "in_vivo" = "in vivo", "database" = "Database", "text-mining" = "Text mining", "bel" = "BEL Large Corpus", "reactome" = "Reactome", "pid" = "NCI-PID", "rlimsp" = "RLIMS-P", "reach" = "REACH", "sparser" = "Sparser", "signor" = "Signor", "psp" = "PhosphoSitePlus")) +
  scale_x_discrete(name = "Source type", labels = c("in_vitro" = "in vitro", "in_vivo" = "in vivo", "database" = "Database", "text-mining" = "Text mining")) +
  scale_y_continuous(name = "Number of kinase-substrate associations")

sources_bp2 <- in_vv_pmapper %>%
  mutate(source_type = fct_rev(fct_infreq(f=source_type))) %>%
  mutate(source = fct_rev(fct_infreq(f=source))) %>%
  ggplot() +
  geom_bar(mapping = aes(x = source_type, fill = source), position = "dodge") +
  theme_classic() +
  coord_flip() +
  theme(
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_viridis(discrete = T, guide = F) +
  scale_x_discrete(name = "", labels = c("in_vitro" = "in vitro", "in_vivo" = "in vivo", "database" = "Database", "text-mining" = "Text mining")) +
  scale_y_continuous(name = "", limits = c(NA, 18000))

sources_bp <- ggdraw(sources_bp1) +
  draw_plot(sources_bp2, .4, .36, .6, .45)


#+ fig.width=6, fig.height=5
sources_bp

ggsave(filename="kinase_substrate_list1.png", plot = sources_bp, path = "./output/plots/kinase_substrate_pairs/", width = 7, height = 4)
ggsave(filename="kinase_substrate_list1.pdf", plot = sources_bp, path = "./output/plots/kinase_substrate_pairs/", width = 7, height = 4)


#' shrink phosphosite table so that each phosphosite only has one row\
#' collapse multiple sources per phosphosite (use comma)\
#' number of phosphosites by number of different sources
sources_n_bp1 <- in_vv_pmapper %>%
  group_by(kinase, substrate, pair, position, residue) %>%
  summarise(source = list(source), source_type = list(source_type), source_n = n()) %>%
  ungroup() %>%
  mutate(source = map_chr(.x = source, .f = ~ str_c(sort(.x), collapse = ","))) %>%
  mutate(source_type = map_chr(.x = source_type, .f = ~ str_c(sort(unique(.x)), collapse = ","))) %>%
  mutate(source_n = as.character(source_n)) %>%
  mutate(source_n = fct_rev(fct_infreq(f=source_n))) %>%
  ggplot() +
  geom_bar(mapping = aes(x = source_n, y = ..count..)) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.text = element_text(size = 10, colour = "black")
  ) +
  scale_x_discrete(name = "Number of different sources") +
  scale_y_continuous(name = "Number of kinase-substrate associations")

sources_n_bp2 <- in_vv_pmapper %>%
  group_by(kinase, substrate, pair, position, residue) %>%
  summarise(source = list(source), source_type = list(source_type), source_n = n()) %>%
  ungroup() %>%
  mutate(source = map_chr(.x = source, .f = ~ str_c(sort(.x), collapse = ","))) %>%
  mutate(source_type = map_chr(.x = source_type, .f = ~ str_c(sort(unique(.x)), collapse = ","))) %>%
  mutate(source_n = as.character(source_n)) %>%
  mutate(source_n = fct_rev(fct_infreq(f=source_n))) %>%
  ggplot() +
  geom_bar(mapping = aes(x = source_n, y = ..count..)) +
  theme_classic() +
  coord_flip() +
  theme(
    axis.text = element_text(size = 10, colour = "black")
  ) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "", limits = c(NA, 4000))

sources_n_bp <- ggdraw(sources_n_bp1) +
  draw_plot(sources_n_bp2, .39, .17, .6, .6)


#+ fig.width=5, fig.height=3
sources_n_bp

ggsave(filename="kinase_substrate_list2.png", plot = sources_n_bp, path = "./output/plots/kinase_substrate_pairs/", width = 5, height = 3)



#' number of kinases/substrates per data source
kn_s1 <- in_vv_pmapper %>%
  select(source_type, kinase) %>%
  distinct() %>%
  group_by(source_type) %>%
  tally() %>%
  ungroup() %>%
  mutate(data = "kinases")

kn_s2 <- in_vv_pmapper %>%
  select(source_type, substrate) %>%
  distinct() %>%
  group_by(source_type) %>%
  tally() %>%
  ungroup() %>%
  mutate(data = "substrates")

kn_s <- bind_rows(
  kn_s1,
  kn_s2) %>%
  mutate(source_type = fct_reorder(source_type, n, mean)) %>%
  ggplot() +
  geom_bar(mapping = aes(x = source_type, y = n, fill = data), stat = "identity", position = "dodge") +
  theme_classic() +
  coord_flip() +
  theme(
    axis.title = element_text(size = 14, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    legend.text = element_text(size = 12, colour = "black"),
    legend.title = element_text(size = 14, colour = "black"),
    legend.position = "bottom") +
  scale_x_discrete(name = "Source type", labels = c("in_vitro" = "in vitro", "in_vivo" = "in vivo", "database" = "Database", "text-mining" = "Text mining")) +
  scale_y_continuous(name = "Number of kinases/substrates") + 
  scale_fill_viridis(discrete = T, name = "Feature", labels = c("kinases" = "kinases", "substrates" = "Substrates"))

#+ fig.width=4, fig.height=2
kn_s

ggsave(filename="kinase_substrate_number.png", plot = kn_s, path = "./output/plots/kinase_substrate_pairs/", width = 5, height = 4)
ggsave(filename="kinase_substrate_number.pdf", plot = kn_s, path = "./output/plots/kinase_substrate_pairs/", width = 5, height = 4)


#' overlap of kinase-substrate associations between data sources
ks_overlap <- in_vv_pmapper %>%
  filter(source_type == "database") %>%
  select(-source_type, -pair) %>%
  unite("psite", kinase, substrate, position, residue) %>%
  group_by(source) %>%
  summarise(psite = list(psite)) %>%
  ungroup()

colorv <- brewer.pal(5, "Pastel2")

venn <- venn.diagram(
  x=ks_overlap$psite,
  category.names = ks_overlap$source,
  filename = "./output/plots/kinase_substrate_pairs/venn_diagram_database_lists.png",
  imagetype = "png",
  fill = colorv)

venn <- venn.diagram(
  x=ks_overlap$psite,
  category.names = ks_overlap$source,
  filename = NULL,
  fill = colorv)

grid.newpage()
grid::grid.draw(venn)


ks_overlap <- in_vv_pmapper %>%
  filter(source_type == "text-mining") %>%
  select(-source_type, -pair) %>%
  unite("psite", kinase, substrate, position, residue) %>%
  group_by(source) %>%
  summarise(psite = list(psite)) %>%
  ungroup()

colorv <- brewer.pal(3, "Pastel2")

venn <- venn.diagram(
  x=ks_overlap$psite,
  category.names = ks_overlap$source,
  filename = "./output/plots/kinase_substrate_pairs/venn_diagram_text_mining_lists.png",
  imagetype = "png",
  fill = colorv)

venn <- venn.diagram(
  x=ks_overlap$psite,
  category.names = ks_overlap$source,
  filename = NULL,
  fill = colorv)

grid.newpage()
grid::grid.draw(venn)


ks_overlap <- in_vv_pmapper %>%
  filter(str_detect(source_type, "in_", negate = T)) %>%
  select(-source, -pair) %>%
  distinct() %>%
  unite("psite", kinase, substrate, position, residue) %>%
  group_by(source_type) %>%
  summarise(psite = list(psite)) %>%
  ungroup()

colorv <- brewer.pal(3, "Pastel2")[-1]

venn <- venn.diagram(
  x=ks_overlap$psite,
  category.names = ks_overlap$source_type,
  filename = "./output/plots/kinase_substrate_pairs/venn_diagram_db_tm_lists.png",
  imagetype = "png",
  fill = colorv)

venn <- venn.diagram(
  x=ks_overlap$psite,
  category.names = ks_overlap$source_type,
  filename = NULL,
  fill = colorv)

grid.newpage()
grid::grid.draw(venn)


ks_overlap <- in_vv_pmapper %>%
  select(-pair, -source) %>%
  distinct() %>%
  unite("psite", kinase, substrate, position, residue) %>%
  group_by(source_type) %>%
  summarise(psite = list(psite)) %>%
  ungroup()

colorv <- brewer.pal(4, "Pastel2")

venn <- venn.diagram(
  x=ks_overlap$psite,
  category.names = ks_overlap$source_type,
  filename = "./output/plots/kinase_substrate_pairs/venn_diagram_db_tm_invivo_invitro_lists.png",
  imagetype = "png",
  fill = colorv)

venn <- venn.diagram(
  x=ks_overlap$psite,
  category.names = ks_overlap$source_type,
  filename = NULL,
  fill = colorv)

grid.newpage()
grid::grid.draw(venn)


files_to_delete <- str_subset(list.files("./output/plots/kinase_substrate_pairs/"), "[0-9]+.log")
file.remove(str_c("./output/plots/kinase_substrate_pairs/", files_to_delete))
