library(RColorBrewer)
library(tidyverse)
library(viridis)
library(cowplot)

source("./src/utils/getSamples.R")


#' load cancer samples for each data type

# mutation
mut <- read_tsv(file = "./output/files/mutations_samples.txt")

# copy-number variation
cnv <- read_tsv(file = "./output/files/cnv_samples.txt")

# trancriptomics
rna <- read_tsv(file = "./output/files/transcriptomics_samples.txt")

# proteomics
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt")

# phopshoproteomics
phospho <- read_tsv(file = "./output/files/phosphoproteomics_samples.txt")

# clinical
clinical <- read_tsv(file = "./output/files/clinical_samples.txt")


#' load multi-omics cancer dataset
multi_omics <- read_tsv(file = "./output/files/all_samples_annotation.txt")


#' make new annotation
cell_lines <- proteomics %>%
  select(sample, study=batch) %>%
  filter(str_detect(study, "cell-lines")) %>%
  mutate(study = str_replace(study, "cell-lines-", ""))

new_annotation <- multi_omics %>%
  select(data, batch, tissue, sample) %>%
  mutate(batch = str_c(batch, tissue, sep = "-")) %>%
  mutate(batch = str_replace(batch, "tcga-brca-breast", "CPTAC-breast\n(TCGA)")) %>%
  mutate(batch = str_replace(batch, "tcga-coread-colorectal", "CPTAC-colorectal\n(TCGA)")) %>%
  mutate(batch = str_replace(batch, "tcga-ov-ovary", "CPTAC-ovarian\n(TCGA)")) %>%
  mutate(batch = str_replace(batch, "cbttc-brain", "CPTAC-brain")) %>%
  mutate(batch = str_replace(batch, "discovery-ccrcc-kidney", "CPTAC-kidney\n(discovery)")) %>%
  mutate(batch = str_replace(batch, "discovery-luad-lung", "CPTAC-lung\n(discovery)")) %>%
  mutate(batch = str_replace(batch, "discovery-ucec-uterus", "CPTAC-uterus\n(discovery)")) %>%
  mutate(batch = str_replace(batch, "colon-opportunities-colorectal", "CPTAC-colon\n(confirmatory)")) %>%
  mutate(batch = str_replace(batch, "eogc-proteogenomics-stomach", "CPTAC-stomach")) %>%
  mutate(batch = str_replace(batch, "hcc-proteogenomics-liver", "CPTAC-liver")) %>%
  mutate(batch = str_replace(batch, "ccle-colorectal", "CCLE-colorectal")) %>%
  mutate(batch = str_replace(batch, "ccle-breast", "CCLE-breast"))

overlap <- new_annotation %>%
  group_by(data, batch, tissue) %>%
  summarise(sample = list(sample)) %>%
  ungroup() %>%
  group_by(batch, tissue) %>%
  summarise(sample = list(sample)) %>%
  ungroup() %>%
  mutate(overlap = map(.x=sample, .f = ~ reduce(.x=.x, .f=intersect))) %>%
  select(-sample) %>%
  unnest() %>%
  rename(sample = overlap) %>%
  mutate(data = "overlap")

new_annotation <- new_annotation %>%
  bind_rows(overlap) %>%
  left_join(cell_lines, by = "sample")
write_tsv(new_annotation, "./output/files/all_samples_annotation2.txt")


# circular barplots
circular_barplot1 <- new_annotation %>%
  bind_rows(tibble(data = NA, batch = as.character(1:4), tissue = NA, sample = NA, study = NA)) %>%
  mutate(data = fct_rev(fct_relevel(data, "mutation", "cnv", "mRNA", "protein", "phosphorylation", "clinical", "overlap"))) %>%
  mutate(batch = fct_rev(fct_infreq(batch))) %>%
  group_by(data, batch, tissue) %>%
  tally() %>%
  ungroup() %>%
  mutate(n = n*-1) %>%
  ggplot(mapping = aes(x = batch, y = n, fill = data)) +
  geom_col(position = "dodge") +
  geom_segment(data = tibble(x = 4.5, y = seq(-100, -1100, by=-200), xend = 16, yend=y), aes(x = x, y = y, xend = xend, yend = yend), color = "grey", alpha = 1, inherit.aes = FALSE, size = 0.3) +
  geom_text(data = tibble(x = 4.4, y = seq(-100, -1100, by=-200), label=y*-1), aes(x = x, y = y, label = label), angle = 0, vjust = c(0.6, 0.4, 0.3, 0.2, 0.1, 0), size = 2.5, color = "black", inherit.aes = FALSE) +
  annotate("text", x = 4, y = -600, label = "Number of samples", size = 3) +
  geom_text(
    data = tibble(x = rev(5:16), y = 120, label = c("CPTAC-breast\n(TCGA)", "CPTAC-brain", "CPTAC-colorectal\n(TCGA)", "CPTAC-ovarian\n(TCGA)", "CPTAC-liver", "CPTAC-kidney\n(discovery)", "CPTAC-colon\n(confirmatory)", "CPTAC-uterus\n(discovery)", "CPTAC-lung\n(discovery)", "CPTAC-stomach", "CCLE-colorectal", "CCLE-breast")),
    mapping = aes(x = x, y = y, label = label),
    angle = c(-10, -35, -55, -75, -100, -125, -145, -170, -192, -212, -237, -257),
    size = 2.5,
    inherit.aes = F) +
  scale_fill_viridis(discrete = T, na.translate = F, name = "Data") +
  scale_y_continuous(limits = c(-1500, 120)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "bottom") +
  coord_polar(direction = -1)

ggsave(filename="circular_barplot1.pdf", plot = circular_barplot1, path = "./output/plots/data_description/", width=7, height=6)
unlink("circular_barplot1.pdf")


circular_barplot2 <- new_annotation %>%
  filter(data == "overlap") %>%
  group_by(data, batch, tissue) %>%
  tally() %>%
  ungroup() %>%
  bind_rows(tibble(data = NA, batch = as.character(c(1:50)), tissue = NA, study = NA, n = NA)) %>%
  mutate(batch = fct_reorder(batch, n, .fun = function(x) x)) %>%
  mutate(n = n*-1) %>%
  ggplot(mapping = aes(x = batch, y = n, fill = data)) +
  geom_col(position = "dodge") +
  geom_segment(data = tibble(x = 12.8, y = seq(-50, -150, by=-50), xend = 0.5, yend=y), aes(x = x, y = y, xend = xend, yend = yend), color = "grey", alpha = 0.8, inherit.aes = FALSE) +
  geom_text(data = tibble(x = 13, y = c(-50,-100,-150), label=y*-1), aes(x = x, y = y, label = label), angle = 0, hjust = -0.1, color = "black", inherit.aes = FALSE) +
  annotate("text", x = 6, y = 40, label = "Overlapping samples", angle = 40) +
  geom_text(
    data = tibble(x = rev(1:12), y = 15, label = c("CPTAC\nliver", "CPTAC\nbrain", "CPTAC\nlung", "CPTAC\nkidney", "CPTAC\ncolon", "CPTAC\nuterus", "CPTAC\nstomach", "CPTAC\nbreast", "CPTAC\ncolo\nrectal", "CCLE\ncolo\nrectal", "CCLE\nbreast", "CPTAC\novarian")),
    mapping = aes(x = x, y = y, label = label),
    angle = c(7, 10, 16, 24, 32, 38, 42, 48, 52, 55, 65, 70),
    size = 1.9,
    inherit.aes = F) +
  scale_fill_viridis(discrete = T) +
  scale_y_continuous(limits = c(-200, 40)) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "none") +
  coord_polar(direction = 1, start = 11.3)

ggsave(filename="circular_barplot2.pdf", plot = circular_barplot2, path = "./output/plots/data_description/", width=8, height=8)
unlink("circular_barplot2.pdf")


# regular barplots
regular_barplot1 <- new_annotation %>%
  mutate(data = fct_rev(fct_infreq(data))) %>%
  mutate(data = fct_relevel(data, "overlap")) %>%
  mutate(data = fct_recode(data, Overlap = "overlap", Phosphorylation = "phosphorylation", Protein = "protein", Clinical = "clinical", Mutation = "mutation", CNV = "cnv")) %>%
  mutate(batch = fct_rev(fct_infreq(batch))) %>%
  ggplot(mapping = aes(x = batch, fill = data)) +
  geom_bar(position = "dodge") +
  #scale_fill_brewer(type = "seq", palette = "PuBuGn", name = "Data", direction = 1) +
  scale_fill_viridis(discrete = T, name = "Data type", direction = 1) +
  scale_y_continuous(limits = c(NA, 1150), name = "Number of samples") +
  theme_classic() +
  theme(
    plot.margin = unit(c(0, 1, 2, 0), "cm"),
    axis.title.y = element_blank(),
    axis.title.x = element_text(colour = "black", size = 18),
    axis.text.y = element_text(colour = "black", size = 14),
    axis.text.x = element_text(colour = "black", size = 16),
    legend.title = element_text(colour = "black", size = 18),
    legend.text = element_text(colour = "black", size = 16),
    legend.direction = "horizontal",
    legend.position = c(0.4,0),
    legend.justification = c(0.5,2)) +
  coord_flip() +
  guides(fill = guide_legend(reverse=T))

regular_barplot2 <- new_annotation %>%
  mutate(batch = fct_rev(fct_infreq(batch))) %>%
  filter(data == "overlap") %>%
  mutate(batch = fct_recode(batch, `CPTAC-lung (discovery)` = "CPTAC-lung\n(discovery)", `CPTAC-uterus (discovery)` = "CPTAC-uterus\n(discovery)", `CPTAC-colon (confirmatory)` = "CPTAC-colon\n(confirmatory)", `CPTAC-kidney (discovery)` = "CPTAC-kidney\n(discovery)", `CPTAC-ovarian (TCGA)` = "CPTAC-ovarian\n(TCGA)", `CPTAC-colorectal (TCGA)` = "CPTAC-colorectal\n(TCGA)", `CPTAC-breast (TCGA)` = "CPTAC-breast\n(TCGA)")) %>%
  ggplot(mapping = aes(x = batch)) +
  geom_bar(fill = viridis_pal()(7)[1]) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(colour = "black", size = 14),
    axis.text.y = element_text(colour = "black", size = 12),
    plot.title = element_text(colour = "black", size = 16)) +
  coord_flip() + 
  labs(title = "Multi-omics overlap")

regular_barplot <- ggdraw(regular_barplot1) +
  draw_plot(regular_barplot2, 0.375, 0.2, 0.58, 0.5)

ggsave(filename="regular_barplot.png", plot = regular_barplot, path = "./output/plots/data_description/", width=8, height=8)
ggsave(filename="regular_barplot.pdf", plot = regular_barplot, path = "./output/plots/data_description/", width=8, height=8)
unlink("regular_barplot.png")
unlink("regular_barplot.pdf")


regular_barplot3 <- new_annotation %>%
  filter(data == "overlap") %>%
  mutate(batch = fct_rev(fct_infreq(batch))) %>%
  mutate(batch = fct_recode(batch, `CPTAC-lung (discovery)` = "CPTAC-lung\n(discovery)", `CPTAC-uterus (discovery)` = "CPTAC-uterus\n(discovery)", `CPTAC-colon (confirmatory)` = "CPTAC-colon\n(confirmatory)", `CPTAC-kidney (discovery)` = "CPTAC-kidney\n(discovery)", `CPTAC-ovarian (TCGA)` = "CPTAC-ovarian\n(TCGA)", `CPTAC-colorectal (TCGA)` = "CPTAC-colorectal\n(TCGA)", `CPTAC-breast (TCGA)` = "CPTAC-breast\n(TCGA)")) %>%
  ggplot(mapping = aes(x = batch)) +
  #geom_bar(fill = viridis_pal()(7)[1]) +
  geom_bar(fill = "#0570b0", color = "black") +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 20),
    axis.text.x = element_text(colour = "black", size = 18),
    axis.text.y = element_text(colour = "black", size = 18)) +
  coord_flip() +
  scale_y_continuous(limits = c(NA, 160), breaks = seq(0, 160, 25)) +
  labs(x = "Study", y = "Number of samples")
ggsave(filename="regular_barplot3.png", plot = regular_barplot3, path = "./output/plots/data_description/", width=8, height=8)
ggsave(filename="regular_barplot3.pdf", plot = regular_barplot3, path = "./output/plots/data_description/", width=8, height=8)


regular_barplot4 <- new_annotation %>%
  mutate(data = fct_rev(fct_infreq(data))) %>%
  mutate(data = fct_relevel(data, "overlap")) %>%
  mutate(data = fct_recode(data, Overlap = "overlap", Phosphorylation = "phosphorylation", Protein = "protein", Clinical = "clinical", Mutation = "mutation", CNV = "cnv")) %>%
  mutate(batch = fct_rev(fct_infreq(batch))) %>%
  mutate(batch = fct_recode(batch, `CPTAC-lung (discovery)` = "CPTAC-lung\n(discovery)", `CPTAC-uterus (discovery)` = "CPTAC-uterus\n(discovery)", `CPTAC-colon (confirmatory)` = "CPTAC-colon\n(confirmatory)", `CPTAC-kidney (discovery)` = "CPTAC-kidney\n(discovery)", `CPTAC-ovarian (TCGA)` = "CPTAC-ovarian\n(TCGA)", `CPTAC-colorectal (TCGA)` = "CPTAC-colorectal\n(TCGA)", `CPTAC-breast (TCGA)` = "CPTAC-breast\n(TCGA)")) %>%
  ggplot(mapping = aes(x = batch, fill = data)) +
  geom_bar(position = "dodge", color = "black") +
  #scale_fill_brewer(type = "div", palette = "RdYlBu", direction = -1) +
  scale_fill_viridis(discrete = T, direction = 1) +
  scale_y_continuous(limits = c(NA, 1200), breaks = seq(0, 1200, 100)) +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 30),
    axis.text = element_text(colour = "black", size = 28),
    legend.position = "bottom",
    legend.title = element_text(colour = "black", size = 30),
    legend.text = element_text(colour = "black", size = 28),
    legend.direction = "horizontal") +
  coord_flip() +
  labs(x = "Study", y = "Number of samples", fill = "Data type") +
  guides(fill = guide_legend(reverse=T))
ggsave(filename="regular_barplot4.png", plot = regular_barplot4, path = "./output/plots/data_description/", width=20, height=14)
ggsave(filename="regular_barplot4.pdf", plot = regular_barplot4, path = "./output/plots/data_description/", width=20, height=14)
