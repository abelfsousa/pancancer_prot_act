#' ---
#' title: "Comparison of the number of kinases regulated across tumour samples/perturbations"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=6, fig.height=6)


#' load R packages
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(viridis)
library(RColorBrewer)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' source R pipeline
source(file = "./src/utils/generalists_specialists_kinases.R")


set.seed(123)


#' load kinase activities from perturbations\
#' David/Danish's second compilation
KA_pertb <- read_tsv(file = "./output/files/KA_esetNR.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  select(-source_type)

#' load kinase activities from tumor samples
KA_tumors <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  select(-source_type)


#' load sample metadata
metadata <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
  filter(data == "phosphorylation") %>%
  select(-data)


#' kinase regulation across all tumour samples\
#' compared against cell lines perturbations
KA_reg <- gs_mapper(KA_tumors, KA_pertb)

KA_reg_plot <- KA_reg %>%
  ggplot() +
  theme_classic() +
  geom_point(mapping = aes(x = pP, y = pT, color = class), alpha = 0.7) +
  geom_smooth(mapping = aes(x = pP, y = pT), method = "lm") +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Regulated kinase") +
  geom_text_repel(data = KA_reg[KA_reg$class != "both", ], mapping = aes(x = pP, y = pT, label = kinase), box.padding = 0.3, point.padding = 0.3, size = 3, colour="black") +
  stat_cor(mapping = aes(x = pP, y = pT), method = "pearson", size = 3) +
  labs(x = "Perturbations (# conditions regulated)", y = "Tumour samples/CCLE (# samples regulated)")

#+ fig.width=8, fig.height=4
KA_reg_plot

ggsave(filename = "kin_activity_number_perturbations_tumours_correlation_all_samples.png", plot = KA_reg_plot, path = "./output/plots/generalists_specialists/", width = 8, height = 4)
ggsave(filename = "kin_activity_number_perturbations_tumours_correlation_all_samples.pdf", plot = KA_reg_plot, path = "./output/plots/generalists_specialists/", width = 8, height = 4)

KA_reg_plot <- KA_reg %>%
  mutate(class = fct_relevel(class, "tumour", "both", "perturbation")) %>%
  ggplot() +
  geom_point(mapping = aes(x = pP, y = pT, color = class, alpha = class, size = class)) +
  geom_smooth(mapping = aes(x = pP, y = pT), method = "lm", se = T, size = 0.5) +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Regulated kinase") +
  scale_alpha_manual(values = c(1,0.5,1), guide = F) +
  scale_size_manual(values = c(1.5,1,1.5), guide = F) +
  #geom_text_repel(data = KA_reg[KA_reg$class != "both", ], mapping = aes(x = pP, y = pT, label = kinase), box.padding = 0.3, point.padding = 0.3, size = 4, colour="black", seed = 123) +
  geom_text_repel(data = KA_reg[KA_reg$class == "tumour" | KA_reg$class == "perturbation" | (KA_reg$class == "both" & KA_reg$pP > 0.15), ], mapping = aes(x = pP, y = pT, label = kinase), box.padding = 0.3, point.padding = 0.3, size = 4, colour="black", seed = 123) +
  stat_cor(mapping = aes(x = pP, y = pT), method = "pearson", size = 4, label.y.npc = 1, label.x.npc = 0, label.sep = "\n") +
  theme_classic() +
  theme(
    plot.margin = unit(c(0.1,0.5,0.1,0.1), "cm"),
    axis.title = element_text(colour = "black", size = 15),
    axis.text = element_text(colour = "black", size = 13),
    legend.title = element_text(colour = "black", size = 16),
    legend.text = element_text(colour = "black", size = 16),
    legend.position = "bottom",
    legend.direction = "horizontal") +
  labs(x = "Kinase regulation (% of perturbations)", y = "Kinase regulation (% of tumours)") +
  guides(color = guide_legend(override.aes = list(alpha=1, size = 3), ncol = 1))

#+ fig.width=4, fig.height=5
KA_reg_plot

ggsave(filename = "kin_activity_number_perturbations_tumours_correlation_all_samples2.png", plot = KA_reg_plot, path = "./output/plots/generalists_specialists/", width = 4, height = 5)
ggsave(filename = "kin_activity_number_perturbations_tumours_correlation_all_samples2.pdf", plot = KA_reg_plot, path = "./output/plots/generalists_specialists/", width = 4, height = 5)


#' kinase regulation by cancer tissue type\
#' compared against cell lines perturbations
KA_reg_tissue <- gs_mapper(KA_tumors, KA_pertb, group = T, group_var = tissue, metadataCancer = metadata) %>%
  mutate(tissue = fct_recode(tissue, Brain = "brain", Breast = "breast", Colorectal = "colorectal", Kidney = "kidney", Liver = "liver", Lung = "lung", Ovary = "ovary", Stomach = "stomach", Uterus = "uterus"))

KA_reg_plot <- KA_reg_tissue %>%
  ggplot() +
  theme_classic() +
  facet_wrap(~ tissue, scales = "free") +
  geom_point(mapping = aes(x = pP, y = pT, color = class), alpha = 0.7, size = 0.2) +
  geom_smooth(mapping = aes(x = pP, y = pT), method = "lm", lwd = 0.1) +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Regulated kinase") +
  geom_text_repel(data = KA_reg_tissue[KA_reg_tissue$class != "both", ], mapping = aes(x = pP, y = pT, label = kinase), box.padding = 0.3, point.padding = 0.3, segment.size = 0.1, size = 2, colour="black") +
  #stat_cor(mapping = aes(x = pP, y = pT), method = "pearson", size = 3) +
  theme(
    strip.background = element_blank(),
    axis.text = element_text(colour = "black", size = 6),
    axis.ticks = element_line(colour = "black", size = 0.5)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  labs(x = "Perturbations (# conditions regulated)", y = "Tumour samples/CCLE (# samples regulated)")

#+ fig.width=8, fig.height=4
KA_reg_plot

ggsave(filename = "kin_activity_number_perturbations_tumours_correlation_by_tissue.png", plot = KA_reg_plot, path = "./output/plots/generalists_specialists/", width = 8, height = 4)
ggsave(filename = "kin_activity_number_perturbations_tumours_correlation_by_tissue.pdf", plot = KA_reg_plot, path = "./output/plots/generalists_specialists/", width = 8, height = 4)

helper <- function(vect, value){
  x <- vect
  #x[x == min(x)] <- min(x)-value
  x[x == max(x)] <- max(x)+value
  
  return(x)
}

dummy <- KA_reg_tissue %>%
  select(tissue, pP, pT) %>%
  group_by(tissue) %>%
  summarise_if(is.numeric, .funs = ~ list(range(.x))) %>%
  ungroup() %>%
  mutate(increment = c(0.05, 0.25, rep(0.4,4), 0.15, 0.14, 0.18)) %>%
  mutate(pT = map2(.x = pT, .y = increment, .f = helper)) %>%
  #mutate_if(is.list, .funs = function(x, y){map2(.x = x, .y = y, .f = helper)}, y = rep(0.1,9)) %>%
  #mutate_if(is.list, .funs = ~ map(.x = ..1, .f = helper, value = 0.1)) %>%
  unnest()

KA_reg_plot <- KA_reg_tissue %>%
  mutate(class = fct_relevel(class, "tumour", "both", "perturbation")) %>%
  ggplot() +
  geom_point(mapping = aes(x = pP, y = pT, color = class, alpha = class, size = class)) +
  geom_smooth(mapping = aes(x = pP, y = pT), method = "lm", se = T, size = 0.25) +
  #geom_blank(data = dummy, mapping = aes(x = pP, y = pT)) +
  facet_wrap(~ tissue, scales = "free") +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Regulated kinase") +
  scale_alpha_manual(values = c(1,0.5,1), guide = F) +
  scale_size_manual(values = c(2,1,2), guide = F) +
  #geom_text_repel(data = KA_reg_tissue[KA_reg_tissue$class != "both", ], mapping = aes(x = pP, y = pT, label = kinase), box.padding = 0.7, point.padding = 0.7, size = 5, segment.size = 0.2, colour="black", seed = 123) +
  geom_text_repel(data = KA_reg_tissue[KA_reg_tissue$class != "both", ], mapping = aes(x = pP, y = pT, label = kinase), box.padding = 0.2, point.padding = 0.7, size = 5, segment.size = 0.2, colour="black", seed = 123) +
  #stat_cor(mapping = aes(x = pP, y = pT), method = "pearson", size = 4.5, label.y.npc = 0.95, label.x.npc = 0.4) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(colour = "black", size = 26),
    legend.text = element_text(colour = "black", size = 24),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 20),
    axis.title = element_text(colour = "black", size = 22),
    axis.text = element_text(colour = "black", size = 18)) +
  guides(color = guide_legend(override.aes = list(alpha=1, size = 5), nrow = 1)) +
  labs(x = "Kinase regulation (% of perturbations)", y = "Kinase regulation (% of tumours)")


#+ fig.width=7, fig.height=5
KA_reg_plot

ggsave(filename = "kin_activity_number_perturbations_tumours_correlation_by_tissue2.png", plot = KA_reg_plot, path = "./output/plots/generalists_specialists/", width = 12, height = 10)
ggsave(filename = "kin_activity_number_perturbations_tumours_correlation_by_tissue2.pdf", plot = KA_reg_plot, path = "./output/plots/generalists_specialists/", width = 12, height = 10)


#' number of specific kinases by tissue
#+ fig.width=5, fig.height=5
KA_reg_tissue %>%
  filter(class != "both") %>%
  ggplot(mapping = aes(x = class, fill = class)) +
  geom_bar(show.legend = F) +
  facet_wrap(~ tissue) +
  scale_y_continuous(breaks = seq(0,10,2)) +
  scale_fill_manual(values = c("#377eb8", "#4daf4a")) +
  theme(axis.text = element_text(colour = "black")) +
  labs(y = "Number of specific kinases", x = "Class")

ggsave(filename = "NspecificKinases_by_tissue.pdf", path = "./output/plots/generalists_specialists/", width = 5, height = 5)
ggsave(filename = "NspecificKinases_by_tissue.png", path = "./output/plots/generalists_specialists/", width = 5, height = 5)

KA_reg_tissue %>%
  filter(class != "both") %>%
  mutate(tissue = fct_infreq(tissue)) %>%
  ggplot(mapping = aes(x = tissue, fill = class)) +
  geom_bar(position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = c("#4daf4a", "#e41a1c"), labels = c("perturbation" = "Perturbation", "tumour" = "Tumour")) +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  theme_classic() +
  theme(
    axis.title.x = element_text(colour = "black", size = 16),
    axis.title.y = element_blank(),
    axis.text = element_text(colour = "black", size = 14),
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", size = 16),
    legend.position = "bottom") +
  labs(y = "Regulated kinases")

ggsave(filename = "NspecificKinases_by_tissue2.pdf", path = "./output/plots/generalists_specialists/", width = 3, height = 6)
ggsave(filename = "NspecificKinases_by_tissue2.png", path = "./output/plots/generalists_specialists/", width = 3, height = 6)


KA_reg_tissue %>%
  filter(class != "both") %>%
  mutate(tissue = fct_infreq(tissue)) %>%
  ggplot(mapping = aes(x = tissue, fill = class)) +
  geom_bar(position = "dodge", show.legend = F) +
  coord_flip() +
  scale_fill_manual(values = c("#4daf4a", "#e41a1c"), labels = c("perturbation" = "Perturbation", "tumour" = "Tumour")) +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  theme_classic() +
  theme(
    axis.title.x = element_text(colour = "black", size = 16),
    axis.title.y = element_text(colour = "black", size = 16),
    axis.text = element_text(colour = "black", size = 14)) +
  labs(x = "Tissues", y = "Regulated kinases")

ggsave(filename = "NspecificKinases_by_tissue3.pdf", path = "./output/plots/generalists_specialists/", width = 3, height = 3)
ggsave(filename = "NspecificKinases_by_tissue3.png", path = "./output/plots/generalists_specialists/", width = 3, height = 3)


KA_reg_tissue %>%
  filter(class != "both") %>%
  mutate(tissue = fct_infreq(tissue)) %>%
  ggplot(mapping = aes(x = tissue, fill = class)) +
  geom_bar(position = "dodge", show.legend = F) +
  #coord_flip() +
  scale_fill_manual(values = c("#4daf4a", "#e41a1c"), labels = c("perturbation" = "Perturbation", "tumour" = "Tumour")) +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  theme_classic() +
  theme(
    axis.title.x = element_text(colour = "black", size = 16),
    axis.title.y = element_text(colour = "black", size = 14),
    axis.text.x = element_text(colour = "black", size = 14, angle = 90, hjust = 1),
    axis.text.y = element_text(colour = "black", size = 14)) +
  labs(x = "Tissues", y = "Number of kinases")

ggsave(filename = "NspecificKinases_by_tissue4.pdf", path = "./output/plots/generalists_specialists/", width = 3, height = 3)
ggsave(filename = "NspecificKinases_by_tissue4.png", path = "./output/plots/generalists_specialists/", width = 3, height = 3)


#' number of tissues by kinase
#+ fig.width=5, fig.height=5
KA_reg_tissue %>%
  filter(class != "both") %>%
  mutate(kinase = fct_infreq(kinase)) %>%
  ggplot(mapping = aes(x = kinase, fill = class)) +
  geom_bar(show.legend = F) +
  coord_flip() +
  facet_wrap(~ class) +
  scale_y_continuous(breaks = c(0:7)) +
  scale_fill_manual(values = c("#377eb8", "#4daf4a")) +
  theme(axis.text = element_text(colour = "black"), axis.text.y = element_text(size = 8)) +
  labs(y = "Number of tissues", x = "Kinase")

ggsave(filename = "specificKinases_numberOfTissues.pdf", path = "./output/plots/generalists_specialists/", width = 5, height = 5)
ggsave(filename = "specificKinases_numberOfTissues.png", path = "./output/plots/generalists_specialists/", width = 5, height = 5)

KA_reg_tissue %>%
  filter(class != "both") %>%
  mutate(class = fct_relevel(class, "tumour", "perturbation")) %>%
  mutate(kinase = fct_rev(fct_infreq(kinase))) %>%
  mutate(kinase = fct_rev(fct_relevel(kinase, "RPS6KA3", "RPS6KA1", "MAPK13", "MAPKAPK2", "RPS6KB1"))) %>%
  ggplot(mapping = aes(x = kinase, fill = class)) +
  geom_bar(show.legend = F, position = "dodge") +
  coord_flip() +
  #facet_wrap(~ class) +
  scale_y_continuous(breaks = c(0:7)) +
  scale_fill_manual(values = c("#e41a1c", "#4daf4a")) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.title = element_text(colour = "black", size = 14),
    axis.text.y = element_text(size = 9, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black")) +
  labs(x = "Kinase", y = "Tissues")

ggsave(filename = "specificKinases_numberOfTissues2.pdf", path = "./output/plots/generalists_specialists/", width = 2, height = 5)
ggsave(filename = "specificKinases_numberOfTissues2.png", path = "./output/plots/generalists_specialists/", width = 2, height = 5)


KA_reg_tissue %>%
  filter(class != "both") %>%
  mutate(class = fct_relevel(class, "tumour", "perturbation")) %>%
  mutate(kinase = fct_rev(fct_infreq(kinase))) %>%
  mutate(kinase = fct_rev(fct_relevel(kinase, "RPS6KA3", "RPS6KA1", "MAPK13", "MAPKAPK2", "RPS6KB1"))) %>%
  ggplot(mapping = aes(x = kinase, fill = class)) +
  geom_bar(show.legend = F, position = "dodge") +
  scale_y_continuous(breaks = c(0:7)) +
  scale_fill_manual(values = c("#e41a1c", "#4daf4a")) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.title.x = element_text(colour = "black", size = 16),
    axis.title.y = element_text(colour = "black", size = 14),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black", angle = 90, hjust = 1)) +
  labs(x = "Kinases", y = "Number of tissues")

ggsave(filename = "specificKinases_numberOfTissues3.pdf", path = "./output/plots/generalists_specialists/", width = 7, height = 3)
ggsave(filename = "specificKinases_numberOfTissues3.png", path = "./output/plots/generalists_specialists/", width = 7, height = 3)


#' number of specific kinases in the tumours by tissue
specific_kin <- KA_reg_tissue %>%
  filter(class == "tumour") %>%
  group_by(tissue) %>%
  summarise(specificKin = n()) %>%
  ungroup() %>%
  mutate(tissue = fct_relabel(tissue, .fun = tolower))

#' number of missing values by tissue
missing_values <- KA_tumors %>%
  inner_join(metadata, by = "sample") %>%
  group_by(tissue) %>%
  summarise(N_NAs = sum(is.na(log10P))) %>%
  ungroup()

#' number of measurements with more than 3 substrates by tissue
kinAct3Subs1 <- KA_tumors %>%
  filter(n >= 3) %>%
  inner_join(metadata, by = "sample") %>%
  group_by(tissue) %>%
  summarise(`NKM_>3Substrates` = n()) %>%
  ungroup()

#' number of kinases identified by tissue
kinAct3Subs2 <- KA_tumors %>%
  filter(n >= 3) %>%
  inner_join(metadata, by = "sample") %>%
  group_by(tissue) %>%
  summarise(Kinases = length(unique(kinase))) %>%
  ungroup()

#' number of samples by tissue
Nsamples <- metadata %>%
  group_by(tissue) %>%
  summarise(Samples = n()) %>%
  ungroup()


#' correlate the number of specific kinases in each tissue with the other measures
dat <- specific_kin %>%
  #inner_join(missing_values, by = "tissue") %>%
  #inner_join(kinAct3Subs1, by = "tissue") %>%
  inner_join(kinAct3Subs2, by = "tissue") %>%
  inner_join(Nsamples, by = "tissue") %>%
  pivot_longer(-c(tissue, specificKin), names_to = "data_type", values_to = "n")

correlations <- dat %>%
  group_by(data_type) %>%
  summarise(x = list(broom::tidy(cor.test(specificKin, n)))) %>%
  ungroup() %>%
  unnest() %>%
  select(data_type, estimate, p.value) %>%
  mutate(string = str_c("r = ", round(estimate, 2), ", P = ", round(p.value, 3)))

scatterplot <- dat %>%
  ggplot(mapping = aes(x = specificKin, y = n, label = tissue)) +
  geom_point(size = 4, alpha = 0.5, shape=21, stroke=1, fill = "black") +
  geom_text_repel(seed = 123, size = 6, color = "#e41a1c") +
  facet_wrap(~ data_type, scales = "fixed", nrow = 2,
             labeller = labeller(data_type = c("Kinases" = str_c("Kinases\n", correlations$string[1]), "Samples" = str_c("Samples\n", correlations$string[2])))) +
  #stat_cor(method = "pearson", label.x.npc = 0.4, label.y.npc = 0.3, size = 4, color = "blue") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 18, colour = "black"),
    axis.title = element_text(colour = "black", size = 20),
    axis.text = element_text(size = 18, colour = "black")) +
  labs(x = "Regulated kinases in tumours", y = "Number of cancer samples and quantified kinases")

#+ fig.width=8, fig.height=4
scatterplot

ggsave(filename = "specificKinases_N_NAs_samples_correlation.pdf", plot = scatterplot, path = "./output/plots/generalists_specialists/", width = 4, height = 8)
ggsave(filename = "specificKinases_N_NAs_samples_correlation.png", plot = scatterplot, path = "./output/plots/generalists_specialists/", width = 4, height = 8)





#' kinase regulation by cancer tissue type\
#' compared against cell lines perturbations\
#' select the perturbations randomly 100 times
KA_reg <- gs_mapper(KA_tumors, KA_pertb, group = T, group_var = tissue, metadataCancer = metadata, perm = T, perm_n = 100)

KA_reg_plot <- KA_reg %>%
  ggplot() +
  theme_classic() +
  facet_wrap(~ tissue, scales = "free") +
  geom_point(mapping = aes(x = nP, y = nT, color = class), alpha = 0.7, size = 0.2) +
  geom_smooth(mapping = aes(x = nP, y = nT), method = "lm", lwd = 0.1) +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Regulated kinase") +
  geom_text_repel(data = KA_reg[KA_reg$class != "both", ], mapping = aes(x = nP, y = nT, label = kinase), box.padding = 0.3, point.padding = 0.3, segment.size = 0.1, size = 2, colour="black") +
  #stat_cor(mapping = aes(x = nP, y = nT), method = "pearson", size = 3) +
  theme(
    strip.background = element_blank(),
    axis.text = element_text(colour = "black", size = 6),
    axis.ticks = element_line(colour = "black", size = 0.5)) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  labs(x = "Perturbations (# conditions regulated)", y = "Tumour samples/CCLE (# samples regulated)")

#+ fig.width=8, fig.height=4
KA_reg_plot

ggsave(filename = "kin_activity_number_perturbations_tumours_correlation_by_tissue_perm100.png", plot = KA_reg_plot, path = "./output/plots/generalists_specialists/", width = 8, height = 4)
ggsave(filename = "kin_activity_number_perturbations_tumours_correlation_by_tissue_perm100.pdf", plot = KA_reg_plot, path = "./output/plots/generalists_specialists/", width = 8, height = 4)
unlink("kin_activity_number_perturbations_tumours_correlation_by_tissue_perm100.png")
unlink("kin_activity_number_perturbations_tumours_correlation_by_tissue_perm100.pdf")


#' number of specific kinases by tissue
#+ fig.width=5, fig.height=5
KA_reg %>%
  filter(class != "both") %>%
  ggplot(mapping = aes(x = class, fill = class)) +
  geom_bar(show.legend = F) +
  facet_wrap(~ tissue) +
  scale_y_continuous(breaks = seq(0,10,2)) +
  scale_fill_manual(values = c("#377eb8", "#4daf4a")) +
  theme(axis.text = element_text(colour = "black")) +
  labs(y = "Number of specific kinases", x = "Class")

ggsave(filename = "NspecificKinases_by_tissue_perm100.pdf", path = "./output/plots/generalists_specialists/", width = 5, height = 5)
ggsave(filename = "NspecificKinases_by_tissue_perm100.png", path = "./output/plots/generalists_specialists/", width = 5, height = 5)
unlink("NspecificKinases_by_tissue_perm100.pdf")
unlink("NspecificKinases_by_tissue_perm100.png")


#' number of tissues by kinase
#+ fig.width=5, fig.height=5
KA_reg %>%
  filter(class != "both") %>%
  mutate(kinase = fct_infreq(kinase)) %>%
  ggplot(mapping = aes(x = kinase, fill = class)) +
  geom_bar(show.legend = F) +
  coord_flip() +
  facet_wrap(~ class) +
  scale_y_continuous(breaks = c(0:7)) +
  scale_fill_manual(values = c("#377eb8", "#4daf4a")) +
  theme(axis.text = element_text(colour = "black"), axis.text.y = element_text(size = 8)) +
  labs(y = "Number of tissues", x = "Kinase")

ggsave(filename = "specificKinases_numberOfTissues_perm100.pdf", path = "./output/plots/generalists_specialists/", width = 5, height = 5)
ggsave(filename = "specificKinases_numberOfTissues_perm100.png", path = "./output/plots/generalists_specialists/", width = 5, height = 5)
unlink("specificKinases_numberOfTissues_perm100.pdf")
unlink("specificKinases_numberOfTissues_perm100.png")

