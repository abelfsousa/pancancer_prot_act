#' ---
#' title: "CNV buffering"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=6, fig.height=6)


#' load R packages and external functions
library(mclust)
library(ggpubr)
library(viridis)
library(tidyverse)
library(cowplot)
library(RColorBrewer)

set.seed(123)

source(file = "src/utils/getSamples.R")

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load CNV data
cnv <- read_tsv(file = "./output/files/cnv.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "cnv")


#' load protein data
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "prot_log2fc")


#' load RNA data
rna <- read_tsv(file = "./output/files/transcriptomics_log2fc.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "rna_log2fc")


#' load metadata
all_samples <- read_tsv("./output/files/all_samples_annotation.txt")

samples_metadata <- getSamples(all_samples, data_types = c("cnv", "mRNA", "protein", "clinical"))

multi_data <- cnv %>%
  inner_join(protein, by = c("gene", "sample")) %>%
  inner_join(rna, by = c("gene", "sample"))
  
samples <- length(unique(multi_data$sample))

#' set a correlation function
corr <- function(df, x, y, method){
  a = df[[x]]
  b = df[[y]]
  cor <- broom::tidy(cor.test(a, b, method = method))
  
  res <- cor %>%
    select(estimate, p.value)
  
  res
}


#' CNV buffering analysis
cnv_buff <- multi_data %>%
  filter(!(is.na(prot_log2fc) | is.na(rna_log2fc))) %>%
  group_by(gene) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > samples*0.1)) %>%
  mutate(cnv_prot = map(.x = data, .f = corr, x = "prot_log2fc", y = "cnv", "pearson")) %>%
  mutate(cnv_rna = map(.x = data, .f = corr, x = "rna_log2fc", y = "cnv", "pearson")) %>%
  select(-data) %>%
  unnest() %>%
  rename(corr_cnv_prot = estimate, pval_cnv_prot = p.value) %>%
  rename(corr_cnv_rna = estimate1, pval_cnv_rna = p.value1) %>%
  mutate(attenuation = corr_cnv_rna - corr_cnv_prot) %>%
  mutate(cluster = as.character(mclust::Mclust(attenuation, G=3)$classification)) %>%
  arrange(desc(attenuation)) %>%
  mutate(class = if_else(cluster==1, "rna-attenuated", if_else(cluster==2, "non-attenuated", "prot-attenuated"))) %>%
  mutate(class = fct_relevel(class, "rna-attenuated", "non-attenuated", "prot-attenuated"))


#' Attenuation plot
correlation_plot <- ggplot(cnv_buff) +
  geom_point(aes(x=corr_cnv_rna, y=corr_cnv_prot, color=class), size=1, alpha = 0.5) +
  geom_density2d(mapping=aes(x=corr_cnv_rna, y=corr_cnv_prot, group=class), color="grey", size = 0.2) +
  geom_abline(slope=1, intercept=0, linetype=2, color="black") +
  geom_vline(xintercept=0, color="grey", size = 0.2, alpha = 0.5) +
  geom_hline(yintercept=0, color="grey", size = 0.2, alpha = 0.5) +
  theme_classic() +
  coord_fixed() +
  theme(axis.title=element_text(colour="black", size=15),
        axis.text=element_text(colour="black", size=12),
        plot.title = element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.position="bottom") +
  scale_x_continuous(limits = c(-0.2, 0.8), breaks=seq(-0.2, 0.8, by = 0.2)) +
  scale_y_continuous(limits = c(-0.2, 0.8), breaks=seq(-0.2, 0.8, by = 0.2)) +
  scale_color_brewer(type = "seq", palette = "Blues", labels = c("rna-attenuated", "non-attenuated", "prot-attenuated"), name = "Attenuation") +
  guides(color = guide_legend(title.position = "top", override.aes = list(size=3, shape=19))) +
  labs(x = "corr(CNV, RNA)", y = "corr(CNV, Protein)")

# Marginal densities along x axis
xdens <- axis_canvas(correlation_plot, axis = "x") +
  geom_density(data = cnv_buff, aes(x = corr_cnv_rna, fill = class), alpha = 0.95, size = 0.6) +
  ggpubr::fill_palette(palette = brewer.pal(3, "Blues"))

# Marginal densities along y axis
# Need to set coord_flip = TRUE, if you plan to use coord_flip()
ydens <- axis_canvas(correlation_plot, coord_flip = TRUE) +
  geom_density(data = cnv_buff, aes(x = corr_cnv_prot, fill = class), alpha = 0.95, size = 0.6) +
  coord_flip() +
  ggpubr::fill_palette(palette = brewer.pal(3, "Blues"))

p1 <- insert_xaxis_grob(correlation_plot, xdens, grid::unit(0.2, "null"), position = "top")
p2 <- insert_yaxis_grob(p1, ydens, grid::unit(0.2, "null"), position = "right")

#+ fig.width=5, fig.height=5
plot(p2)

ggsave(filename="cnv_buffering_scatterplot.png", plot = p2, path = "./output/plots/cnv_buffering/", width=5, height=5)
unlink("cnv_buffering_scatterplot.png")


#' load cnv buffering data from MCP paper
mcp <- read_tsv("./data/mcp_paper/protein_attenuation.txt") %>%
  mutate(class = str_replace(class, "-protein", ""))

buff <- cnv_buff %>%
  select(gene, attenuation1=attenuation, class1=class) %>%
  inner_join(mcp[, c("gene", "attenuation", "class")], by = c("gene")) %>%
  rename(attenuation2 = attenuation, class2 = class) %>%
  mutate(same = if_else(class1 == class2, "same", "different"))


#' Attenuation correlation between MCP paper and current data
corr_attenuation <- ggplot(data = buff, mapping = aes(x = attenuation1, y = attenuation2)) +
  geom_point() +
  stat_cor() +
  labs(x = "Current data", y = "MCP paper data")

#+ fig.width=5, fig.height=5
corr_attenuation

ggsave(filename="cnv_buffering_corr_mcp.png", plot = corr_attenuation, path = "./output/plots/cnv_buffering/", width=5, height=5)
unlink("cnv_buffering_corr_mcp.png")


#' Attenuation plot (genes changing attenuation category are highlighted)
# cnv_buff2 <- cnv_buff %>%
#   inner_join(mcp[, c("gene", "attenuation", "class")], by = c("gene")) %>%
#   rename(attenuationNow = attenuation.x, classNow = class.x, attenuationMCP = attenuation.y, classMCP = class.y) %>%
#   mutate(same = if_else(classNow == classMCP, "same", "different"))
# 
# correlation_plot2 <- ggplot(cnv_buff2) +
#   geom_point(aes(x=corr_cnv_rna, y=corr_cnv_prot, fill=same, color=classNow), size=1, shape = 21) +
#   geom_abline(slope=1, intercept=0, linetype=2, color="black") +
#   geom_vline(xintercept=0, color="grey", size = 0.2, alpha = 0.5) +
#   geom_hline(yintercept=0, color="grey", size = 0.2, alpha = 0.5) +
#   theme_classic() +
#   coord_fixed() +
#   theme(axis.title=element_text(colour="black", size=15),
#         axis.text=element_text(colour="black", size=12),
#         plot.title = element_blank(),
#         legend.text=element_text(size=12),
#         legend.title=element_text(size=14),
#         legend.position="bottom",
#         legend.box = "horizontal") +
#   scale_x_continuous(limits = c(-0.2, 0.8), breaks=seq(-0.2, 0.8, by = 0.2)) +
#   scale_y_continuous(limits = c(-0.2, 0.8), breaks=seq(-0.2, 0.8, by = 0.2)) +
#   scale_color_brewer(type = "seq", palette = "Blues", labels = c("non-attenuated", "lowly-attenuated", "highly-attenuated"), name = "Attenuation", guide=F) +
#   scale_fill_manual(values = c("#b2182b", "#ffffff")) +
#   guides(fill = guide_legend(title.position = "left", override.aes = list(size=2))) +
#   labs(x = "corr(CNV, RNA)", y = "corr(CNV, Protein)")
# 
# #+ fig.width=5, fig.height=5
# correlation_plot2
# 
# ggsave(filename="cnv_buffering_scatterplot2.png", plot = correlation_plot2, path = "./output/plots/cnv_buffering/", width=5, height=5)
# unlink("cnv_buffering_corr_mcp2.png")
