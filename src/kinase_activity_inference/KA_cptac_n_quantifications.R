#' ---
#' title: "Description of CPTAC kinase activity inference data - number of quantifications"
#' author: "abelsousa"
#' ---

#+ setup, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)



#' load R packages
library(tidyverse)

source("./src/utils/getSamples.R")


#' load all samples gathered in this study
all_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")


#' get samples with phosphorylation, protein and mRNA data
samples <- getSamples(all_samples, c("phosphorylation", "mRNA", "protein")) %>%
  select(sample)


#' import CPTAC kinase-activity inference data\
files <- c("phosphoQ_Protregout_allSamp_withZtransf", "phosphoQ_RNAregout_allSamp_withZtransf", "phosphoQ")
ka <- map(.x = files, .f = ~ read_tsv(str_c("./output/files/CPTAC_KA_", .x, ".txt.gz", sep = "")))
ka <- map2(.x = ka, .y = str_split(files, "_all", simplify = T)[,1], .f = ~ .x %>% mutate(data = .y))
ka <- bind_rows(ka)


#' get number of kinases infered with more than 3, 6, 10 and 15 substrates in each sample\
#' use database/tex-mining kinase-substrate pairs
kaf <- map(.x=c(3,6,10,15), .f = ~ ka[ka$source_type == "DB_text-mining", ] %>% group_by(data, sample) %>% summarise(n = sum(n > .x)) %>% ungroup())
kaf <- map2(.x = kaf, .y=c(3,6,10,15), .f = ~ .x %>% mutate(filter = .y))
kaf <- bind_rows(kaf)


#' # make boxplots from the number of quantifications
boxpl_ka <- kaf %>%
  filter(sample %in% samples$sample) %>%
  mutate(data = fct_reorder(.f=data, .x=n, .fun = median)) %>%
  mutate(data = fct_recode(.f=data, prot_reg="phosphoQ_Protregout", rna_reg="phosphoQ_RNAregout", not_reg="phosphoQ")) %>%
  mutate(filter = fct_recode(.f=as.character(filter), `>3`="3", `>6`="6", `>10`="10", `>15`="15")) %>%
  mutate(filter = fct_inorder(filter)) %>%
  ggplot(mapping = aes(x = data, y = n, fill = data)) +
  geom_boxplot(notch = TRUE, outlier.size = 0.1, show.legend = FALSE) +
  facet_wrap(~ filter, scales = "fixed") +
  coord_flip() +
  theme_classic() +
  theme(
    plot.title = element_text(colour = "black", size = 14),
    axis.title = element_text(colour = "black", size = 14),
    axis.text = element_text(colour = "black", size = 12),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 14)) +
  labs(x = "Phosphorylation data", y = "Number of kinases with more than # substrates in each sample", title = "DB/text-mining lists")

#+ fig.width=8, fig.height=3
boxpl_ka

ggsave(filename="KA_cptac_kinase_Nsub_samples_boxpl.png", plot = boxpl_ka, path = "./output/plots/kinase_activity_inference/", width=8, height=3)
unlink("KA_cptac_kinase_Nsub_samples_boxpl.png")



#' # make barplot from the total number of quantifications
barpl_ka <- kaf %>%
  filter(sample %in% samples$sample) %>%
  group_by(filter, data) %>%
  summarise(n = sum(n)) %>%
  ungroup() %>%
  mutate(data = fct_reorder(.f=data, .x=n, .fun = median)) %>%
  mutate(data = fct_recode(.f=data, prot_reg="phosphoQ_Protregout", rna_reg="phosphoQ_RNAregout", not_reg="phosphoQ")) %>%
  mutate(filter = fct_recode(.f=as.character(filter), `>3`="3", `>6`="6", `>10`="10", `>15`="15")) %>%
  mutate(filter = fct_inorder(filter)) %>%
  ggplot(mapping = aes(x = data, y = n, fill = data)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ filter, scales = "fixed") +
  coord_flip() +
  theme_classic() +
  theme(
    plot.title = element_text(colour = "black", size = 14),
    axis.title = element_text(colour = "black", size = 14),
    axis.text = element_text(colour = "black", size = 12),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 14),
    panel.spacing = unit(1, "lines")) +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  labs(x = "Phosphorylation data", y = "Total number of quantifications with more than # substrates", title = "DB/text-mining lists")

#+ fig.width=8, fig.height=3
barpl_ka

ggsave(filename="KA_cptac_kinase_Nsub_total.png", plot = barpl_ka, path = "./output/plots/kinase_activity_inference/", width=8, height=3)
unlink("KA_cptac_kinase_Nsub_total.png")
