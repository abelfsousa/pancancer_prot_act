#' ---
#' title: "Distribution of the kinase-activity scores using CPTAC phosphorylation data"
#' author: "abelsousa"
#' ---

#+ setup, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)



#' load R packages
library(RColorBrewer)
library(tidyverse)
library(viridis)


#' import CPTAC kinase-activity inference data\
files <- c("phosphoQ_Protregout_allSamp_withZtransf.txt.gz", "phosphoQ_RNAregout_allSamp_withZtransf.txt.gz", "phosphoQ.txt.gz")
ka <- map(.x = files, .f = ~ read_tsv(paste0("./output/files/CPTAC_KA_", .x)))
ka <- map2(.x = ka, .y = str_split(files, "\\.txt", simplify = T)[,1], .f = ~ .x %>% mutate(data = .y))
ka <- bind_rows(ka)



#' plot the distribution of the kinase-activity scores
plot <- ka %>%
  mutate(data = fct_recode(data, `phosphoQ protein reg-out` = "phosphoQ_Protregout_allSamp_withZtransf", `phosphoQ rna reg-out` = "phosphoQ_RNAregout_allSamp_withZtransf")) %>%
  filter(data == "phosphoQ protein reg-out") %>%
  ggplot(mapping = aes(x = source_type, y = log10P, fill = source_type)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ data, scales = "fixed") +
  coord_flip() +
  theme_classic() +
  theme(
    axis.title = element_text(colour = "black", size = 14),
    axis.text = element_text(colour = "black", size = 12),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 14)) +
  scale_fill_discrete(guide = F) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  labs(x = "Kinase-substrate list", y = "Kinase-activity score")

#+ fig.width=8, fig.height=3
plot

ggsave(filename="KA_cptac_distribution.png", plot = plot, path = "./output/plots/kinase_activity_inference/", width=8, height=3)
unlink("KA_cptac_distribution.png")
