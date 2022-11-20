#' ---
#' title: "Kinase Activities Quantiles"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=6, fig.height=6)


#' load R packages
library(tidyverse)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load kinase activities from perturbations

## David/Danish's second compilation
KA_pertb <- read_tsv(file = "./output/files/KA_esetNR.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  select(-source_type) %>%
  filter(n>=3) %>%
  select(-n)


#' load kinase activities form tumor samples
KA_tumors <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  select(-source_type) %>%
  filter(n>=3) %>%
  select(-n)


#' merge perturbation and tumour kinase activities\
#' calculate kinase activities quantiles
KA_Quantiles <- KA_pertb %>%
  mutate(class = "perturbations") %>%
  bind_rows(mutate(KA_tumors, class = "tumours")) %>%
  select(class, everything()) %>%
  group_by(class, sample) %>%
  summarise(q = list(broom::tidy(quantile(log10P)))) %>%
  ungroup() %>%
  unnest() %>%
  group_by(class, names) %>%
  summarise(mean = mean(x), sd = sd(x)) %>%
  ungroup() %>%
  mutate(names = fct_relevel(names, "0%", "25%", "50%", "75%", "100%"))


#' plot kinase activities quantiles
KA_Quantiles_plot <- KA_Quantiles %>%
  filter(!names %in% c("0%", "100%")) %>%
  ggplot(mapping = aes(x = names, y = mean, fill = names)) +
  geom_col() +
  geom_text(mapping = aes(x = names, y = sign(mean)*(abs(mean)+0.05), label = round(mean, 3)), size = 2) +
  geom_errorbar(mapping = aes(ymin = mean-(sd*2), ymax = mean+(sd*2), color = names), size = 0.2, width = 0.3, position = "dodge") +
  facet_wrap(~ class) +
  scale_fill_discrete(guide = F) +
  scale_color_discrete(guide = F) +
  theme(
    strip.text = element_text(size = 8, colour = "black"),
    axis.ticks = element_line(size = 0.1, colour = "black"),
    axis.text = element_text(size = 6, colour = "black"),
    axis.title = element_text(size = 8, colour = "black")) +
  labs(x = "Quantiles", y = "Mean")


#+ fig.width=3, fig.height=4
KA_Quantiles_plot


ggsave(filename = "kin_activities_quantiles_tumor_perturbations.png", plot = KA_Quantiles_plot, path = "./output/plots/generalists_specialists/", width = 3, height = 4)
unlink("kin_activities_quantiles_tumor_perturbations.png")
