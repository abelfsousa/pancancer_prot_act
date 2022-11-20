#' ---
#' title: "Correlate kinase activities between perturbation sets"
#' author: "Abel Sousa"
#' ---


#' load R packages
library(tidyverse)
library(ggpubr)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load kinase activities from perturbations

# David's first compilation
KA_pertb1 <- read_tsv(file = "./output/files/BK_phospho_KA.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  select(-source_type, -n) %>%
  rename(act1 = log10P)

# David/Danish's second compilation
KA_pertb2 <- read_tsv(file = "./output/files/KA_esetNR.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  select(-source_type, -n) %>%
  rename(act2 = log10P)


#' correlate kinase activities profiles for the same samples
pertb <- KA_pertb1 %>%
  inner_join(KA_pertb2, by = c("sample", "kinase")) %>%
  group_by(sample) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(na.exclude(.x)) >= 5)) %>%
  pull(data) %>%
  map_dfr(~ broom::tidy(cor.test(.x$act1, .x$act2))) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = factor(0), y = estimate)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "Perturbations", y = "Pearson correlation")

#+ fig.width=4, fig.height=4
pertb


pertb <- KA_pertb1 %>%
  inner_join(KA_pertb2, by = c("sample", "kinase")) %>%
  group_by(sample) %>%
  filter(sum(!(is.na(act1) | is.na(act2))) > 5) %>%
  summarise(act1 = median(act1, na.rm = T), act2 = median(act2, na.rm = T)) %>%
  ungroup() %>%
  ggplot(mapping = aes(x = act1, y = act2)) +
  geom_point() +
  stat_cor() +
  labs(x = "Fisrt compilation (sample kinase activity median)", y = "Second compilation (sample kinase activity median)")

#+ fig.width=4, fig.height=4
pertb


#' correlate kinase activities across samples
pertb <- KA_pertb1 %>%
  inner_join(KA_pertb2, by = c("sample", "kinase")) %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(na.exclude(.x)) >= 5)) %>%
  pull(data) %>%
  map_dfr(~ broom::tidy(cor.test(.x$act1, .x$act2))) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = factor(0), y = estimate)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x = "Kinases", y = "Pearson correlation")

#+ fig.width=4, fig.height=4
pertb


pertb <- KA_pertb1 %>%
  inner_join(KA_pertb2, by = c("sample", "kinase")) %>%
  group_by(kinase) %>%
  filter(sum(!(is.na(act1) | is.na(act2))) > 5) %>%
  summarise(act1 = median(act1, na.rm = T), act2 = median(act2, na.rm = T)) %>%
  ungroup() %>%
  filter(abs(act1) < 2 & abs(act2) < 2) %>%
  ggplot(mapping = aes(x = act1, y = act2)) +
  geom_point() +
  stat_cor() +
  labs(x = "Fisrt compilation (kinase activity median across samples)", y = "Second compilation (kinase activity median across samples)")

#+ fig.width=4, fig.height=4
pertb

