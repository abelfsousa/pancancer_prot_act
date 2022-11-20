#' ---
#' title: "Correlation of kinase activity measured by substrates and kinase regulatory phosphosites"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=6, fig.height=6)


#' load R packages
library(tidyverse)
library(ggpubr)

set.seed(123)

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


kin_subN <- 3
kin_siteN <- 3

ka_kinSites <- read_tsv(file = "./output/files/CPTAC_KA_kin_psites_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  rename(kinase=gene, kinSiteActv=log10P) %>%
  filter(n >= kin_siteN) %>%
  select(-n)

ka_kinSub <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  select(-source_type) %>%
  rename(kinSubActv=log10P) %>%
  filter(n >= kin_subN) %>%
  select(-n)

ka_cor <- inner_join(ka_kinSites, ka_kinSub, by = c("sample", "kinase")) %>% 
  ggplot(mapping = aes(x = kinSiteActv, y = kinSubActv)) +
  geom_point() +
  geom_smooth(method = "lm") +
  stat_cor() + 
  theme_classic()

#+ fig.width=4, fig.height=4
ka_cor

ggsave(filename = "ka_kinPsites_kinSub_cor_scatter.png", plot = ka_cor, path = "./output/plots/kinase_activity_substrate_psite_correlation/", height = 4, width = 4)
unlink("ka_kinPsites_kinSub_cor_scatter.png")


ka_corByKin <- inner_join(ka_kinSites, ka_kinSub, by = c("sample", "kinase")) %>%
  group_by(kinase) %>%
  mutate(n = n()) %>%
  filter(n >= 5) %>%
  select(-n) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x=data, .f=~broom::tidy(cor.test(.x$kinSiteActv, .x$kinSubActv, method = "pearson")))) %>%
  select(-data) %>%
  unnest() %>%
  select(kinase, cor = estimate) %>%
  mutate(class = "True kinase")

shuffle_tibble <- function(tib){
  tib_shuf <- tib %>%
    mutate(index = sample(seq_len(nrow(tib)), nrow(tib), replace = FALSE)) %>%
    mutate(sample = sample[index]) %>%
    select(-index)
}

ka_kinSub_shuff <- ka_kinSub %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(data2 = map(.x = data, .f = shuffle_tibble)) %>%
  select(-data) %>%
  unnest()

ka_kinSites_shuff <- ka_kinSites %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(data2 = map(.x = data, .f = shuffle_tibble)) %>%
  select(-data) %>%
  unnest()

ka_corByKin_shuff <- inner_join(ka_kinSites_shuff, ka_kinSub_shuff, by = c("sample", "kinase")) %>%
  group_by(kinase) %>%
  mutate(n = n()) %>%
  filter(n >= 5) %>%
  select(-n) %>%
  nest() %>%
  ungroup() %>%
  mutate(cor = map(.x=data, .f=~broom::tidy(cor.test(.x$kinSiteActv, .x$kinSubActv, method = "pearson")))) %>%
  select(-data) %>%
  unnest() %>%
  select(kinase, cor = estimate) %>%
  mutate(class = "Shuffled kinase")


ka_cor_dist <- bind_rows(ka_corByKin, ka_corByKin_shuff) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = class, y = cor), notch = TRUE, outlier.shape = NA) +
  geom_jitter(mapping = aes(x = class, y = cor, color = class), alpha = 0.7, width = 0.1) +
  #geom_dotplot(mapping = aes(x = class, y = cor, color = class), binaxis="y", stackdir="center", dotsize = 0.5) +
  stat_compare_means(mapping = aes(x = class, y = cor), label.y = -0.75, label.x = 2.5) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12)
  ) +
  scale_color_discrete(guide=F) +
  scale_y_continuous(limits = c(-1,1)) +
  labs(x = "", y = "Pearson's r (substrate vs reg. phosphosite kinase activity)")

#+ fig.width=6, fig.height=2
ka_cor_dist

ggsave(filename = "ka_kinPsites_kinSub_cor_boxplot.png", plot = ka_cor_dist, path = "./output/plots/kinase_activity_substrate_psite_correlation/", height = 2, width = 6)
unlink("ka_kinPsites_kinSub_cor_boxplot.png")

