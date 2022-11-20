#' ---
#' title: "Correlation of kinase activity measured by substrates and kinase regulatory phosphosites"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=6, fig.height=6)


#' load R packages
library(tidyverse)
library(viridis)
library(ggpubr)
library(Hmisc)

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


#' correlate phosphosite-based and substrate-based activities of the same kinase
ka_corByKin <- inner_join(ka_kinSites, ka_kinSub, by = c("sample", "kinase")) %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(cor = map(.x=data, .f=~broom::tidy(cor.test(.x$kinSiteActv, .x$kinSubActv, method = "pearson")))) %>%
  select(-data) %>%
  unnest() %>%
  select(kinase, pearson_r = estimate, p_value = p.value) %>%
  mutate(class = "Same kinase") %>%
  rename(kinaseA = kinase) %>%
  mutate(kinaseB = kinaseA) %>%
  select(kinaseA, kinaseB, everything())


#' correlate phosphosite-based and substrate-based activities for all kinase pairs
ka_corByKin_all <- ka_kinSites %>%
  rename(kinaseA = kinase) %>%
  inner_join(ka_kinSub %>% rename(kinaseB = kinase), by = c("sample")) %>%
  group_by(kinaseA, kinaseB) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(.x = data, .f = ~ nrow(.x) > 10)) %>%
  mutate(cor = map(.x=data, .f=~broom::tidy(cor.test(.x$kinSiteActv, .x$kinSubActv, method = "pearson")))) %>%
  select(-data) %>%
  unnest() %>%
  select(kinaseA, kinaseB, pearson_r = estimate, p_value = p.value) %>%
  mutate(class = "All pairs") %>%
  filter(kinaseA != kinaseB)

# set up a function to convert a squared matrix into a data frame (tibble)
matrix_to_tb <- function(x, y, kA, kB){
  
  mat <- x
  mat_data <- y
  KIN_A <- kA
  KIN_B <- kB
  
  mat <- mat[kA,kB]
  df <- mat %>%
    as.data.frame() %>%
    rownames_to_column(var="kinaseA") %>%
    as_tibble() %>%
    pivot_longer(-c("kinaseA"), names_to="kinaseB", values_to=mat_data) %>%
    mutate(kinaseA = str_replace(kinaseA, "A_", ""), kinaseB = str_replace(kinaseB, "B_", ""))
  
  df
}

# set up a function to join the rcorr output into a single tibble
joinMatrices <- function(x, kinA, kinB){
  
  rcorr_matrices <- x
  rcorr_names <- names(rcorr_matrices)
  
  KIN_A <- kinA
  KIN_B <- kinB
  
  dfs <- map2(.x = rcorr_matrices, .y = rcorr_names, .f = matrix_to_tb, kA = KIN_A, kB = KIN_B)
  
  tb <- reduce(.x = dfs, .f = ~ inner_join(.x, .y, by = c("kinaseA", "kinaseB")))
  
  tb
}

ka_kinSites_mat <- ka_kinSites %>%
  filter(sample %in% intersect(unique(ka_kinSites$sample), unique(ka_kinSub$sample))) %>%
  arrange(sample) %>%
  mutate(kinase = str_c("A", kinase, sep = "_")) %>%
  pivot_wider(names_from = "sample", values_from = "kinSiteActv") %>%
  as.data.frame() %>%
  column_to_rownames(var = "kinase") %>%
  as.matrix() %>%
  t()

ka_kinSub_mat <- ka_kinSub %>%
  filter(sample %in% intersect(unique(ka_kinSites$sample), unique(ka_kinSub$sample))) %>%
  arrange(sample) %>%
  mutate(kinase = str_c("B", kinase, sep = "_")) %>%
  pivot_wider(names_from = "sample", values_from = "kinSubActv") %>%
  as.data.frame() %>%
  column_to_rownames(var = "kinase") %>%
  as.matrix() %>%
  t()

kinasesA <- colnames(ka_kinSites_mat)
kinasesB <- colnames(ka_kinSub_mat)

ka_corByKin_all2 <- rcorr(ka_kinSites_mat, ka_kinSub_mat) %>%
  joinMatrices(kinasesA, kinasesB) %>%
  filter(n>10) %>%
  select(-n, pearson_r=r, p_value=P) %>%
  mutate(class = "All pairs")

inner_join(ka_corByKin_all, ka_corByKin_all2, by = c("kinaseA", "kinaseB", "class")) %>%
  summarise(x = cor(pearson_r.x, pearson_r.y))


#' plot the distributions
ka_cor_dist <- bind_rows(ka_corByKin, ka_corByKin_all) %>%
  mutate(class = fct_recode(class, Others = "All pairs")) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = class, y = pearson_r, fill = class), notch = TRUE, color = "grey60", alpha = 0.8) +
  #geom_jitter(mapping = aes(x = class, y = pearson_r, color = class), alpha = 0.2, width = 0.1) +
  stat_compare_means(method = "t.test", mapping = aes(x = class, y = pearson_r), label.x = 2.4, label.y = 0.6, size = 3.5) +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()) +
  scale_fill_manual(values = c("#bdbdbd", "#67a9cf"), guide = F) +
  scale_y_continuous(limits = c(-1,1)) +
  labs(x = "", y = "Pearson's r (substrate vs phosphosite-based kinase activity)")

#+ fig.width=6, fig.height=2
ka_cor_dist

ggsave(filename = "kinase_activities_substrate_psite_based_correlation_box.png", plot = ka_cor_dist, path = "./output/plots/kinase_activity_substrate_psite_correlation/", height = 2, width = 6)
ggsave(filename = "kinase_activities_substrate_psite_based_correlation_box.pdf", plot = ka_cor_dist, path = "./output/plots/kinase_activity_substrate_psite_correlation/", height = 2, width = 6)
unlink("kinase_activities_substrate_psite_based_correlation_box.png")
unlink("kinase_activities_substrate_psite_based_correlation_box.pdf")

ka_cor_dist <- bind_rows(ka_corByKin, ka_corByKin_all) %>%
  mutate(class = fct_recode(class, Others = "All pairs")) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = class, y = pearson_r, fill = class), notch = TRUE, color = "grey60", alpha = 0.8) +
  #geom_jitter(mapping = aes(x = class, y = pearson_r, color = class), alpha = 0.2, width = 0.1) +
  stat_compare_means(method = "t.test", mapping = aes(x = class, y = pearson_r), label.x = 1, label.y = 1.1, size = 4) +
  #coord_flip() +
  theme_classic() +
  theme(
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(color = "black", size = 12)) +
  #scale_fill_manual(values = c("#bdbdbd", "#67a9cf"), guide = F) +
  scale_fill_viridis(discrete = T, guide = F, direction = -1) +
  scale_y_continuous(limits = c(-1,1.1)) +
  labs(x = "Kinase-kinase pair", y = "Pearson's r (substrate vs phosphosite-based kinase activity)")

#+ fig.width=3, fig.height=6
ka_cor_dist

ggsave(filename = "kinase_activities_substrate_psite_based_correlation_box2.png", plot = ka_cor_dist, path = "./output/plots/kinase_activity_substrate_psite_correlation/", height = 6, width = 3)
ggsave(filename = "kinase_activities_substrate_psite_based_correlation_box2.pdf", plot = ka_cor_dist, path = "./output/plots/kinase_activity_substrate_psite_correlation/", height = 6, width = 3)
unlink("kinase_activities_substrate_psite_based_correlation_box2.png")
unlink("kinase_activities_substrate_psite_based_correlation_box2.pdf")


ka_cor_dist <- bind_rows(ka_corByKin, ka_corByKin_all) %>%
  mutate(class = fct_recode(class, Others = "All pairs")) %>%
  ggplot() +
  geom_histogram(mapping = aes(x = pearson_r, y = stat(density), fill = class), position = "identity", alpha = 0.8) +
  theme_classic() +
  scale_fill_manual(values = c("#bdbdbd", "#67a9cf")) +
  theme(
    axis.text = element_text(color = "black", size = 8),
    legend.text = element_text(color = "black", size = 8),
    axis.title = element_text(color = "black", size = 10),
    legend.title = element_text(color = "black", size = 10)) +
  labs(x = "Pearson's r (substrate vs phosphosite-based kinase activity)", y = "Density", fill = "Class")

#+ fig.width=6, fig.height=2
ka_cor_dist

ggsave(filename = "kinase_activities_substrate_psite_based_correlation_hist.png", plot = ka_cor_dist, path = "./output/plots/kinase_activity_substrate_psite_correlation/", height = 2, width = 6)
ggsave(filename = "kinase_activities_substrate_psite_based_correlation_hist.pdf", plot = ka_cor_dist, path = "./output/plots/kinase_activity_substrate_psite_correlation/", height = 2, width = 6)
unlink("kinase_activities_substrate_psite_based_correlation_hist.png")
unlink("kinase_activities_substrate_psite_based_correlation_hist.pdf")
