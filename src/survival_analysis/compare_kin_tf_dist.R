#' load R packages
library(tidyverse)


#' load kinase-activity inference data\
#' (quantile-normalized protein regressed-out phosphorylation data)\
#' select quantifications with more than 3 substrates
k_subN <- 3
kin_act <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  select(protein=kinase, sample, activity=log10P) %>%
  mutate(protein_type = "Kinase")


#' load transcription factor activities
tf_act <- data.table::fread("./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  as_tibble() %>%
  rename(tf=V1) %>%
  pivot_longer(-tf, names_to = "sample", values_to = "tf_activity") %>%
  select(protein=tf, sample, activity=tf_activity) %>%
  mutate(protein_type = "TF")


prot_act <- bind_rows(kin_act, tf_act)


stats <- prot_act %>%
  group_by(protein_type, sample) %>%
  summarise(x = list(as_tibble(skimr::skim(activity)))) %>%
  ungroup() %>%
  unnest(cols = c(x)) %>%
  select(protein_type, sample, contains(c("p25", "p50", "p75"))) %>%
  pivot_longer(-c(protein_type, sample), names_to = "statistic", values_to = "value") %>%
  group_by(protein_type, statistic) %>%
  summarise(MEAN = mean(value), SD = sd(value)) %>%
  ungroup() %>%
  mutate(statistic = str_replace(statistic, "numeric.", ""))

# stats <- prot_act %>%
#   group_by(protein_type, sample) %>%
#   summarise(x = list(broom::tidy(quantile(activity)))) %>%
#   ungroup() %>%
#   unnest(cols = c(x)) %>%
#   rename(statistic = names, value = x) %>%
#   filter(!statistic %in% c("0%", "100%")) %>%
#   group_by(protein_type, statistic) %>%
#   summarise(MEAN = mean(value), SD = sd(value)) %>%
#   ungroup() %>%
#   mutate(statistic = str_replace(statistic, "%", "")) %>%
#   mutate(statistic = str_c("p", statistic, sep = ""))


qplot <- stats %>%
  ggplot(mapping = aes(x = statistic, y = MEAN, fill = protein_type)) + 
  geom_col(position = "dodge") +
  #geom_text(mapping = aes(x = statistic, y = sign(MEAN)*(abs(MEAN)+0.05), label = round(MEAN, 3)), position = position_dodge(width = 1), size = 3) +
  #geom_errorbar(mapping = aes(ymin = MEAN-(SD), ymax = MEAN+(SD)), size = 0.5, width = 0.2, position = position_dodge(width = 0.9), show.legend = F) +
  geom_pointrange(mapping = aes(ymin = MEAN-(SD), ymax = MEAN+(SD)), size = 0.1, position = position_dodge(width = 0.9), show.legend = F) +
  theme_classic() +
  theme(
    legend.position = "bottom") +
  scale_x_discrete(labels = c("p25" = "Q1", "p50" = "Q2", "p75" = "Q3"), name = "Quantile") +
  scale_y_continuous(name = "Mean", limits = c(-2,2)) +
  scale_fill_discrete(name = "Protein")

ggsave(filename = "kin_tf_activities_quantiles_cptac.png", plot = qplot, path = "./output/plots/survival_analysis/", width = 3, height = 4)
ggsave(filename = "kin_tf_activities_quantiles_cptac.pdf", plot = qplot, path = "./output/plots/survival_analysis/", width = 3, height = 4)


kinases <- prot_act %>%
  filter(protein_type == "Kinase") %>%
  pull(protein) %>%
  unique()

TFs <- prot_act %>%
  filter(protein_type == "TF") %>%
  pull(protein) %>%
  unique()


all_dist <- prot_act %>%
  mutate(x = if_else(protein_type == "Kinase", 1, 2)) %>%
  mutate(protein = fct_relevel(.f=protein, kinases)) %>%
  ggplot(mapping = aes(x = protein, fill = as.character(protein_type), y = activity)) +
  geom_boxplot(outlier.size = 0.05, size = 0.1) + 
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom") +
  labs(x = "Protein", y = "Activity", fill = "Protein type")

ggsave(filename = "kin_tf_activities_all_dists_cptac.png", plot = all_dist, path = "./output/plots/survival_analysis/", width = 10, height = 4)
ggsave(filename = "kin_tf_activities_all_dists_cptac.pdf", plot = all_dist, path = "./output/plots/survival_analysis/", width = 10, height = 4)

