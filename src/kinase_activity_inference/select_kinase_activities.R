# select kinase activities

library(tidyverse)

# load kinase activity inference data
k_subN <- 3

KA1 <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-n, -source_type) %>%
  rename(kin_activity = log10P)
write_tsv(KA1, "./output/files/cptac_kin_activities.txt")

KA2 <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_ProtBatchRegout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-n, -source_type) %>%
  rename(kin_activity = log10P)
write_tsv(KA2, "./output/files/cptac_kin_activities_BatchRegOut.txt")

KA3 <- data.table::fread("./data/Danish/kinaseActMatImputed.tsv") %>%
  as_tibble() %>%
  rename(kinase = V1) %>%
  pivot_longer(-kinase, names_to = "sample", values_to = "kin_activity")

KA4 <- read_tsv("./output/files/kinaseActMatImputed_batch_regOut.txt") %>%
  pivot_longer(-kinase, names_to = "sample", values_to = "kin_activity")

KA5 <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withoutZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-n, -source_type) %>%
  rename(kin_activity = log10P)


inner_join(KA2, KA4, by = c("sample", "kinase")) %>%
  group_by(kinase) %>%
  summarise(corr = cor(kin_activity.x, kin_activity.y), n = n()) %>%
  ungroup() %>%
  filter(n>5) %>%
  ggplot(mapping = aes(x = factor(0), y = corr)) +
  geom_boxplot()
