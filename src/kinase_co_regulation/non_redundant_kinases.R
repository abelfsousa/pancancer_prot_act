library(tidyverse)

source(file = "./data/Danish/filterRedundantTerms.R")


# load kinase-substrate list
ks <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz") %>%
  select(-pair, -source) %>%
  distinct() %>%
  mutate(source_type = str_replace(source_type, "database|text-mining", "DB_text-mining")) %>%
  distinct() %>%
  filter(source_type == "DB_text-mining") %>%
  select(kinase, substrate) %>%
  distinct() %>%
  group_by(kinase) %>%
  summarise(substrate = list(substrate), n = n()) %>%
  ungroup() %>%
  filter(n >= 3) %>%
  select(-n)

ji <- function(x, y){
  ji <- length(intersect(x,y))/length(union(x,y))
  ji
}

k_pairs <- expand_grid(Ka = ks$kinase, Kb = ks$kinase) %>%
  inner_join(ks, by = c("Ka" = "kinase")) %>%
  rename(Ka_sub = substrate) %>%
  inner_join(ks, by = c("Kb" = "kinase")) %>%
  rename(Kb_sub = substrate) %>%
  mutate(ji = map2_dbl(.x = Ka_sub, .y = Kb_sub, .f = ji)) %>%
  select(-Ka_sub, -Kb_sub)

k_pairs_mat <- k_pairs %>%
  pivot_wider(names_from = "Kb", values_from = "ji") %>%
  column_to_rownames(var = "Ka") %>%
  as.matrix()

k_pairs_cluster <- hclust(as.dist(1-k_pairs_mat))
plot(k_pairs_cluster, cex = 0.5)


# -- use Danish method to remove kinase redundancies
# load kinase-substrate list
ks2 <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz") %>%
  select(-pair, -source) %>%
  distinct() %>%
  mutate(source_type = str_replace(source_type, "database|text-mining", "DB_text-mining")) %>%
  distinct() %>%
  filter(source_type == "DB_text-mining") %>%
  select(kinase, substrate) %>%
  distinct() %>%
  group_by(kinase) %>%
  summarise(substrate = list(substrate), n = n()) %>%
  ungroup() %>%
  filter(n >= 3) %>%
  select(-n) %>%
  unnest(cols = "substrate") %>%
  rename(ont=kinase, gene=substrate) %>%
  as.matrix()


# select non-redundant kinases using different tree height cutoffs
nr_kinases_0_8 <- pruneGOTerms(ks2, category = "broad", numTerms = 1, treeHeight = 0.8) %>%
  as_tibble() %>%
  rename(kinase = value)

nr_kinases_0_9 <- pruneGOTerms(ks2, category = "broad", numTerms = 1, treeHeight = 0.9) %>%
  as_tibble() %>%
  rename(kinase = value)

nr_kinases_0_99 <- pruneGOTerms(ks2, category = "broad", numTerms = 1, treeHeight = 0.99) %>%
  as_tibble() %>%
  rename(kinase = value)


write_tsv(nr_kinases_0_8, "./output/files/non_redundant_kinases_cutoff_0_8.txt")
write_tsv(nr_kinases_0_9, "./output/files/non_redundant_kinases_cutoff_0_9.txt")
write_tsv(nr_kinases_0_99, "./output/files/non_redundant_kinases_cutoff_0_99.txt")
