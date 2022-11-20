library(tidyverse)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


sortPrt <- function(x, y){
  sorted <- sort(c(x,y))
  tibble(a=sorted[1],b=sorted[2])
}


# load string network gene to protein mapping IDs
string_net_map <- read_tsv(file = "./data/string/9606.protein.info.v11.5.txt.gz") %>%
  rename(protein_external_id = `#string_protein_id`)


# load string network
string_net <- data.table::fread(file = "./data/string/9606.protein.physical.links.v11.5.txt.gz") %>%
  as_tibble() %>%
  inner_join(string_net_map[, c("protein_external_id", "preferred_name")], by = c("protein1" = "protein_external_id")) %>%
  select(-protein1) %>%
  rename(a=preferred_name) %>%
  inner_join(string_net_map[, c("protein_external_id", "preferred_name")], by = c("protein2" = "protein_external_id")) %>%
  select(-protein2) %>%
  rename(b=preferred_name) %>%
  mutate(x = map2(a, b, sortPrt)) %>%
  select(-a, -b) %>%
  unnest() %>%
  select(a, b, weight=combined_score) %>%
  distinct()

write_tsv(string_net, "./output/files/string_network_sorted_pairs_physical.txt.gz")
