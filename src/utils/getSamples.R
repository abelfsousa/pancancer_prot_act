library(tidyverse)

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy

getSamples <- function(all_samples, data_types){
  
  l = length(data_types)
  
  samples <- all_samples %>%
    filter(data %in% data_types) %>%
    rename(info = data) %>%
    group_by(info, batch, cancer, tissue) %>%
    summarise(sample = list(sample)) %>%
    ungroup() %>%
    group_by(batch, cancer, tissue) %>%
    summarise(sample = list(sample)) %>%
    ungroup() %>%
    filter(map_dbl(.x=sample, .f=length) == l) %>%
    mutate(overlap = map(.x=sample, .f = ~ reduce(.x=.x, .f=intersect))) %>%
    select(-sample) %>%
    unnest() %>%
    rename(sample = overlap)
  
  return(samples)
  
}


# all_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")
# 
# samples <- all_samples %>%
#   filter(data %in% c("phosphorylation", "protein")) %>%
#   rename(info = data) %>%
#   group_by(info, batch, cancer, tissue) %>%
#   summarise(sample = list(sample)) %>%
#   ungroup() %>%
#   group_by(batch, cancer, tissue) %>%
#   summarise(sample = list(sample)) %>%
#   ungroup() %>%
#   filter(map_dbl(.x=sample, .f=length) == 2) %>%
#   mutate(overlap = map(.x=sample, .f = ~ reduce(.x=.x, .f=intersect))) %>%
#   select(-sample) %>%
#   unnest() %>%
#   rename(sample = overlap)
# 
# samples <- all_samples %>%
#   filter(data %in% c("phosphorylation", "protein")) %>%
#   rename(info = data) %>%
#   group_by(info, batch, cancer, tissue) %>%
#   summarise(sample = list(sample)) %>%
#   ungroup() %>%
#   pivot_wider(names_from = "info", values_from = "sample") %>%
#   mutate(overlap = map2(.x=phosphorylation, .y=protein, .f = ~ intersect(.x,.y))) %>%
#   select(-phosphorylation, -protein) %>%
#   unnest() %>%
#   rename(sample = overlap)
