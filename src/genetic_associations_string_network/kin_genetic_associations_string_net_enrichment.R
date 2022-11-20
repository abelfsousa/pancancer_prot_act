library(tidyverse)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


sortPrt <- function(x, y){
  sorted <- sort(c(x,y))
  tibble(a=sorted[1],b=sorted[2])
}


# load mutation-kinase activity associations
associations <- read_tsv(file = "./output/files/ka_mutStatus_allGenes.txt.gz") %>%
  mutate(x = map2(kinase, gene, sortPrt)) %>%
  unnest() %>%
  select(a, b, everything())


# load string network (vertices names sorted for each edge; duplicates removed)
string_net <- read_tsv("./output/files/string_network_sorted_pairs.txt.gz")


# calculate string network enrichment across multiple edge weight cutoffs
enrichment <- function(cutoff, net, associations, beta, pval){
  print(cutoff)
  
  net <- net %>%
    filter(weight >= cutoff)
  
  interest_set <- associations %>%
    filter(p.adjust < pval & abs(estimate) > beta) %>%
    nrow()
  
  not_interest_set <- associations %>%
    filter(p.adjust > pval | abs(estimate) < beta) %>%
    nrow()
  
  overlap <- associations %>%
    semi_join(net, by = c("a", "b"))
  
  interest_inNet <- overlap %>%
    filter(p.adjust < pval & abs(estimate) > beta) %>%
    nrow()
  
  not_interest_inNet <- overlap %>%
    filter(p.adjust > pval | abs(estimate) < beta) %>%
    nrow()
  
  mat <- data.frame(
    interest = c(interest_inNet, interest_set-interest_inNet),
    not_interest = c(not_interest_inNet, not_interest_set-not_interest_inNet),
    row.names = c("in_net", "out_net"))
  
  ft <- fisher.test(mat, alternative = "greater")
  
  p <- ft$p.value
  
  p
}


net_enrichment <- tibble(weight = seq(150, 999, 1)) %>%
  mutate(pvalue = map_dbl(
    .x = weight,
    .f = enrichment,
    net = string_net,
    associations = associations,
    beta = 0,
    pval = 0.1))

plot <- net_enrichment %>%
  ggplot(mapping = aes(x = weight, y = -log10(pvalue))) +
  geom_line(size = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.3) +
  theme_minimal() +
  scale_y_continuous(limits = c(0,4)) +
  scale_x_continuous(breaks = seq(150, 999, 100)) +
  labs(x = "String network edge weight cutoff")

ggsave(filename = "kin_genetic_associations_string_enrichment_across_weights.png", plot = plot, path = "./output/plots/genetic_associations_string_network/", width = 4, height = 2)
ggsave(filename = "kin_genetic_associations_string_enrichment_across_weights.pdf", plot = plot, path = "./output/plots/genetic_associations_string_network/", width = 4, height = 2)
unlink("kin_genetic_associations_string_enrichment_across_weights.png")
unlink("kin_genetic_associations_string_enrichment_across_weights.pdf")
