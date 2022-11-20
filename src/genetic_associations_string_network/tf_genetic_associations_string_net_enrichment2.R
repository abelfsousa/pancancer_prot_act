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
associations <- read_tsv(file = "./output/files/tf_mutStatus_allGenes.txt.gz") %>%
  filter(tf != gene) %>%
  mutate(p.adjust = p.adjust(p.value, "BH")) %>%
  mutate(x = map2(tf, gene, sortPrt)) %>%
  unnest() %>%
  select(a, b, everything())


# load string network (vertices names sorted for each edge; duplicates removed)
string_net <- read_tsv("./output/files/string_network_sorted_pairs.txt.gz")


# perform a fisher test on a specific association cutoff
ftest <- function(pairs, network, pval, beta){
  print(str_c("pvalue:", pval, "beta:", beta, sep = " "))
  
  interest_set <- pairs %>%
    filter(p.adjust < pval & abs(estimate) > beta) %>%
    nrow()
  
  not_interest_set <- pairs %>%
    filter(p.adjust > pval | abs(estimate) < beta) %>%
    nrow()
  
  overlap <- pairs %>%
    semi_join(network, by = c("a", "b"))
  
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


# filter the string network using interaction cutoffs
# calculate enrichment pvalues across multiple association cutoffs
enrichment <- function(wt, network, cutoffs, associations, by, fixed_value){
  print(wt)
  net <- network %>%
    filter(weight >= wt)
  
  pvalues <- tibble(cutoff = cutoffs)
  
  if(by == "pvalue"){
    pvalues <- pvalues  %>%
      mutate(pvalue = map_dbl(
        .x = cutoff,
        .f = ftest,
        pairs = associations,
        network = net,
        beta = fixed_value))
  } else if(by == "beta"){
    pvalues <- pvalues  %>%
      mutate(pvalue = map_dbl(
        .x = cutoff,
        .f = ftest,
        pairs = associations,
        network = net,
        pval = fixed_value))
  } else {
    stop("'by' not recognized. Please use 'pvalue' or 'beta'")
  }
  pvalues
}


# set of P-value cutoffs
x <- 1e-10
while(x[length(x)] < 1e-3){
  v <- x[length(x)]
  v <- v*10
  x <- c(x, v)
}
ctoffs <- c(rev(seq(0.1,1,0.1)), rev(seq(0.01,0.09,0.01)), rev(x))

net_enrichment <- tibble(weight = seq(150, 950, 100)) %>%
  mutate(pvalue = map(
    .x = weight,
    .f = enrichment,
    network = string_net,
    cutoffs = ctoffs,
    associations = associations,
    by = "pvalue",
    fixed_value = 0)) %>%
  unnest()

plot <- net_enrichment %>%
  ggplot(mapping = aes(x = -log10(cutoff), y = -log10(pvalue))) +
  geom_line(size = 0.3) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.3) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", size = 0.3) +
  facet_wrap(~ weight) +
  theme_minimal() +
  theme(panel.spacing.x = unit(0.5, "lines")) +
  labs(x = "Association adjusted P-value cutoff (-log10)", y = "Fisher-test P-value (-log10)")

ggsave(filename = "tf_genetic_associations_string_enrichment_across_pvalues.png", plot = plot, path = "./output/plots/genetic_associations_string_network/", width = 8, height = 4)
ggsave(filename = "tf_genetic_associations_string_enrichment_across_pvalues.pdf", plot = plot, path = "./output/plots/genetic_associations_string_network/", width = 8, height = 4)
unlink("tf_genetic_associations_string_enrichment_across_pvalues.png")
unlink("tf_genetic_associations_string_enrichment_across_pvalues.pdf")


# set of P-value cutoffs
x <- 1e-10
while(x[length(x)] < 1e-3){
  v <- x[length(x)]
  v <- v*10
  x <- c(x, v)
}
y = seq(1e-05, 1e-04, by=1e-05)
#ctoffs <- c(rev(seq(0.2,1,0.2)), rev(seq(0.02,0.08,0.02)), rev(x))
#ctoffs <- c(1, 0.5, 0.05, 0.01, rev(x))
#ctoffs <- c(rev(seq(0.1,1,0.1)), rev(seq(0.01,0.08,0.01)), rev(x), y[c(-1, -length(y))])
ctoffs <- c(rev(seq(0.1,1,0.1)), rev(seq(0.01,0.08,0.01)), rev(x))

net_enrichment2 <- tibble(weight = seq(150, 950, 100)) %>%
  mutate(pvalue = map(
    .x = weight,
    .f = enrichment,
    network = string_net,
    cutoffs = ctoffs,
    associations = associations,
    by = "pvalue",
    fixed_value = 0)) %>%
  unnest()

plot <- net_enrichment2 %>%
  mutate(weight = weight/1000) %>%
  filter(weight == 0.85) %>%
  filter(cutoff >= 1e-05) %>%
  ggplot(mapping = aes(x = -log10(cutoff), y = -log10(pvalue))) +
  geom_line(size = 0.3, color = "#e31a1c") +
  geom_point(size = 2, color = "#e31a1c") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.2) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", size = 0.2) +
  #facet_wrap(~ weight) +
  theme_classic() +
  theme(
    plot.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = 14),
    axis.text = element_text(colour = "black", size = 16),
    axis.title = element_text(colour = "black", size = 18),
    panel.spacing.x = unit(0.4, "lines")) +
  labs(x = "Association adjusted P-value cutoff", y = "Fisher-test P-value")

ggsave(filename = "tf_genetic_associations_string_enrichment_across_pvalues2.png", plot = plot, path = "./output/plots/genetic_associations_string_network/", width = 5, height = 3)
ggsave(filename = "tf_genetic_associations_string_enrichment_across_pvalues2.pdf", plot = plot, path = "./output/plots/genetic_associations_string_network/", width = 5, height = 3)


# set of absolute beta cutoffs
ctoffs <- seq(0,6.5,0.1)

net_enrichment <- tibble(weight = seq(150, 950, 100)) %>%
  mutate(pvalue = map(
    .x = weight,
    .f = enrichment,
    network = string_net,
    cutoffs = ctoffs,
    associations = associations,
    by = "beta",
    fixed_value = 0.01)) %>%
  unnest()

plot <- net_enrichment %>%
  ggplot(mapping = aes(x = cutoff, y = -log10(pvalue))) +
  geom_line(size = 0.3) +
  geom_point(size = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", size = 0.3) +
  facet_wrap(~ weight) +
  theme_minimal() +
  theme(panel.spacing.x = unit(1, "lines")) +
  labs(x = expression("Association |" ~ beta ~ "| cutoff (mutation status)"), y = "Fisher-test P-value (-log10)")

ggsave(filename = "tf_genetic_associations_string_enrichment_across_betas.png", plot = plot, path = "./output/plots/genetic_associations_string_network/", width = 6, height = 4)
ggsave(filename = "tf_genetic_associations_string_enrichment_across_betas.pdf", plot = plot, path = "./output/plots/genetic_associations_string_network/", width = 6, height = 4)
unlink("tf_genetic_associations_string_enrichment_across_betas.png")
unlink("tf_genetic_associations_string_enrichment_across_betas.pdf")
