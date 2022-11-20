# assess the protein abundance distribution between different mutation types

# load R packages
library(ggpubr)
library(rstatix)
library(viridis)
library(tidyverse)

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load NCI60/CRC65 cell lines
nci60_crc65 <- read_tsv(file = "./output/files/nci60_crc65_cell_lines.txt")

shared_cell_lines <- nci60_crc65 %>%
  group_by(batch) %>%
  summarise(cell_line = list(cell_line)) %>%
  pull(cell_line) %>%
  reduce(intersect)


# load NCI60/CRC65 mutations
# select only the gene, sample and type of mutation
# remove duplicates (same gene in the same sample can have multiple mutations of same type)
nci60_crc65_mutations <- read_tsv(file = "./output/files/nci60_crc65_mutations.txt.gz") %>%
  select(sample, gene, mutation_type) %>%
  distinct()


# load protein abundance inference data from NCI60/CRC65 dataset
protein <- read_tsv(file = "./output/files/nci60_crc65_protein_log2fc.txt.gz") %>%
  filter(!(batch == "CRC65" & sample %in% shared_cell_lines)) %>%
  select(sample, gene, log2fc) %>%
  filter(!is.na(log2fc))


# filter the mutations that occurred in proteins
# select only the pairs sample-protein that are unique to the same type of mutation
# this was done to prevent the assignment of the same protein abundance to different types of mutations
mutations <- nci60_crc65_mutations %>%
  semi_join(protein, by = c("gene", "sample")) %>%
  group_by(mutation_type) %>%
  nest() %>%
  ungroup() %>%
  mutate(data2 = (function(x) map2(.x=1:length(x),.y=x,.f=function(a,b)reduce(.x=c(list(b),x[-a]),.f=setdiff))) (data)) %>%
  select(-data) %>%
  unnest()

# define a background to compare with (no known mutation for the protein in the sample)
background <- protein %>%
  anti_join(nci60_crc65_mutations, by = c("gene", "sample")) %>%
  mutate(mutation_type = "Background (no mutation)")


protein_mutation <- protein %>%
  inner_join(mutations, by = c("sample", "gene")) %>%
  bind_rows(background)

filter_abundances <- function(vect, lims){
  x <- vect
  x <- x[x >= min(lims) & x <= max(lims)]
  
  return(x)
}

protein_mutation2 <- protein_mutation %>%
  group_by(mutation_type) %>%
  summarise(log2fc = list(log2fc)) %>%
  ungroup() %>%
  mutate(limits = map(.x = log2fc, .f = ~ c(quantile(.x, 0.25)-(1.5*IQR(.x)), quantile(.x, 0.75)+(1.5*IQR(.x))))) %>%
  mutate(log2fc = map2(.x = log2fc, .y = limits, .f = filter_abundances)) %>%
  select(-limits) %>%
  unnest()

pvalues_ <- protein_mutation2 %>%
  t_test(log2fc ~ mutation_type, ref.group = "Background (no mutation)", p.adjust.method = "none") %>%
  add_xy_position(x = "mutation_type")

pvalues <- protein_mutation %>%
  t_test(log2fc ~ mutation_type, ref.group = "Background (no mutation)", p.adjust.method = "none") %>%
  add_xy_position(x = "mutation_type") %>%
  inner_join(pvalues_[, c("group2", "y.position")], by = c("group2")) %>%
  rename(y.position = y.position.y)

N <- protein_mutation %>%
  group_by(mutation_type) %>%
  summarise(n = n()) %>%
  ungroup()

N_ <- protein_mutation2 %>%
  group_by(mutation_type) %>%
  summarise(min = min(log2fc)) %>%
  mutate(min = min(min)+0.15) %>%
  ungroup()

N <- N %>% inner_join(N_)


distribution <- protein_mutation2 %>%
  ggplot(mapping = aes(x=mutation_type, y = log2fc)) +
  geom_boxplot(mapping = aes(fill = mutation_type), alpha = 0.8, show.legend = F, outlier.shape = NA) +
  stat_pvalue_manual(pvalues, label = "p", tip.length = 0.01, coord.flip = TRUE, size = 6, step.increase = 0.02) +
  geom_text(data = N, mapping = aes(x = mutation_type, y = min, label = n), size = 7) +
  coord_flip() +
  scale_fill_viridis(discrete = T) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 24, colour = "black", hjust = 0.5),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 24, colour = "black"),
    axis.text.x = element_text(size = 22, colour = "black"),
    axis.text.y = element_text(size = 18, colour = "black")) +
  labs(y = "Protein abundance (log2fc)", title = "NCI60/CRC65 cancer cell lines")
distribution

ggsave(filename = "nci60_crc65_prot_abundance_mutation_distribution.png", plot = distribution, path = "./output/plots/nci60_crc65/", height = 6, width = 12)
ggsave(filename = "nci60_crc65_prot_abundance_mutation_distribution.pdf", plot = distribution, path = "./output/plots/nci60_crc65/", height = 6, width = 12)
