library(tidyverse)
library(estimate)
library(cowplot)


samples_metadata <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
  filter(data == "protein") %>%
  select(-data)


rna <- read_tsv(file = "./output/files/transcriptomics_log2fpkm.txt.gz")

rna <- rna %>%
  rename(NAME = gene) %>%
  mutate(Description = NAME) %>%
  select(NAME, Description, everything())

write.table(rna, "./output/files/transcriptomics_log2fpkm_ESTIMATE_input.gct", sep="\t", quote=F, row.names=F)


estimateScore(
  input.ds = "./output/files/transcriptomics_log2fpkm_ESTIMATE_input.gct",
  output.ds = "./output/files/samples_purity.gct",
  platform = "affymetrix")

purity <- read_tsv("./output/files/samples_purity.gct", skip = 2) %>%
  select(-Description) %>%
  set_names(nm = c("NAME", colnames(rna)[3:ncol(rna)])) %>%
  pivot_longer(-NAME, names_to = "sample", values_to = "score") %>%
  rename(metric = NAME) %>%
  filter(metric == "TumorPurity")

write_tsv(purity, "./output/files/samples_purity.txt")


purity <- purity %>%
  inner_join(samples_metadata, by = "sample")


purity_density1 <- ggplot(data = purity, mapping = aes(x = score)) +
  geom_density() +
  labs(x = "Tumour purity score", title = "All samples")

purity_density2 <- ggplot(data = purity, mapping = aes(x = score)) +
  geom_density() +
  facet_wrap(~ batch, scales = "fixed") +
  theme(axis.text = element_text(size = 7)) +
  labs(x = "Tumour purity score", title = "By study")

ggsave(filename = "purity_score_distribution_all_samples.png", plot = purity_density1, path = "./output/plots/purity_analysis/", height = 4, width = 4)
unlink("purity_score_distribution_all_samples.png")

ggsave(filename = "purity_score_distribution_by_study.png", plot = purity_density2, path = "./output/plots/purity_analysis/", height = 6, width = 6)
unlink("purity_score_distribution_by_study.png")

