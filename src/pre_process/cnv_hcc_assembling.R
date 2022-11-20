library(tidyverse)


# assembling of cnv data from "Integrated Proteogenomic Characterization of HBV-Related Hepatocellular Carcinoma"
# GISTIC2 scores from NODE



# cancer samples with proteomics measurements
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt") %>%
  filter(str_detect(batch, "hcc"))


# protein coding genes
protein_coding <- read_tsv("./output/files/gencode.v19.annotation.txt") %>%
  select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  select(-gene_type) %>%
  group_by(gene_name) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)


hcc <- read_tsv(file = "./data/dna/cnv/node/hcc/all_thresholded.by_genes.txt") %>%
  select(-c(`Locus ID`, Cytoband)) %>%
  rename(gene=`Gene Symbol`) %>%
  filter(gene %in% protein_coding$gene_name)


# samples
samples <- tibble(sample = colnames(hcc)[-c(1)], batch = "hcc", cancer = "hcc")
write.table(samples, "./output/files/cnv_samples_hcc.txt", sep="\t", quote=F, row.names=F)


hcc <- hcc %>%
  select_if(.predicate = colnames(.) %in% c("gene", proteomics$sample))

gz1 <- gzfile("./output/files/cnv_hcc.txt.gz", "w")
write.table(hcc, gz1, sep="\t", quote=F, row.names=F)
close(gz1)

