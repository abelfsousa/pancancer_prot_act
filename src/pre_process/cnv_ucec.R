library(tidyverse)


# assembling of cnv data from "Proteogenomic Characterization of Endometrial Carcinoma"
# log2 copy-number ratios from linkedomics


# cancer samples with proteomics measurements
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt") %>%
  filter(str_detect(batch, "ucec"))


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


ucec <- read_tsv(file = "./data/dna/cnv/linkedomics/ucec/all_thresholded.by_genes.txt") %>%
  rename(gene=`Gene Symbol`) %>%
  select(-`Locus ID`, -`Cytoband`) %>%
  rename_all(.funs = ~ str_replace(.x, "-", ".")) %>%
  filter(gene %in% protein_coding$gene_name)


# samples
samples <- tibble(sample = colnames(ucec)[-c(1)], batch = "discovery-ucec", cancer = "ucec")
write.table(samples, "./output/files/cnv_samples_ucec.txt", sep="\t", quote=F, row.names=F)


ucec <- ucec %>%
  select_if(.predicate = colnames(.) %in% c("gene", proteomics$sample))

gz1 <- gzfile("./output/files/cnv_ucec.txt.gz", "w")
write.table(ucec, gz1, sep="\t", quote=F, row.names=F)
close(gz1)
