library(tidyverse)


# assembling of cnv data from tcga samples
# GISTIC2 scores from cbioportal



# cancer samples with proteomics measurements
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt") %>%
  filter(str_detect(batch, "tcga"))


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


# BRCA
brca_cnv <- read_tsv("./data/dna/cnv/cbioportal/brca_tcga_pan_can_atlas_2018/data_CNA.txt") %>%
  rename(gene=Hugo_Symbol) %>%
  select(-Entrez_Gene_Id) %>%
  filter(gene %in% protein_coding$gene_name) %>%
  group_by(gene) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n) %>%
  rename_at(.vars = 2:ncol(.), .funs = ~ substring(.x, 1, 12)) %>%
  rename_at(.vars = 2:ncol(.), .funs = ~ str_replace_all(.x, "-", "."))
  #select_if(.predicate = colnames(.) %in% c("gene", proteomics$sample))

# samples
brca_samples <- tibble(sample = colnames(brca_cnv)[-c(1)], batch = "tcga-brca", cancer = "brca")



coread_cnv <- read_tsv("./data/dna/cnv/cbioportal/coadread_tcga_pan_can_atlas_2018/data_CNA.txt") %>%
  rename(gene=Hugo_Symbol) %>%
  select(-Entrez_Gene_Id) %>%
  filter(gene %in% protein_coding$gene_name) %>%
  group_by(gene) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n) %>%
  rename_at(.vars = 2:ncol(.), .funs = ~ substring(.x, 1, 12)) %>%
  rename_at(.vars = 2:ncol(.), .funs = ~ str_replace_all(.x, "-", "."))
  #select_if(.predicate = colnames(.) %in% c("gene", proteomics$sample))

# samples
coread_samples <- tibble(sample = colnames(coread_cnv)[-c(1)], batch = "tcga-coread", cancer = "coread")



ov_cnv <- read_tsv("./data/dna/cnv/cbioportal/ov_tcga_pan_can_atlas_2018/data_CNA.txt") %>%
  rename(gene=Hugo_Symbol) %>%
  select(-Entrez_Gene_Id) %>%
  filter(gene %in% protein_coding$gene_name) %>%
  rename_at(.vars = 2:ncol(.), .funs = ~ substring(.x, 1, 12)) %>%
  rename_at(.vars = 2:ncol(.), .funs = ~ str_replace_all(.x, "-", "."))
  #select_if(.predicate = colnames(.) %in% c("gene", proteomics$sample))

# samples
ov_samples <- tibble(sample = colnames(ov_cnv)[-c(1)], batch = "tcga-hgsc", cancer = "hgsc")



# join sample batch and cancer together
tcga_samples <- bind_rows(
  brca_samples,
  coread_samples,
  ov_samples
)

write.table(tcga_samples, "./output/files/cnv_samples_tcga.txt", sep="\t", quote=F, row.names=F)



# join all datasets together by gene
# keep only the genes in common to the three datasets

tcga_cnv <- inner_join(brca_cnv, coread_cnv, by = c("gene")) %>%
  inner_join(ov_cnv, by = c("gene")) %>%
  select_if(.predicate = colnames(.) %in% c("gene", proteomics$sample))

gz1 <- gzfile("./output/files/cnv_tcga.txt.gz", "w")
write.table(tcga_cnv, gz1, sep="\t", quote=F, row.names=F)
close(gz1)

