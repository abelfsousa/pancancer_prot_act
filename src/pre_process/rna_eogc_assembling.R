library(tidyverse)

source("src/utils/merge_files.R")
source("src/utils/remove_dup_genes.R")


# assembling of EOGC rna-seq data from "Proteogenomic Characterization of Human Early-Onset Gastric Cancer"
# quantile normalized RSEM expression count values


# eogc samples with protein measurements
eogc_proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt") %>%
  filter(batch == "cptac-GC") %>%
  mutate(proteomic_ID = str_extract(sample, "T[0-9]{1,4}")) %>%
  select(sample, proteomic_ID, batch, cancer)


# gene list
gene_list <- read_tsv("/Volumes/G-DRIVE/data/geo/Proteogenomics_of_Gastric_Cancer/GSE122401_RAW/gene_ids.txt")


# files dir
eogc_dir <- "/Volumes/G-DRIVE/data/geo/Proteogenomics_of_Gastric_Cancer/GSE122401_RAW/"


# sample IDs mapping
sample_map <- read_tsv("/Volumes/G-DRIVE/data/cptac/Proteogenomics_of_Gastric_Cancer/clinical/sample_IDs_mapping.txt") %>%
  filter(str_detect(proteomic_ID, "T")) %>%
  inner_join(eogc_proteomics[, c("sample", "proteomic_ID")], by = c("proteomic_ID"))


# file names
eogc_expr_files <- read_tsv("/Volumes/G-DRIVE/data/geo/Proteogenomics_of_Gastric_Cancer/GSE122401_RAW/gene_expr_files_tumours.txt", col_names = F) %>%
  rename(file = X1) %>%
  mutate(geo_ID = str_split(string = file, pattern = c("_"), n = 2, simplify = T)[,1])


# build gene expression table
eogc_rsem_counts <- merge_files_gc(eogc_expr_files, eogc_dir, gene_list) %>%
  column_to_rownames(var = "gene_id") %>%
  round() %>%
  rownames_to_column(var = "gene_id") %>%
  as_tibble() %>%
  rename(gene = gene_id)



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

eogc_rsem_counts <- eogc_rsem_counts %>%
  filter(gene %in% protein_coding$gene_id)


# obtain gene lengths
#gene_lengths <- EDASeq::getGeneLengthAndGCContent(eogc_rsem_counts$gene, "hsa", mode=c("biomart")) %>%
#  as.data.frame() %>%
#  rownames_to_column(var = "gene") %>%
#  as_tibble() %>%
#  filter(!is.na(length))

#http://www.genemine.org/gtftools.php
#https://www.biorxiv.org/content/10.1101/263517v1
gene_lengths <- read_tsv("./output/files/gencode.v19.gene_lengths.txt") %>%
  mutate(gene = str_split_fixed(gene, "\\.", 2)[,1]) %>%
  select(gene, length = merged)



# calculate FPKM
eogc_rsem_fpkm <- eogc_rsem_counts %>%
  inner_join(gene_lengths, by = c("gene")) %>%
  dplyr::select(gene, length, everything())

eogc_rsem_fpkm <- edgeR::rpkm(y = column_to_rownames(eogc_rsem_fpkm[, -c(2)], var = "gene"), gene.length = eogc_rsem_fpkm$length, log = F) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()


# get gene symbol IDs
eogc_rsem_fpkm <- eogc_rsem_fpkm %>%
  inner_join(protein_coding, by = c("gene" = "gene_id")) %>%
  select(-gene) %>%
  select(gene = gene_name, everything())


# change sample IDs
eogc_rsem_fpkm <- eogc_rsem_fpkm %>%
  gather(key = "geo_ID", value = "fpkm", -gene) %>%
  inner_join(sample_map[, c("geo_ID", "sample")], by = "geo_ID") %>%
  select(-geo_ID) %>%
  spread(key = "sample", value = "fpkm")


gz1 <- gzfile("./output/files/transcriptomics_eogc_fpkm.txt.gz", "w")
write.table(eogc_rsem_fpkm, gz1, sep="\t", quote=F, row.names=F)
close(gz1)


#eogc samples
eogc_samples <- tibble(sample = colnames(eogc_rsem_fpkm)[-c(1)], batch = "GC", cancer = "eogc")

write.table(eogc_samples, "./output/files/transcriptomics_samples_eogc.txt", sep="\t", quote=F, row.names=F)
