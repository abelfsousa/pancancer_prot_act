library(tidyverse)
library(data.table)

source("src/utils/merge_files.R")
source("src/utils/remove_dup_genes.R")



# assembling of rna-seq data from tcga studies


# gene list
gene_list <- read_tsv("/Volumes/G-DRIVE/data/gdc/tcga/rnaseq/gene_annotation/gene_ids.txt")



# gene annotation
gene_annotation <- read_tsv("/Volumes/G-DRIVE/data/gdc/tcga/rnaseq/gene_annotation/geneAnnot.gencode.v22.txt")


# protein coding genes
protein_coding <- gene_annotation %>%
  filter(gene_type == "protein_coding") %>%
  filter(str_detect(gene_id, "ENSGR", negate = T)) %>%
  select(1:2) %>%
  group_by(gene_name) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)


# cancer samples with proteomics measurements
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt")


# cancer samples with phosphoproteomics measurements
phospho <- read_tsv(file = "./output/files/phosphoproteomics_samples.txt")



# BRCA

brca_dir <- "/Volumes/G-DRIVE/data/gdc/tcga/rnaseq/brca/files/"

# HTSeq counts
brca_file_htseq_counts <- fread(file = "/Volumes/G-DRIVE/data/gdc/tcga/rnaseq/brca/metadata/gdc_sample_sheet.tumours.htseq_counts.tsv", check.names = T) %>%
  as_tibble()

brca_htseq_counts <- merge_files(brca_file_htseq_counts, brca_dir, gene_list)

#keep duplicated samples (coming from multiple files for the same sample)
brca_htseq_counts <- brca_htseq_counts %>%
  inner_join(protein_coding[, c("gene_id", "gene_name")], by =  c("gene" = "gene_id")) %>%
  select(-gene) %>%
  select(gene = gene_name, everything())


# HTSeq FPKM
brca_file_htseq_fpkm <- fread(file = "/Volumes/G-DRIVE/data/gdc/tcga/rnaseq/brca/metadata/gdc_sample_sheet.tumours.htseq_fpkm.tsv", check.names = T) %>%
  as_tibble()

brca_htseq_fpkm <- merge_files(brca_file_htseq_fpkm, brca_dir, gene_list)

#average duplicated samples by gene (coming from multiple files for the same sample)
brca_htseq_fpkm <- brca_htseq_fpkm %>%
  inner_join(protein_coding[, c("gene_id", "gene_name")], by =  c("gene" = "gene_id")) %>%
  select(-gene) %>%
  select(gene = gene_name, everything()) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  mutate(sample = str_replace(sample, ".x|.y", "")) %>%
  group_by(sample, gene) %>%
  summarise(fpkm = mean(fpkm)) %>%
  ungroup() %>%
  spread(key = "sample", value = "fpkm")

#brca samples
brca_samples <- tibble(sample = colnames(brca_htseq_fpkm)[-c(1)], batch = "tcga-brca", cancer = "brca")




# COAD

coad_dir <- "/Volumes/G-DRIVE/data/gdc/tcga/rnaseq/coad/files/"

# HTSeq counts
coad_file_htseq_counts <- fread(file = "/Volumes/G-DRIVE/data/gdc/tcga/rnaseq/coad/metadata/gdc_sample_sheet.tumours.htseq_counts.tsv", check.names = T) %>%
  as_tibble()

coad_htseq_counts <- merge_files(coad_file_htseq_counts, coad_dir, gene_list)

#keep duplicated samples (coming from multiple files for the same sample)
coad_htseq_counts <- coad_htseq_counts %>%
  inner_join(protein_coding[, c("gene_id", "gene_name")], by =  c("gene" = "gene_id")) %>%
  select(-gene) %>%
  select(gene = gene_name, everything())


# HTSeq FPKM
coad_file_htseq_fpkm <- fread(file = "/Volumes/G-DRIVE/data/gdc/tcga/rnaseq/coad/metadata/gdc_sample_sheet.tumours.htseq_fpkm.tsv", check.names = T) %>%
  as_tibble()

coad_htseq_fpkm <- merge_files(coad_file_htseq_fpkm, coad_dir, gene_list)

#average duplicated samples by gene (coming from multiple files for the same sample)
coad_htseq_fpkm <- coad_htseq_fpkm %>%
  inner_join(protein_coding[, c("gene_id", "gene_name")], by =  c("gene" = "gene_id")) %>%
  select(-gene) %>%
  select(gene = gene_name, everything()) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  mutate(sample = str_replace(sample, ".x|.y", "")) %>%
  group_by(sample, gene) %>%
  summarise(fpkm = mean(fpkm)) %>%
  ungroup() %>%
  spread(key = "sample", value = "fpkm")

#coad samples
coad_samples <- tibble(sample = colnames(coad_htseq_fpkm)[-c(1)], batch = "tcga-coread", cancer = "coread")




# READ

read_dir <- "/Volumes/G-DRIVE/data/gdc/tcga/rnaseq/read/files/"

# HTSeq counts
read_file_htseq_counts <- fread(file = "/Volumes/G-DRIVE/data/gdc/tcga/rnaseq/read/metadata/gdc_sample_sheet.tumours.htseq_counts.tsv", check.names = T) %>%
  as_tibble()

read_htseq_counts <- merge_files(read_file_htseq_counts, read_dir, gene_list)

read_htseq_counts <- read_htseq_counts %>%
  inner_join(protein_coding[, c("gene_id", "gene_name")], by =  c("gene" = "gene_id")) %>%
  select(-gene) %>%
  select(gene = gene_name, everything())


# HTSeq FPKM
read_file_htseq_fpkm <- fread(file = "/Volumes/G-DRIVE/data/gdc/tcga/rnaseq/read/metadata/gdc_sample_sheet.tumours.htseq_fpkm.tsv", check.names = T) %>%
  as_tibble()

read_htseq_fpkm <- merge_files(read_file_htseq_fpkm, read_dir, gene_list)

read_htseq_fpkm <- read_htseq_fpkm %>%
  inner_join(protein_coding[, c("gene_id", "gene_name")], by =  c("gene" = "gene_id")) %>%
  select(-gene) %>%
  select(gene = gene_name, everything())

#read samples
read_samples <- tibble(sample = colnames(read_htseq_fpkm)[-c(1)], batch = "tcga-coread", cancer = "coread")




# OV

ov_dir <- "/Volumes/G-DRIVE/data/gdc/tcga/rnaseq/ov/files/"

# HTSeq counts
ov_file_htseq_counts <- fread(file = "/Volumes/G-DRIVE/data/gdc/tcga/rnaseq/ov/metadata/gdc_sample_sheet.tumours.htseq_counts.tsv", check.names = T) %>%
  as_tibble()

ov_htseq_counts <- merge_files(ov_file_htseq_counts, ov_dir, gene_list)

ov_htseq_counts <- ov_htseq_counts %>%
  inner_join(protein_coding[, c("gene_id", "gene_name")], by =  c("gene" = "gene_id")) %>%
  select(-gene) %>%
  select(gene = gene_name, everything())


# HTSeq FPKM
ov_file_htseq_fpkm <- fread(file = "/Volumes/G-DRIVE/data/gdc/tcga/rnaseq/ov/metadata/gdc_sample_sheet.tumours.htseq_fpkm.tsv", check.names = T) %>%
  as_tibble()

ov_htseq_fpkm <- merge_files(ov_file_htseq_fpkm, ov_dir, gene_list)

ov_htseq_fpkm <- ov_htseq_fpkm %>%
  inner_join(protein_coding[, c("gene_id", "gene_name")], by =  c("gene" = "gene_id")) %>%
  select(-gene) %>%
  select(gene = gene_name, everything())

#ov samples
ov_samples <- tibble(sample = colnames(ov_htseq_fpkm)[-c(1)], batch = "tcga-hgsc", cancer = "hgsc")



# join sample batch and cancer together
tcga_rna_samples <- bind_rows(
  brca_samples,
  coad_samples,
  read_samples,
  ov_samples
)

write.table(tcga_rna_samples, "./output/files/transcriptomics_samples_tcga.txt", sep="\t", quote=F, row.names=F)



# join all datasets together by gene
# keep only the genes in common to the three datasets
# HTSeq FPKM

tcga_fpkm <- inner_join(brca_htseq_fpkm, coad_htseq_fpkm, by = c("gene")) %>%
  inner_join(read_htseq_fpkm, by = c("gene")) %>%
  inner_join(ov_htseq_fpkm, by = c("gene")) %>%
  select_if(.predicate = colnames(.) %in% c("gene", proteomics$sample))

gz1 <- gzfile("./output/files/transcriptomics_tcga_fpkm.txt.gz", "w")
write.table(tcga_fpkm, gz1, sep="\t", quote=F, row.names=F)
close(gz1)



# join all datasets together by gene
# keep only the genes in common to the three datasets
# HTSeq counts

tcga_counts <- inner_join(brca_htseq_counts, coad_htseq_counts, by = c("gene")) %>%
  inner_join(read_htseq_counts, by = c("gene")) %>%
  inner_join(ov_htseq_counts, by = c("gene")) %>%
  select_if(.predicate = str_sub(colnames(.),1,12) %in% c("gene", proteomics$sample))

gz2 <- gzfile("./output/files/transcriptomics_tcga_counts.txt.gz", "w")
write.table(tcga_counts, gz2, sep="\t", quote=F, row.names=F)
close(gz2)

