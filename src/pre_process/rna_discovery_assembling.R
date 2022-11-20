library(tidyverse)
library(data.table)

source("src/utils/merge_files.R")
source("src/utils/remove_dup_genes.R")



# assembling of rna-seq data from discovery cptac studies


# cancer samples with proteomics measurements
protein_samples <- read_tsv(file = "./output/files/proteomics_samples.txt")


# gene list
gene_list <- read_tsv("/Volumes/G-DRIVE/data/gdc/CPTAC-3/rnaseq/gene_annotation/gene_ids.txt")


# gene annotation
gene_annotation <- read_tsv("/Volumes/G-DRIVE/data/gdc/CPTAC-3/rnaseq/gene_annotation/geneAnnot.gencode.v22.txt")


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




# CCRCC

ccrcc_dir <- "/Volumes/G-DRIVE/data/gdc/CPTAC-3/rnaseq/ccrcc/files/"

# HTSeq counts
ccrcc_file_htseq_counts <- fread(file = "/Volumes/G-DRIVE/data/gdc/CPTAC-3/rnaseq/ccrcc/metadata/gdc_sample_sheet.tumours.htseq_counts.tsv", check.names = T) %>%
  as_tibble() %>%
  separate(col = Case.ID, into = c("Case.ID1", "Case.ID2"), sep = ",") %>%
  select(-Case.ID2, Case.ID = Case.ID1) %>%
  # C3L-00908 has two separate files corresponding to two diferent tumors. The protein abundance was calculated with C3L-00908-03 and not C3L-00908-01
  filter(Sample.ID != "C3L-00908-01")

ccrcc_htseq_counts <- merge_files(ccrcc_file_htseq_counts, ccrcc_dir, gene_list)

ccrcc_htseq_counts <- ccrcc_htseq_counts %>%
  inner_join(protein_coding[, c("gene_id", "gene_name")], by =  c("gene" = "gene_id")) %>%
  select(-gene) %>%
  select(gene = gene_name, everything())


# HTSeq FPKM
ccrcc_file_htseq_fpkm <- fread(file = "/Volumes/G-DRIVE/data/gdc/CPTAC-3/rnaseq/ccrcc/metadata/gdc_sample_sheet.tumours.htseq_fpkm.tsv", check.names = T) %>%
  as_tibble() %>%
  separate(col = Case.ID, into = c("Case.ID1", "Case.ID2"), sep = ",") %>%
  select(-Case.ID2, Case.ID = Case.ID1) %>%
  # C3L-00908 has two separate files corresponding to two diferent tumors. The protein abundance was calculated with C3L-00908-03 and not C3L-00908-01
  filter(Sample.ID != "C3L-00908-01")

ccrcc_htseq_fpkm <- merge_files(ccrcc_file_htseq_fpkm, ccrcc_dir, gene_list)

ccrcc_htseq_fpkm <- ccrcc_htseq_fpkm %>%
  inner_join(protein_coding[, c("gene_id", "gene_name")], by =  c("gene" = "gene_id")) %>%
  select(-gene) %>%
  select(gene = gene_name, everything())

#ccrcc samples
ccrcc_samples <- tibble(sample = colnames(ccrcc_htseq_fpkm)[-c(1)], batch = "discovery-ccrcc", cancer = "ccrcc")




# LUAD

luad_dir <- "/Volumes/G-DRIVE/data/gdc/CPTAC-3/rnaseq/luad/files/"

# HTSeq counts
luad_file_htseq_counts <- fread(file = "/Volumes/G-DRIVE/data/gdc/CPTAC-3/rnaseq/luad/metadata/gdc_sample_sheet.tumours.htseq_counts.tsv", check.names = T) %>%
  as_tibble() %>%
  separate(col = Case.ID, into = c("Case.ID1", "Case.ID2"), sep = ",") %>%
  select(-Case.ID2, Case.ID = Case.ID1)

luad_htseq_counts <- merge_files(luad_file_htseq_counts, luad_dir, gene_list)

#keep duplicated samples (coming from multiple files for the same sample)
luad_htseq_counts <- luad_htseq_counts %>%
  inner_join(protein_coding[, c("gene_id", "gene_name")], by =  c("gene" = "gene_id")) %>%
  select(-gene) %>%
  select(gene = gene_name, everything())


# HTSeq FPKM
luad_file_htseq_fpkm <- fread(file = "/Volumes/G-DRIVE/data/gdc/CPTAC-3/rnaseq/luad/metadata/gdc_sample_sheet.tumours.htseq_fpkm.tsv", check.names = T) %>%
  as_tibble() %>%
  separate(col = Case.ID, into = c("Case.ID1", "Case.ID2"), sep = ",") %>%
  select(-Case.ID2, Case.ID = Case.ID1)

luad_htseq_fpkm <- merge_files(luad_file_htseq_fpkm, luad_dir, gene_list)

#average duplicated samples by gene (coming from multiple files for the same sample)
luad_htseq_fpkm <- luad_htseq_fpkm %>%
  inner_join(protein_coding[, c("gene_id", "gene_name")], by =  c("gene" = "gene_id")) %>%
  select(-gene) %>%
  select(gene = gene_name, everything()) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  mutate(sample = str_replace(sample, ".x|.y", "")) %>%
  group_by(sample, gene) %>%
  summarise(fpkm = mean(fpkm)) %>%
  ungroup() %>%
  spread(key = "sample", value = "fpkm")

#luad samples
luad_samples <- tibble(sample = colnames(luad_htseq_fpkm)[-c(1)], batch = "discovery-luad", cancer = "luad")




#UCEC

ucec_dir <- "/Volumes/G-DRIVE/data/gdc/CPTAC-3/rnaseq/ucec/files/"

# HTSeq counts
ucec_file_htseq_counts <- fread(file = "/Volumes/G-DRIVE/data/gdc/CPTAC-3/rnaseq/ucec/metadata/gdc_sample_sheet.tumours.htseq_counts.tsv", check.names = T) %>%
  as_tibble() %>%
  separate(col = Case.ID, into = c("Case.ID1", "Case.ID2"), sep = ",") %>%
  select(-Case.ID2, Case.ID = Case.ID1) %>%
  # C3N-01825 has two separate files corresponding to two diferent tumors. The protein abundance was calculated with C3L-00908-01 and not C3L-00908-03
  filter(Sample.ID != "C3N-01825-03")

ucec_htseq_counts <- merge_files(ucec_file_htseq_counts, ucec_dir, gene_list)

ucec_htseq_counts <- ucec_htseq_counts %>%
  inner_join(protein_coding[, c("gene_id", "gene_name")], by =  c("gene" = "gene_id")) %>%
  select(-gene) %>%
  select(gene = gene_name, everything())


# HTSeq FPKM
ucec_file_htseq_fpkm <- fread(file = "/Volumes/G-DRIVE/data/gdc/CPTAC-3/rnaseq/ucec/metadata/gdc_sample_sheet.tumours.htseq_fpkm.tsv", check.names = T) %>%
  as_tibble() %>%
  separate(col = Case.ID, into = c("Case.ID1", "Case.ID2"), sep = ",") %>%
  select(-Case.ID2, Case.ID = Case.ID1) %>%
  # C3N-01825 has two separate files corresponding to two diferent tumors. The protein abundance was calculated with C3L-00908-01 and not C3L-00908-03
  filter(Sample.ID != "C3N-01825-03")

ucec_htseq_fpkm <- merge_files(ucec_file_htseq_fpkm, ucec_dir, gene_list)

ucec_htseq_fpkm <- ucec_htseq_fpkm %>%
  inner_join(protein_coding[, c("gene_id", "gene_name")], by =  c("gene" = "gene_id")) %>%
  select(-gene) %>%
  select(gene = gene_name, everything())

#ucec samples
ucec_samples <- tibble(sample = colnames(ucec_htseq_fpkm)[-c(1)], batch = "discovery-ucec", cancer = "ucec")



# join sample batch and cancer together
discovery_rna_samples <- bind_rows(
  ccrcc_samples,
  luad_samples,
  ucec_samples
)

write.table(discovery_rna_samples, "./output/files/transcriptomics_samples_discovery.txt", sep="\t", quote=F, row.names=F)



# join all datasets together by gene
# keep only the genes in common to the three datasets
# HTSeq FPKM

discovery_fpkm <- inner_join(ccrcc_htseq_fpkm, luad_htseq_fpkm, by = c("gene")) %>%
  inner_join(ucec_htseq_fpkm, by = c("gene")) %>%
  select_if(.predicate = colnames(.) %in% c("gene", protein_samples$sample))

gz1 <- gzfile("./output/files/transcriptomics_discovery_fpkm.txt.gz", "w")
write.table(discovery_fpkm, gz1, sep="\t", quote=F, row.names=F)
close(gz1)



# join all datasets together by gene
# keep only the genes in common to the three datasets
# HTSeq counts

discovery_counts <- inner_join(ccrcc_htseq_counts, luad_htseq_counts, by = c("gene")) %>%
  inner_join(ucec_htseq_counts, by = c("gene")) %>%
  select_if(.predicate = str_split_fixed(colnames(.), ".x|.y", 2)[,1] %in% c("gene", protein_samples$sample))

gz2 <- gzfile("./output/files/transcriptomics_discovery_counts.txt.gz", "w")
write.table(discovery_counts, gz2, sep="\t", quote=F, row.names=F)
close(gz2)

