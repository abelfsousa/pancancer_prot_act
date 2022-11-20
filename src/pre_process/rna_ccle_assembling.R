library(tidyverse)

source("src/utils/remove_dup_genes.R")

# assembling of CCLE rna-seq data
# FPKM RSEM expression values


# cancer samples with proteomics measurements
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt") %>%
  filter(str_detect(batch, "cell-lines"))


# get protein coding genes
protein_coding <- read_tsv("./output/files/gencode.v19.annotation.txt") %>%
  select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  select(-gene_type) %>%
  group_by(gene_name) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)


# ccle cell lines metadata
cell_lines_meta <- read_tsv(file = "/Volumes/G-DRIVE/data/broad/ccle/Cell_lines_annotations_20181226.txt") %>%
  select(CCLE_ID, tissue=type_refined, cancer=tcga_code) %>%
  mutate(cancer = tolower(cancer)) %>%
  mutate(cancer = str_replace(cancer, "coad/read", "coread")) %>%
  filter(cancer %in% c("coread", "brca") & tissue %in% c("colorectal", "breast")) %>%
  separate(col = "CCLE_ID", into = c("sample", "tissue2"), sep = "_", extra = "merge") %>%
  select(-tissue2) %>%
  mutate(sample = str_replace(sample, "COLO320", "COLO320HSR")) %>%
  mutate(sample = str_replace(sample, "MDAMB175VII", "MDAMB175")) %>%
  mutate(sample = str_replace(sample, "MDAMB134VI", "MDAMB134"))

#cell_lines %>% filter(!(sample %in% cell_lines_meta$cell))
#cell_lines_meta  %>% filter(str_detect(cell, "COLO320"))
#cell_lines_meta  %>% filter(str_detect(cell, "MDAMB"))


# expression matrix
ccle_rsem_fpkm <- data.table::fread(file = "/Volumes/G-DRIVE/data/broad/ccle/CCLE_RNAseq_genes_rpkm_20180929.gct.gz") %>%
  as_tibble() %>%
  select(-TT_OESOPHAGUS, -TT_THYROID, -KMH2_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE) %>%
  rename_at(.vars = 3:length(.), .funs = ~ str_split_fixed(.x, "_", n=2)[,1]) %>%
  rename(gene = Name, gene_name = Description) %>%
  rename(COLO320HSR = COLO320, MDAMB175 = MDAMB175VII, MDAMB134 = MDAMB134VI) %>%
  mutate(gene = str_split_fixed(gene, "\\.", 2)[,1]) %>%
  select_if(.predicate = colnames(.) %in% c("gene", "gene_name", cell_lines_meta$sample)) %>%
  semi_join(protein_coding, by = c("gene" = "gene_id", "gene_name")) %>%
  select(-gene) %>%
  rename(gene = gene_name)


# ccle brca/coread cell lines with mRNA data
ccle_samples <- tibble(sample = colnames(ccle_rsem_fpkm)[-c(1)], batch = "ccle") %>%
  inner_join(cell_lines_meta[, c("sample", "cancer")], by = c("sample"))

write.table(ccle_samples, "./output/files/transcriptomics_samples_ccle.txt", sep="\t", quote=F, row.names=F)

  

# ccle brca/coread cell lines with mRNA and protein data
ccle_rsem_fpkm <- ccle_rsem_fpkm %>%
  select_if(.predicate = colnames(.) %in% c("gene", proteomics$sample))

gz1 <- gzfile("./output/files/transcriptomics_ccle_fpkm.txt.gz", "w")
write.table(ccle_rsem_fpkm, gz1, sep="\t", quote=F, row.names=F)
close(gz1)

