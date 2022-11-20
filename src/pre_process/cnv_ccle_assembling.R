library(tidyverse)


# assembling of cnv data from ccle cell lines
# GISTIC2 scores from cbioportal



# cancer samples with proteomics measurements
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt") %>%
  filter(str_detect(batch, "cell-lines"))


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



# CNV data
ccle <- read_tsv(file = "./data/dna/cnv/cbioportal/cellline_ccle_broad/data_CNA.txt") %>%
  rename(gene=Hugo_Symbol) %>%
  select(-Entrez_Gene_Id) %>%
  filter(gene %in% protein_coding$gene_name) %>%
  group_by(gene) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n) %>%
  select(-c(TT_OESOPHAGUS, TT_THYROID)) %>%
  rename_at(.vars = 2:ncol(.), .funs = ~ str_split(.x, "_", 2, simplify = T)[,1]) %>%
  rename(COLO320HSR = COLO320, MDAMB175 = MDAMB175VII, MDAMB134 = MDAMB134VI) %>%
  select_if(.predicate = colnames(.) %in% c("gene", cell_lines_meta$sample))



# samples
ccle_samples <- tibble(sample = colnames(ccle)[-c(1)], batch = "ccle") %>%
  inner_join(cell_lines_meta[, c("sample", "cancer")], by = "sample")

write.table(ccle_samples, "./output/files/cnv_samples_ccle.txt", sep="\t", quote=F, row.names=F)


ccle <- ccle %>%
  select_if(.predicate = colnames(.) %in% c("gene", proteomics$sample))

gz1 <- gzfile("./output/files/cnv_ccle.txt.gz", "w")
write.table(ccle, gz1, sep="\t", quote=F, row.names=F)
close(gz1)
