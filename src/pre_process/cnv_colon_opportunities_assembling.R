library(tidyverse)


# assembling of cnv data from "Colon Cancer Therapeutic Opportunities"
# GISTIC2 scores from linkedomics



# cancer samples with proteomics measurements
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt") %>%
  filter(str_detect(batch, "colon-opportunities"))


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


colon <- read_tsv(file = "./data/dna/cnv/linkedomics/colon/Human__CPTAC_COAD__VU__SCNA__ExomeSeq__01_28_2016__BCM__Gene__BCM_CopyWriteR_GISTIC2_threshold.cct") %>%
  rename(gene=attrib_name) %>%
  filter(gene %in% protein_coding$gene_name)


# samples
samples <- tibble(sample = colnames(colon)[-c(1)], batch = "colon-opportunities", cancer = "coread")
write.table(samples, "./output/files/cnv_samples_colon_oppt.txt", sep="\t", quote=F, row.names=F)


colon <- colon %>%
  select_if(.predicate = colnames(.) %in% c("gene", proteomics$sample))

gz1 <- gzfile("./output/files/cnv_colon_oppt.txt.gz", "w")
write.table(colon, gz1, sep="\t", quote=F, row.names=F)
close(gz1)

