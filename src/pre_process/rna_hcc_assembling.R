library(tidyverse)


# assembling of HCC rna-seq data
# FPKM RSEM expression values


# cancer samples with proteomics measurements
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt") %>%
  filter(cancer == "hcc")


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


rna_hcc_counts <- read_tsv(file = "./data/rna/hcc/HCC_UQ_RSEM.tsv.gz") %>%
  select(-X344) %>%
  rename(gene = protein) %>%
  mutate(gene = str_split_fixed(gene, "\\.", 2)[,1]) %>%
  inner_join(protein_coding, by = c("gene" = "gene_id")) %>%
  select(gene, gene_name, everything()) %>%
  rename_at(.vars = c(3:ncol(.)), .funs = ~ str_c("T", .x, sep = ""))



#http://www.genemine.org/gtftools.php
#https://www.biorxiv.org/content/10.1101/263517v1
gene_lengths <- read_tsv("./output/files/gencode.v19.gene_lengths.txt") %>%
  mutate(gene = str_split_fixed(gene, "\\.", 2)[,1]) %>%
  select(gene, length = merged)



# calculate FPKM
rna_hcc_fpkm <- rna_hcc_counts %>%
  inner_join(gene_lengths, by = c("gene")) %>%
  select(-gene) %>%
  select(gene = gene_name, length, everything())

rna_hcc_fpkm <- edgeR::rpkm(y = column_to_rownames(rna_hcc_fpkm[, -c(2)], var = "gene"), gene.length = rna_hcc_fpkm$length, log = F) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()



# samples
samples <- tibble(sample = colnames(rna_hcc_fpkm)[-c(1)], batch = "hcc", cancer = "hcc")
write.table(samples, "./output/files/transcriptomics_hcc_samples.txt", sep="\t", quote=F, row.names=F)



# plot <- rna_hcc_fpkm %>%
#   gather(-gene, key = "sample", value = "fpkm") %>%
#   mutate(log2fpkm = log2(fpkm + 1)) %>%
#   ggplot(mapping = aes(x = sample, y = log2fpkm)) +
#   theme_classic() +
#   geom_boxplot(outlier.size = 0.1, lwd = 0.1) +
#   theme(
#     axis.title=element_text(colour="black", size=22),
#     axis.text.y=element_text(colour="black", size=18),
#     axis.text.x=element_blank(),
#     axis.ticks.x=element_blank(),
#     strip.background=element_blank(),
#     strip.text.x=element_text(colour="black", size=22),
#     legend.position = "bottom",
#     legend.title = element_text(colour="black", size=22),
#     legend.text = element_text(colour="black", size=18)) +
#   labs(x = "Sample", y = "log2 FPKM")
# plot


rna_hcc_fpkm <- rna_hcc_fpkm %>%
  select_if(.predicate = colnames(.) %in% c("gene", proteomics$sample))


gz1 <- gzfile("./output/files/transcriptomics_hcc_fpkm_rsem.txt.gz", "w")
write.table(rna_hcc_fpkm, gz1, sep="\t", quote=F, row.names=F)
close(gz1)

