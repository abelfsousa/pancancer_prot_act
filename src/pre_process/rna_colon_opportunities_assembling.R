library(tidyverse)
library(ggpubr)



# assembling of rna-seq data from "Colon Cancer Therapeutic Opportunities"

# cancer samples with proteomics measurements
protein_samples <- read_tsv(file = "./output/files/proteomics_samples.txt")


# load HTSeq counts
rna_counts <- read_tsv("./data/rna/colon/htseq-gene-counts.txt") %>%
  rename(gene = idx)


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

rna_counts <- rna_counts %>%
  filter(gene %in% protein_coding$gene_name)


# calculate log2CPM
rna_cpm <- edgeR::cpm(column_to_rownames(rna_counts, var = "gene"), log = T, prior.count = 1) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  gather(key = "sample", value = "log2CPM", -gene) %>%
  group_by(gene) %>%
  filter(mean(log2CPM) > 0) %>%
  ungroup() %>%
  spread(key = "sample", value = "log2CPM")




# obtain gene lengths

#human_mart <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
#gene_ids <- biomaRt::getBM(attributes=c("hgnc_symbol", "ensembl_gene_id"), filters="hgnc_symbol", values=rna_counts$gene, mart=human_mart) %>%
#  as_tibble()

gene_ids <- rna_counts %>%
  select(gene_name = gene) %>%
  inner_join(protein_coding, by = "gene_name")

#gene_lengths <- EDASeq::getGeneLengthAndGCContent(gene_ids$gene_id, "hsa", mode=c("biomart")) %>%
#  as.data.frame() %>%
#  rownames_to_column(var = "gene") %>%
#  as_tibble() %>%
#  filter(!is.na(length))

#http://www.genemine.org/gtftools.php
#https://www.biorxiv.org/content/10.1101/263517v1
gene_lengths <- read_tsv("./output/files/gencode.v19.gene_lengths.txt") %>%
  mutate(gene = str_split_fixed(gene, "\\.", 2)[,1]) %>%
  select(gene_id = gene, length = merged)


gene_ids <- gene_ids %>%
  inner_join(gene_lengths, by = "gene_id") %>%
  select(-gene_id) %>%
  rename(gene = gene_name)


# calculate FPKM
rna_fpkm <- rna_counts %>%
  inner_join(gene_ids, by = "gene") %>%
  dplyr::select(gene, length, everything())

rna_fpkm <- edgeR::rpkm(y = column_to_rownames(rna_fpkm[, -c(2)], var = "gene"), gene.length = rna_fpkm$length, log = F) %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble()


# samples
samples <- tibble(sample = colnames(rna_fpkm)[-c(1)], batch = "colon-opportunities", cancer = "coread")
write.table(samples, "./output/files/transcriptomics_colon_oppt_samples.txt", sep="\t", quote=F, row.names=F)


rna_fpkm <- rna_fpkm %>%
  select(c("gene", colnames(.)[(colnames(.) %in% protein_samples$sample)]))

gz1 <- gzfile("./output/files/transcriptomics_colon_oppt_fpkm_htseq.txt.gz", "w")
write.table(rna_fpkm, gz1, sep="\t", quote=F, row.names=F)
close(gz1)




# load RSEM fpkm
rna_fpkm_rsem <- read_tsv("./data/rna/colon/rsem-fpkm-gene.txt") %>%
  dplyr::rename(gene = idx)


# correlate both FPKM values and plot a boxplot
fpkm_boxplot <- rna_fpkm %>%
  gather(key = "sample", value = "fpkm_htseq", -gene) %>%
  inner_join(rna_fpkm_rsem %>% gather(key = "sample", value = "fpkm_rsem", -gene),
             by = c("gene", "sample")) %>%
  group_by(gene) %>%
  summarise(cor = cor(fpkm_htseq, fpkm_rsem)) %>%
  ungroup() %>%
  ggplot(mapping = aes(x = factor(0), y = cor)) +
  geom_boxplot(outlier.shape = NA) +
  theme_classic()
fpkm_boxplot


# calculate FPKM spearman correlation
fpkm_cor <- rna_fpkm %>%
  gather(key = "sample", value = "fpkm_htseq", -gene) %>%
  inner_join(rna_fpkm_rsem %>% gather(key = "sample", value = "fpkm_rsem", -gene),
             by = c("gene", "sample")) %>%
  nest(data = c(sample, fpkm_htseq, fpkm_rsem)) %>%
  #mutate(cor = map(.x = data, .f = function(x) broom::tidy(cor.test(x$fpkm_htseq, x$fpkm_rsem, method = "spearman")))) %>%
  mutate(cor = map(.x = data, .f = ~ broom::tidy(cor.test(.x$fpkm_htseq, .x$fpkm_rsem, method = "spearman")))) %>%
  dplyr::select(-data) %>%
  unnest()
median(fpkm_cor$estimate, na.rm = T)


# plot both FPKM median values for each gene
fpkm_plot <- rna_fpkm %>%
  gather(key = "sample", value = "fpkm_htseq", -gene) %>%
  inner_join(rna_fpkm_rsem %>% gather(key = "sample", value = "fpkm_rsem", -gene),
             by = c("gene", "sample")) %>%
  mutate(fpkm_htseq = log2(fpkm_htseq), fpkm_rsem = log2(fpkm_rsem)) %>%
  group_by(gene) %>%
  summarise(fpkm_htseq = median(fpkm_htseq), fpkm_rsem = median(fpkm_rsem)) %>%
  ungroup() %>%
  ggplot(mapping = aes(x = fpkm_htseq, y = fpkm_rsem)) +
  geom_point() +
  theme_classic() +
  stat_cor(method = "pearson")
fpkm_plot



rna_fpkm_rsem <- rna_fpkm_rsem %>%
  filter(gene %in% protein_coding$gene_name) %>%
  select(c("gene", colnames(.)[(colnames(.) %in% protein_samples$sample)]))

gz2 <- gzfile("./output/files/transcriptomics_colon_oppt_fpkm_rsem.txt.gz", "w")
write.table(rna_fpkm_rsem, gz2, sep="\t", quote=F, row.names=F)
close(gz2)


