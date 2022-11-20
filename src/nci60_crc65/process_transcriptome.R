library(tidyverse)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# -- RNA-seq data


# load cell lines
nci60_crc65 <- read_tsv("./output/files/nci60_crc65_cell_lines.txt")


# NCI60 cell lines
nci60_rna <- data.table::fread("./data/nci60_crc65/transcriptome/nci60/nci60_rna.txt") %>%
  as_tibble() %>%
  select(-2:-6) %>%
  rename(gene=`Gene name`) %>%
  pivot_longer(-gene, names_to = "sample", values_to = "log2fpkm") %>%
  mutate(sample = str_split_fixed(sample, ":", 2)[,2]) %>%
  mutate(sample = toupper(str_replace_all(sample, "-| ", ""))) %>%
  mutate(fpkm = (2^log2fpkm)-1) %>%
  #mutate(log2fpkm_ = log2(fpkm+1)) %>%
  group_by(gene) %>%
  mutate(log2fc = log2(fpkm+1)-log2(median(fpkm+1))) %>%
  #mutate(log2fc = log2((fpkm+1)/median(fpkm+1))) %>%
  ungroup() %>%
  mutate(batch = "NCI60") %>%
  select(batch, sample, everything())


# CRC65 cell lines

# rna seq data (FPKM) from https://portals.broadinstitute.org/ccle/data
crc65_rna <- data.table::fread(file = "./data/nci60_crc65/transcriptome/crc65/CCLE_RNAseq_genes_rpkm_20180929.gct.gz") %>%
  as_tibble() %>%
  select(-Name, gene = Description) %>%
  group_by(gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  filter(gene %in% unique(nci60_rna$gene)) %>%
  pivot_longer(-gene, names_to = "sample", values_to = "fpkm") %>%
  mutate(sample = str_split_fixed(sample, "_", 2)[,1]) %>%
  filter(sample %in% nci60_crc65[nci60_crc65$batch == "CRC65", "cell_line", drop=T]) %>%
  #group_by(gene, sample) %>%
  #summarise(fpkm = mean(fpkm)) %>%
  #ungroup() %>%
  mutate(log2fpkm = log2(fpkm+1)) %>%
  group_by(gene) %>%
  mutate(log2fc = log2(fpkm+1)-log2(median(fpkm+1))) %>%
  #mutate(log2fc = log2((fpkm+1)/median(fpkm+1))) %>%
  ungroup() %>%
  mutate(batch = "CRC65") %>%
  select(batch, sample, everything())


# rna seq data (log2TPM) from https://depmap.org/portal/download/
crc65_rna2_info <- data.table::fread(file = "./data/nci60_crc65/transcriptome/crc65/sample_info.csv") %>%
  as_tibble()

crc65_rna2 <- data.table::fread(file = "./data/nci60_crc65/transcriptome/crc65/CCLE_expression.csv.gz") %>%
  as_tibble() %>%
  inner_join(crc65_rna2_info[,1:2], by = c("V1" = "DepMap_ID")) %>%
  select(stripped_cell_line_name, everything()) %>%
  select(-V1, sample = stripped_cell_line_name) %>%
  pivot_longer(-sample, names_to = "gene", values_to = "log2tpm") %>%
  mutate(gene = str_split_fixed(gene, " ", 2)[,1]) %>%
  filter(sample %in% nci60_crc65[nci60_crc65$batch == "CRC65", "cell_line", drop=T]) %>%
  filter(gene %in% unique(nci60_rna$gene)) %>%
  mutate(tpm = (2^log2tpm)-1) %>%
  #mutate(log2tpm_ = log2(tpm+1)) %>%
  group_by(gene) %>%
  mutate(log2fc = log2(tpm+1)-log2(median(tpm+1))) %>%
  #mutate(log2fc = log2((tpm+1)/median(tpm+1))) %>%
  ungroup() %>%
  mutate(batch = "CRC65") %>%
  select(batch, sample, everything())


# correlate CRC65 rna seq data using the shared samples and genes between datasets

# across genes
corr_genes <- crc65_rna[, c("sample", "gene", "log2fc")] %>%
  inner_join(crc65_rna2[, c("sample", "gene", "log2fc")], by = c("sample", "gene")) %>%
  group_by(gene) %>%
  summarise(cor = cor(log2fc.x, log2fc.y)) %>%
  ungroup() %>%
  ggplot(mapping = aes(x = cor)) +
  geom_histogram(bins = 10)

# across samples
corr_samples <- crc65_rna[, c("sample", "gene", "log2fc")] %>%
  inner_join(crc65_rna2[, c("sample", "gene", "log2fc")], by = c("sample", "gene")) %>%
  group_by(sample) %>%
  summarise(cor = cor(log2fc.x, log2fc.y)) %>%
  ungroup() %>%
  ggplot(mapping = aes(x = cor)) +
  geom_histogram(bins = 10)


# in NCI60 select genes common to CRC65
nci60_rna <- nci60_rna %>%
  filter(gene %in% unique(crc65_rna2$gene))


rna <- bind_rows(nci60_rna[, c("batch", "sample", "gene", "log2fc")], crc65_rna2[, c("batch", "sample", "gene", "log2fc")])
write_tsv(rna, "./output/files/nci60_crc65_rna_log2fc.txt.gz")

rna_boxplot <- rna %>%
  mutate(batch = fct_relevel(batch, "NCI60", "CRC65")) %>%
  ggplot(mapping = aes(x = sample, y = log2fc, fill = batch)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_wrap(~ batch, ncol = 2, scales = "free_x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    strip.text = element_text(size=16),
    legend.text = element_text(size=14),
    legend.title = element_text(size=16))

ggsave(filename = "nci60_crc65_rna_boxplot.png", plot = rna_boxplot, path = "./output/plots/nci60_crc65/", width = 12, height = 8)
unlink("nci60_crc65_rna_boxplot.png")



# save expression data in matrices format for each cell line group
nci60_log2fc_mat <- nci60_rna %>%
  select(sample, gene, log2fc) %>%
  pivot_wider(names_from = "sample", values_from = "log2fc")

nci60_log2fpkm_mat <- nci60_rna %>%
  select(sample, gene, log2fpkm) %>%
  pivot_wider(names_from = "sample", values_from = "log2fpkm")

write_tsv(nci60_log2fc_mat, "./output/files/nci60_matrix_rna_log2fc.txt")
write_tsv(nci60_log2fpkm_mat, "./output/files/nci60_matrix_rna_log2fpkm.txt")


crc65_log2fc_mat <- crc65_rna2 %>%
  select(sample, gene, log2fc) %>%
  pivot_wider(names_from = "sample", values_from = "log2fc")

crc65_log2tpm_mat <- crc65_rna2 %>%
  select(sample, gene, log2tpm) %>%
  pivot_wider(names_from = "sample", values_from = "log2tpm")

write_tsv(crc65_log2fc_mat, "./output/files/crc65_matrix_rna_log2fc.txt")
write_tsv(crc65_log2tpm_mat, "./output/files/crc65_matrix_rna_log2tpm.txt")
