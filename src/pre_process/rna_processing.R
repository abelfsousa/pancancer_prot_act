#' ---
#' title: "Transcriptomics data"
#' author: "Abel Sousa"
#' date: "October 1st, 2019"
#' ---

#+ setup, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)


#' load R packages
suppressMessages(library(limma))
suppressMessages(library(viridis))
suppressMessages(library(tidyverse))


#' cancer samples with **proteomics** measurements
protein_samples <- read_tsv(file = "./output/files/proteomics_samples.txt") %>%
  mutate(batch = str_replace(batch, "cell-lines-law|cell-lines-lpk|cell-lines-rmlt", "ccle"))
nrow(protein_samples)


#' cancer samples with **phosphoproteomics** measurements
phospho_samples <- read_tsv(file = "./output/files/phosphoproteomics_samples.txt") %>%
  mutate(batch = str_replace(batch, "cell-lines-rmlt", "ccle"))
nrow(phospho_samples)


#' cancer samples with **transcriptomics** measurements
rna_tcga <- read_tsv("./output/files/transcriptomics_samples_tcga.txt")
rna_discovery <- read_tsv("./output/files/transcriptomics_samples_discovery.txt")
rna_colon <- read_tsv("./output/files/transcriptomics_colon_oppt_samples.txt")
rna_eogc <- read_tsv("./output/files/transcriptomics_samples_eogc.txt")
rna_ccle <- read_tsv("./output/files/transcriptomics_samples_ccle.txt")
rna_cbttc <- read_tsv("./output/files/transcriptomics_samples_cbttc.txt")
rna_hcc <- read_tsv("./output/files/transcriptomics_hcc_samples.txt")

rna_samples <- bind_rows(rna_tcga, rna_discovery, rna_colon, rna_eogc, rna_ccle, rna_cbttc, rna_hcc)
nrow(rna_samples)

write.table(rna_samples, "./output/files/transcriptomics_samples.txt", sep="\t", quote=F, row.names=F)

#rna_samples <- rna_samples %>%
#  mutate(batch = str_c("cptac", batch, sep = "-")) %>%
#  mutate(batch = str_replace(batch, "cptac-ccle", "ccle"))


#' cancer samples **with protein and mRNA**
multi_samples1 <- protein_samples %>%
  filter(sample %in% rna_samples$sample)
nrow(multi_samples1)

write.table(multi_samples1, "./output/files/multi_samples_prot_rna.txt", sep="\t", quote=F, row.names=F)


#' cancer samples **with protein, phosphorylation and mRNA**
multi_samples2 <- multi_samples1 %>%
  filter(sample %in% phospho_samples$sample)
nrow(multi_samples2)

write.table(multi_samples2, "./output/files/multi_samples_prot_phos_rna.txt", sep="\t", quote=F, row.names=F)


#' ## number of samples in each experimental batch
#' ### cancer samples with mRNA and protein abundance
samples_barplot <- multi_samples1 %>%
  mutate(batch = fct_rev(fct_infreq(batch))) %>%
  ggplot() +
  geom_bar(mapping = aes(x = batch, fill = cancer), stat = "count") +
  coord_flip() +
  theme_classic() +
  scale_fill_viridis(discrete = T) +
  theme(
    axis.title = element_text(color = "black", size = 16),
    legend.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 14)
  ) +
  scale_y_continuous(limits = c(NA, 200)) +
  labs(x = "Batch", y = "Samples", fill = "Cancer")

#+ fig.width=10, fig.height=5
samples_barplot

ggsave(filename="rna_samples.png", plot = samples_barplot, path = "./output/plots/pre_process/", width = 10, height = 5)
unlink("rna_samples.png")




#' load transcriptomics data
#' FPKM values
pattern <- c("Atypical Teratoid Rhabdoid Tumor \\(ATRT\\)|Craniopharyngioma|Ependymoma|Ganglioglioma|High-grade glioma/astrocytoma|Low-grade glioma/astrocytoma|Medulloblastoma")
rna <- data.table::fread(input = "./output/files/transcriptomics_tcga_fpkm.txt.gz") %>%
  as_tibble() %>%
  inner_join(data.table::fread(input = "./output/files/transcriptomics_discovery_fpkm.txt.gz"), by = "gene") %>%
  inner_join(data.table::fread(input = "./output/files/transcriptomics_colon_oppt_fpkm_rsem.txt.gz"), by = "gene") %>%
  inner_join(data.table::fread(input = "./output/files/transcriptomics_eogc_fpkm.txt.gz"), by = "gene") %>%
  inner_join(data.table::fread(input = "./output/files/transcriptomics_ccle_fpkm.txt.gz"), by = "gene") %>%
  inner_join(data.table::fread(input = "./output/files/transcriptomics_cbttc_fpkm.txt.gz"), by = "gene") %>%
  inner_join(data.table::fread(input = "./output/files/transcriptomics_hcc_fpkm_rsem.txt.gz"), by = "gene") %>%
  gather(key = "sample", value = "fpkm", -gene) %>%
  inner_join(multi_samples1, by = "sample") %>%
  mutate(cancer2 = str_replace(cancer, pattern, "cbttc")) %>%
  select(gene, sample, cancer, cancer2, batch, fpkm) %>%
  group_by(gene, batch, cancer2) %>%
  filter((sum(fpkm > 0)) > (n()*0.5)) %>%
  ungroup() %>%
  mutate(log2fpkm = log2(fpkm + 1)) %>%
  group_by(gene, batch, cancer2) %>%
  mutate(log2fc = log2(fpkm+1)-log2(median(fpkm+1))) %>%
  ungroup() %>%
  select(-cancer2)




#' ## log2FPKM distribution by sample and experimental batch
rna_log2fpkm_distribution <- ggplot(data = rna, mapping = aes(x = sample, y = log2fpkm, fill = cancer)) +
  theme_classic() +
  geom_boxplot(outlier.size = 0.1, lwd = 0.1) +
  facet_wrap(~ batch, ncol = 3, scales = "free") +
  theme(
    axis.title=element_text(colour="black", size=22),
    axis.text.y=element_text(colour="black", size=18),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.background=element_blank(),
    strip.text.x=element_text(colour="black", size=22),
    legend.position = "bottom",
    legend.title = element_text(colour="black", size=22),
    legend.text = element_text(colour="black", size=18)) +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "Sample", y = "log2 FPKM") +
  guides(fill = guide_legend(nrow = 2))

#+ fig.width=30, fig.height=20
rna_log2fpkm_distribution

ggsave(filename="rna_log2fpkm_distribution_1.png", plot = rna_log2fpkm_distribution, path = "./output/plots/pre_process/", width = 30, height = 20)
unlink("rna_log2fpkm_distribution_1.png")


gz1 <- gzfile("./output/files/transcriptomics_log2fpkm.txt.gz", "w")
write.table(spread(rna[, c("gene", "sample", "log2fpkm")], key = "sample", value = "log2fpkm"), gz1, sep="\t", quote=F, row.names=F)
close(gz1)



#' ## log2FC distribution by sample and experimental batch
rna_log2fc_distribution <- ggplot(data = rna, mapping = aes(x = sample, y = log2fc, fill = cancer)) +
  theme_classic() +
  geom_boxplot(outlier.size = 0.1, lwd = 0.1) +
  facet_wrap(~ batch, ncol = 3, scales = "free") +
  theme(
    axis.title=element_text(colour="black", size=22),
    axis.text.y=element_text(colour="black", size=18),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.background=element_blank(),
    strip.text.x=element_text(colour="black", size=22),
    legend.position = "bottom",
    legend.title = element_text(colour="black", size=22),
    legend.text = element_text(colour="black", size=18)) +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "Sample", y = "log2 FC") +
  guides(fill = guide_legend(nrow = 2))

#+ fig.width=30, fig.height=20
rna_log2fc_distribution

ggsave(filename="rna_log2fc_distribution_1.png", plot = rna_log2fc_distribution, path = "./output/plots/pre_process/", width = 30, height = 20)
unlink("rna_log2fc_distribution_1.png")


gz2 <- gzfile("./output/files/transcriptomics_log2fc.txt.gz", "w")
write.table(spread(rna[, c("gene", "sample", "log2fc")], key = "sample", value = "log2fc"), gz2, sep="\t", quote=F, row.names=F)
close(gz2)



