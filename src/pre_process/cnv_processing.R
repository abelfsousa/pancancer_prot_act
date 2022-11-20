library(tidyverse)



# cnv data processing


# cancer samples with proteomics measurements
protein_samples <- read_tsv(file = "./output/files/proteomics_samples.txt") %>%
  mutate(batch = str_replace(batch, "cell-lines-law|cell-lines-lpk|cell-lines-rmlt", "ccle"))


# cancer samples with phosphoproteomics measurements
phospho_samples <- read_tsv(file = "./output/files/phosphoproteomics_samples.txt") %>%
  mutate(batch = str_replace(batch, "cell-lines-rmlt", "ccle"))


# cancer samples with mRNA measurements
rna_samples <- read_tsv(file = "./output/files/transcriptomics_samples.txt")


# cancer samples with protein and mRNA
multi_samples_prt_rna <- read_tsv(file = "./output/files/multi_samples_prot_rna.txt")


# cancer samples with protein, phospho and mRNA
multi_samples_prt_pho_rna <- read_tsv(file = "./output/files/multi_samples_prot_phos_rna.txt")


# all samples with cnv data
all_samples <- bind_rows(
  read_tsv(file = "./output/files/cnv_samples_tcga.txt"),
  read_tsv(file = "./output/files/cnv_samples_cbttc.txt"),
  read_tsv(file = "./output/files/cnv_samples_colon_oppt.txt"),
  read_tsv(file = "./output/files/cnv_samples_ccle.txt"),
  read_tsv(file = "./output/files/cnv_samples_hcc.txt"),
  read_tsv(file = "./output/files/cnv_samples_ccrcc.txt"),
  read_tsv(file = "./output/files/cnv_samples_ucec.txt"),
)

write.table(all_samples, "./output/files/cnv_samples.txt", sep="\t", quote=F, row.names=F)


# all cnv from all datasets (with protein data)
all_cnv <- read_tsv(file = "./output/files/cnv_tcga.txt.gz") %>%
  inner_join(read_tsv(file = "./output/files/cnv_ccle.txt.gz"), by = "gene") %>%
  inner_join(read_tsv(file = "./output/files/cnv_hcc.txt.gz"), by = "gene") %>%
  inner_join(read_tsv(file = "./output/files/cnv_colon_oppt.txt.gz"), by = "gene") %>%
  inner_join(read_tsv(file = "./output/files/cnv_cbttc.txt.gz"), by = "gene") %>%
  inner_join(read_tsv(file = "./output/files/cnv_ccrcc.txt.gz"), by = "gene") %>%
  inner_join(read_tsv(file = "./output/files/cnv_ucec.txt.gz"), by = "gene")

gz1 <- gzfile("./output/files/cnv.txt.gz", "w")
write.table(all_cnv, gz1, sep="\t", quote=F, row.names=F)
close(gz1)

