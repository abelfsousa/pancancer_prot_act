library(tidyverse)


# mutation data processing


# cancer samples with proteomics measurements
protein_samples <- read_tsv(file = "./output/files/proteomics_samples.txt") %>%
  mutate(batch = str_replace(batch, "cell-lines-law|cell-lines-lpk|cell-lines-rmlt", "ccle"))


# cancer samples with phosphoproteomics measurements
phospho_samples <- read_tsv(file = "./output/files/phosphoproteomics_samples.txt") %>%
  mutate(batch = str_replace(batch, "cell-lines-rmlt", "ccle"))


# cancer samples with protein and mRNA
multi_samples_prt_rna <- read_tsv(file = "./output/files/multi_samples_prot_rna.txt")


# cancer samples with protein, phospho and mRNA
multi_samples_prt_pho_rna <- read_tsv(file = "./output/files/multi_samples_prot_phos_rna.txt")



# all samples with mutation data
all_samples <- c("mutations_samples_tcga.txt", "mutations_samples_cbttc.txt", "mutations_samples_discovery.txt", "mutations_colon_oppt_samples.txt",
                 "mutations_samples_ccle.txt", "mutations_eogc_samples.txt", "mutations_hcc_samples.txt")

all_samples <- map_dfr(.x=all_samples, .f = ~ read_tsv(paste0("./output/files/", .x)))

write.table(all_samples, "./output/files/mutations_samples.txt", sep="\t", quote=F, row.names=F)



# all mutations from all datasets (with protein data)
all_mutations <- c("mutations_tcga.txt.gz", "mutations_cbttc.txt.gz", "mutations_discovery.txt.gz", "mutations_colon_oppt.txt.gz",
                 "mutations_ccle.txt.gz", "mutations_eogc.txt.gz", "mutations_hcc.txt.gz")

all_mutations <- map_dfr(.x=all_mutations, .f = ~ read_tsv(paste0("./output/files/", .x))) %>%
  distinct() %>%
  select(-consequence, -biotype)

gz1 <- gzfile("./output/files/mutations.txt.gz", "w")
write.table(all_mutations, gz1, sep="\t", quote=F, row.names=F)
close(gz1)



# gene_symbols, transcript IDs, protein IDs and sequences
gene_seq <- all_mutations %>%
  select(gene_symbol, transcript_id, protein_id, seq) %>%
  distinct()
write.table(gene_seq, "./output/files/all_mutations_gene_trpt_prot_seq.txt", sep="\t", quote=F, row.names=F)
