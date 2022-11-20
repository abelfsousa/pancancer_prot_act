# protein sequences of the canonical uniprot transcripts


library(tidyverse)


# get canonical uniprot transcripts
canonical_trpt <- read_tsv(file = "./data/uniprot/canonical_trpt/isoform_overrides_uniprot.txt", skip = 1)


# protein sequences
prot_seqs <- read_tsv(file = "./output/files/protein_sequences_GRCh37.txt") %>%
  filter(state == "known" & gene_biotype == "protein_coding" & trpt_biotype == "protein_coding")


# protein sequences for canonical uniprot transcripts
gene_prot_seq <- canonical_trpt %>%
  select(gene_name, trpt_id=enst_id) %>%
  inner_join(prot_seqs[, c("trpt_id", "protein_id", "seq")], by = "trpt_id")

write.table(gene_prot_seq, file = "./output/files/uniprot_canTrpt_protSeqGRCh37.txt", sep = "\t", row.names = F, quote = FALSE)
