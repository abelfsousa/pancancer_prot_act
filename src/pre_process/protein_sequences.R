# protein sequences

# GRCh37
p_sequences <- Biostrings::readAAStringSet(filepath = "./data/protein/protein_sequences/Homo_sapiens.GRCh37.pep.all.fa.gz")
p_sequences <- tibble(ids = names(p_sequences), seq = BiocGenerics::paste(p_sequences)) %>%
  separate(col = ids, into = c("protein_id", "state", "chr", "gene_id", "trpt_id", "gene_biotype", "trpt_biotype"), sep = " ") %>%
  mutate(
    state = str_split_fixed(state, "\\:",2)[,2],
    chr = str_split_fixed(chr, "\\:",5)[,3],
    gene_id = str_split_fixed(gene_id, "\\:",2)[,2],
    trpt_id = str_split_fixed(trpt_id, "\\:",2)[,2],
    gene_biotype = str_split_fixed(gene_biotype, "\\:",2)[,2],
    trpt_biotype = str_split_fixed(trpt_biotype, "\\:",2)[,2]) %>%
  mutate(
    protein_id = str_split_fixed(protein_id, "\\.",2)[,1],
    gene_id = str_split_fixed(gene_id, "\\.",2)[,1],
    trpt_id = str_split_fixed(trpt_id, "\\.",2)[,1])

write.table(p_sequences, file = "./output/files/protein_sequences_GRCh37.txt", sep = "\t", row.names = F, quote = FALSE)


# GRCh38
p_sequences <- Biostrings::readAAStringSet(filepath = "./data/protein/protein_sequences/Homo_sapiens.GRCh38.pep.all.fa.gz")
p_sequences <- tibble(ids = names(p_sequences), seq = BiocGenerics::paste(p_sequences)) %>%
  separate(col = ids, into = c("protein_id", "state", "chr", "gene_id", "trpt_id", "gene_biotype", "trpt_biotype", "gene_symbol"), sep = " ", extra = "drop") %>%
  select(-state) %>%
  mutate(
    chr = str_split_fixed(chr, "\\:",5)[,3],
    gene_id = str_split_fixed(gene_id, "\\:",2)[,2],
    trpt_id = str_split_fixed(trpt_id, "\\:",2)[,2],
    gene_biotype = str_split_fixed(gene_biotype, "\\:",2)[,2],
    trpt_biotype = str_split_fixed(trpt_biotype, "\\:",2)[,2],
    gene_symbol = str_split_fixed(gene_symbol, "\\:",2)[,2]) %>%
  mutate(
    protein_id = str_split_fixed(protein_id, "\\.",2)[,1],
    gene_id = str_split_fixed(gene_id, "\\.",2)[,1],
    trpt_id = str_split_fixed(trpt_id, "\\.",2)[,1])

write.table(p_sequences, file = "./output/files/protein_sequences_GRCh38.txt", sep = "\t", row.names = F, quote = FALSE)
