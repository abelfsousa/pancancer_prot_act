library(tidyverse)


# assembling of mutation data from tcga studies
# cbioportal datasets


source("./src/utils/mutation_prot_seq_match.R")


# cancer samples with proteomics measurements
protein_samples <- read_tsv(file = "./output/files/proteomics_samples.txt")


# mutations to select
mutations_accepted <- c("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins", "Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site")


# protein sequences
prot_seqs <- read_tsv(file = "./output/files/protein_sequences_GRCh37.txt") %>%
  filter(state == "known" & gene_biotype == "protein_coding" & trpt_biotype == "protein_coding")


# amino acids code
aa_code <- read_tsv(file = "./data/protein/protein_sequences/aminoacids_code.txt", skip = c(1))


# canonical uniprot transcripts
canonical_trpt <- read_tsv(file = "./data/uniprot/canonical_trpt/isoform_overrides_uniprot.txt", skip = 1)
  #filter(!is.na(refseq_id)) %>%
  #mutate(refseq_id = str_split_fixed(refseq_id, "\\.", 2)[,1])


# BRCA
cbio_brca <- data.table::fread(file = "./data/dna/mutation/cbioportal/brca_tcga_pan_can_atlas_2018/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  select(sample = Tumor_Sample_Barcode, consequence=Consequence, variant_class=Variant_Classification, variant_type=Variant_Type, CHROM=Chromosome, POS=Start_Position, REF=Reference_Allele, ALT=Tumor_Seq_Allele2, gene_symbol=Hugo_Symbol, biotype=BIOTYPE, gene_id=Gene, transcript_id=Feature, protein_id=ENSP, HGVSp, t_ref_count, t_alt_count) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(sample = str_sub(sample, 1, 12)) %>%
  mutate(VAF = t_alt_count/(t_alt_count+t_ref_count)) %>%
  select(-t_alt_count, -t_ref_count) %>%
  filter(biotype == "protein_coding") %>%
  filter(variant_class %in% mutations_accepted) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "enst_id"))

cbio_brca2 <- cbio_brca %>%
  filter(sample %in% protein_samples$sample) %>%
  inner_join(prot_seqs[, c("protein_id", "seq")], by = "protein_id") %>%
  group_by(sample) %>%
  mutate(id_var = 1:n()) %>%
  ungroup() %>%
  select(sample, id_var, everything()) %>%
  group_by(sample, id_var) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = cbio_datasets, code = aa_code)) %>%
  unnest(cols = data) %>%
  select(-id_var)

exclude <- cbio_brca2 %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

cbio_brca2 <- cbio_brca2 %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)

# brca samples
brca_samples <- tibble(sample = unique(cbio_brca$sample), batch = "tcga-brca", cancer = "brca")




# COREAD
cbio_coread <- data.table::fread(file = "./data/dna/mutation/cbioportal/coadread_tcga_pan_can_atlas_2018/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  select(sample = Tumor_Sample_Barcode, consequence=Consequence, variant_class=Variant_Classification, variant_type=Variant_Type, CHROM=Chromosome, POS=Start_Position, REF=Reference_Allele, ALT=Tumor_Seq_Allele2, gene_symbol=Hugo_Symbol, biotype=BIOTYPE, gene_id=Gene, transcript_id=Feature, protein_id=ENSP, HGVSp, t_ref_count, t_alt_count) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(sample = str_sub(sample, 1, 12)) %>%
  mutate(VAF = t_alt_count/(t_alt_count+t_ref_count)) %>%
  select(-t_alt_count, -t_ref_count) %>%
  filter(biotype == "protein_coding") %>%
  filter(variant_class %in% mutations_accepted) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "enst_id"))

cbio_coread2 <- cbio_coread %>%
  filter(sample %in% protein_samples$sample) %>%
  inner_join(prot_seqs[, c("protein_id", "seq")], by = "protein_id") %>%
  group_by(sample) %>%
  mutate(id_var = 1:n()) %>%
  ungroup() %>%
  select(sample, id_var, everything()) %>%
  group_by(sample, id_var) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = cbio_datasets, code = aa_code)) %>%
  unnest(cols = data) %>%
  select(-id_var)

exclude <- cbio_coread2 %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

cbio_coread2 <- cbio_coread2 %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)


# coread samples
coread_samples <- tibble(sample = unique(cbio_coread$sample), batch = "tcga-coread", cancer = "coread")




# OV
cbio_ov <- data.table::fread(file = "./data/dna/mutation/cbioportal/ov_tcga_pan_can_atlas_2018/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  select(sample = Tumor_Sample_Barcode, consequence=Consequence, variant_class=Variant_Classification, variant_type=Variant_Type, CHROM=Chromosome, POS=Start_Position, REF=Reference_Allele, ALT=Tumor_Seq_Allele2, gene_symbol=Hugo_Symbol, biotype=BIOTYPE, gene_id=Gene, transcript_id=Feature, protein_id=ENSP, HGVSp, t_ref_count, t_alt_count) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(sample = str_sub(sample, 1, 12)) %>%
  mutate(VAF = t_alt_count/(t_alt_count+t_ref_count)) %>%
  select(-t_alt_count, -t_ref_count) %>%
  filter(biotype == "protein_coding") %>%
  filter(variant_class %in% mutations_accepted) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "enst_id"))

cbio_ov2 <- cbio_ov %>%
  filter(sample %in% protein_samples$sample) %>%
  inner_join(prot_seqs[, c("protein_id", "seq")], by = "protein_id") %>%
  group_by(sample) %>%
  mutate(id_var = 1:n()) %>%
  ungroup() %>%
  select(sample, id_var, everything()) %>%
  group_by(sample, id_var) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = cbio_datasets, code = aa_code)) %>%
  unnest(cols = data) %>%
  select(-id_var)

exclude <- cbio_ov2 %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

cbio_ov2 <- cbio_ov2 %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)


# ov samples
ov_samples <- tibble(sample = unique(cbio_ov$sample), batch = "tcga-hgsc", cancer = "hgsc")



# join all samples with mutations together
tcga_samples <- bind_rows(brca_samples, coread_samples) %>%
  bind_rows(ov_samples)

write.table(tcga_samples, "./output/files/mutations_samples_tcga.txt", sep="\t", quote=F, row.names=F)



# join all mutations together
# keep only the samples with protein abundance
tcga_mutations <- bind_rows(cbio_brca2, cbio_coread2) %>%
  bind_rows(cbio_ov2)

gz1 <- gzfile("./output/files/mutations_tcga.txt.gz", "w")
write.table(tcga_mutations, gz1, sep="\t", quote=F, row.names=F)
close(gz1)


# write all variant classification - consequence pairs
consequences <- cbio_brca2 %>%
  select(variant_type, variant_class, consequence) %>%
  distinct() %>%
  bind_rows(cbio_coread2 %>% select(variant_type, variant_class, consequence) %>% distinct()) %>%
  bind_rows(cbio_ov2 %>% select(variant_type, variant_class, consequence) %>% distinct()) %>%
  distinct() %>%
  arrange(variant_type, variant_class, consequence) %>%
  mutate(consequence = str_replace_all(consequence, ",", "&"), variant_type = str_replace_all(variant_type, "ONP|TNP", "SNP")) %>%
  distinct()

write.table(consequences, file = "./output/files/tcga_mutations_var_class_conseq.txt", quote = F, row.names = F, sep = "\t")


