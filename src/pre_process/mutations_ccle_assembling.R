library(tidyverse)


source("./src/utils/mutation_prot_seq_match.R")


# assembling of mutation data from ccle cell lines
# depmap dataset



# cancer samples with proteomics measurements
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt")

cell_lines <- proteomics %>%
  filter(str_detect(batch, "cell-lines"))


# protein sequences
prot_seqs <- read_tsv(file = "./output/files/protein_sequences_GRCh37.txt") %>%
  filter(state == "known" & gene_biotype == "protein_coding" & trpt_biotype == "protein_coding")


# canonical uniprot transcripts
canonical_trpt <- read_tsv(file = "./data/uniprot/canonical_trpt/isoform_overrides_uniprot.txt", skip = 1)
  #filter(!is.na(refseq_id)) %>%
  #mutate(refseq_id = str_split_fixed(refseq_id, "\\.", 2)[,1])


# mutations to select
mutations_accepted <- c(
  "Missense_Mutation",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Nonsense_Mutation",
  "Nonstop_Mutation",
  "Splice_Site")



# depmap metadata
depmap_meta <- data.table::fread(file = "./data/dna/mutation/depmap/sample_info.csv") %>%
  as_tibble() %>%
  select(DepMap_ID, stripped_cell_line_name, `CCLE Name`, disease) %>%
  filter(disease %in% c("Breast Cancer", "Colon/Colorectal Cancer")) %>%
  mutate(disease = if_else(disease == "Breast Cancer", "brca", "coread")) %>%
  mutate(stripped_cell_line_name = str_replace(stripped_cell_line_name, "MDAMB175VII", "MDAMB175")) %>%
  mutate(stripped_cell_line_name = str_replace(stripped_cell_line_name, "MDAMB134VI", "MDAMB134"))

#cell_lines %>% filter(!(sample %in% depmap_meta$stripped_cell_line_name))
#depmap_meta %>% filter(str_detect(stripped_cell_line_name, "MDAMB1"))


depmap_mutations <- data.table::fread(file = "./data/dna/mutation/depmap/CCLE_mutations.csv.gz") %>%
  as_tibble() %>%
  select(DepMap_ID, variant_class=Variant_Classification, variant_type=Variant_Type, CHROM=Chromosome, POS=Start_position, REF=Reference_Allele, ALT=Tumor_Seq_Allele1, gene_symbol=Hugo_Symbol, transcript_id=Annotation_Transcript, HGVSp=Protein_Change) %>%
  filter(variant_class %in% mutations_accepted) %>%
  mutate(transcript_id = str_split_fixed(transcript_id, "\\.", 2)[,1]) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "enst_id")) %>%
  inner_join(depmap_meta[, c(1:2)], by = "DepMap_ID") %>%
  select(DepMap_ID, stripped_cell_line_name, everything()) %>%
  mutate(HGVSp = str_replace_na(HGVSp, ""))


# cell lines
ccle_samples <- tibble(sample = unique(depmap_mutations$stripped_cell_line_name), batch = "ccle") %>%
  inner_join(depmap_meta[, c("stripped_cell_line_name", "disease")], by = c("sample" = "stripped_cell_line_name")) %>%
  rename(cancer = disease)

write.table(ccle_samples, "./output/files/mutations_samples_ccle.txt", sep="\t", quote=F, row.names=F)


depmap_mutations2 <- depmap_mutations %>%
  filter(stripped_cell_line_name %in% cell_lines$sample) %>%
  select(-DepMap_ID) %>%
  rename(sample = stripped_cell_line_name) %>%
  inner_join(prot_seqs[, c("gene_id", "trpt_id", "protein_id", "seq")], by = c("transcript_id" = "trpt_id")) %>%
  group_by(sample) %>%
  mutate(id_var = 1:n()) %>%
  ungroup() %>%
  select(sample, id_var, everything()) %>%
  group_by(sample, id_var) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = ccle_dataset)) %>%
  unnest(cols = data) %>%
  filter(!is.na(same)) %>%
  select(-id_var)

exclude <- depmap_mutations2 %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

depmap_mutations2 <- depmap_mutations2 %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)



gz1 <- gzfile("./output/files/mutations_ccle.txt.gz", "w")
write.table(depmap_mutations2, gz1, sep="\t", quote=F, row.names=F)
close(gz1)



