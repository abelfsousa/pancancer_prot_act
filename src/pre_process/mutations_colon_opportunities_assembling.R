library(tidyverse)

source("./src/utils/mutation_prot_seq_match.R")

# assembling of mutation data from "Colon Cancer Therapeutic Opportunities"
# linkedomics dataset


# cancer samples with proteomics measurements
protein_samples <- read_tsv(file = "./output/files/proteomics_samples.txt")


# protein sequences
prot_seqs <- read_tsv(file = "./output/files/protein_sequences_GRCh37.txt") %>%
  filter(state == "known" & gene_biotype == "protein_coding" & trpt_biotype == "protein_coding")


# amino acids code
aa_code <- read_tsv(file = "./data/protein/protein_sequences/aminoacids_code.txt", skip = c(1))


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
  "Splice_Site",
  "De_novo_Start_OutOfFrame",
  "De_novo_Start_InFrame",
  "Start_Codon_SNP",
  "Start_Codon_Del",
  "Start_Codon_Ins",
  "Stop_Codon_Del",
  "Stop_Codon_Ins")

mutations_accepted <- c(
  "Missense_Mutation",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "Nonsense_Mutation",
  "Nonstop_Mutation",
  "Splice_Site")


mutations_colon <- data.table::fread(file = "./data/dna/mutation/linkedomics/colon_cancer/final_oncotator_all.maf.gz") %>%
  as_tibble() %>%
  select(sample = Tumor_Sample_Barcode, consequence=Ensembl_so_term, variant_class=Variant_Classification, variant_type=Variant_Type, CHROM=Chromosome, POS=Start_position, REF=Reference_Allele, ALT=Tumor_Seq_Allele2, gene_symbol=Hugo_Symbol, biotype=gencode_transcript_type, transcript_id=Annotation_Transcript, HGVSp=HGVS_protein_change) %>%
  filter(biotype == "protein_coding") %>%
  filter(variant_class %in% mutations_accepted) %>%
  mutate(transcript_id = str_split_fixed(transcript_id, "\\.", 2)[,1]) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "enst_id")) %>%
  inner_join(prot_seqs[, c("trpt_id", "gene_id")], by = c("transcript_id" = "trpt_id"))

mutations_colon1 <- mutations_colon %>%
  filter(HGVSp == "") %>%
  inner_join(prot_seqs[, c("protein_id", "trpt_id")], by = c("transcript_id" = "trpt_id"))

mutations_colon2 <- mutations_colon %>%
  filter(!(HGVSp == "" | HGVSp == "Exception_encountered")) %>%
  separate(col = HGVSp, into = c("protein_id", "HGVSp"), sep = "\\:")

mutations_colon <- bind_rows(mutations_colon2, mutations_colon1)

mutations_colon2 <- mutations_colon %>%
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

exclude <- mutations_colon2 %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

mutations_colon2 <- mutations_colon2 %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)



# colon samples with mutation info
colon_samples <- tibble(sample = unique(mutations_colon$sample), batch = "colon-opportunities", cancer = "coread")
write.table(colon_samples, "./output/files/mutations_colon_oppt_samples.txt", sep="\t", quote=F, row.names=F)



gz1 <- gzfile("./output/files/mutations_colon_oppt.txt.gz", "w")
write.table(mutations_colon2, gz1, sep="\t", quote=F, row.names=F)
close(gz1)





