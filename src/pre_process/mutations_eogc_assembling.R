library(tidyverse)

source("./src/utils/mutation_prot_seq_match.R")


# assembling of mutation data from eogc cptac study: "Proteogenomic Characterization of Human Early-Onset Gastric Cancer"


# eogc samples with protein measurements
eogc_proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt") %>%
  filter(batch == "cptac-GC") %>%
  mutate(sample2 = str_extract(sample, "^T[0-9]{1,4}")) %>%
  mutate(sample2 = str_c(str_extract(sample2, "[0-9]+$"), "T", sep = "")) %>%
  select(sample1=sample, sample2, batch, cancer)


# code to replace original mutation classification
var_class_type = tibble(
  variant_class1 = c("synonymous SNV", "nonsynonymous SNV", "stopgain SNV", "stoploss SNV", "frameshift deletion", "frameshift insertion", "nonframeshift deletion", "nonframeshift insertion", "splicing"),
  variant_class = c("Silent", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site"),
  variant_type = c("SNP", "SNP", "SNP", "SNP", "DEL", "INS", "DEL", "INS", "SNP"))


# protein sequences
prot_seqs <- read_tsv(file = "./output/files/protein_sequences_GRCh37.txt") %>%
  filter(state == "known" & gene_biotype == "protein_coding" & trpt_biotype == "protein_coding")


# canonical uniprot transcripts
canonical_trpt <- read_tsv(file = "./data/uniprot/canonical_trpt/isoform_overrides_uniprot.txt", skip = 1) %>%
  filter(!is.na(refseq_id)) %>%
  mutate(refseq_id = str_split_fixed(refseq_id, "\\.", 2)[,1])


# mutation data
eogc_mutations1 <- read_tsv(file = "./data/dna/mutation/geo/eogc/WES_somaticVariants.txt") %>%
  filter(str_detect(MutationType, "splicing", negate = T)) %>%
  #mutate(GenomicAnnotation = str_split(GenomicAnnotation, ",")) %>%
  #mutate(GenomicAnnotation = map(.x = GenomicAnnotation, .f = function(x) as.character(str_split(x, ",", simplify = T)) )) %>%
  mutate(GenomicAnnotation = map(.x = GenomicAnnotation, ~ as.character(str_split(.x, ",", simplify = T)) )) %>%
  unnest(cols = GenomicAnnotation) %>%
  filter(GenomicAnnotation != "") %>%
  separate(col = GenomicAnnotation, into = c("gene", "transcript_id", "exon_n", "HGVSc", "HGVSp"), sep = ":") %>%
  select(-gene, -exon_n, -HGVSc) %>%
  select(sample = Patient, variant_class1 = MutationType, CHROM = CHR, POS, REF, ALT, gene_symbol = Gene, transcript_id, HGVSp) %>%
  inner_join(var_class_type, by="variant_class1") %>%
  select(-variant_class1) %>%
  select(sample, variant_class, variant_type, everything()) %>%
  filter(variant_class != "Silent") %>%
  #group_by(sample, variant_class, variant_type, CHROM, POS, REF, ALT, gene_symbol) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id", "refseq_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "refseq_id")) %>%
  #ungroup() %>%
  select(-transcript_id) %>%
  rename(transcript_id = enst_id) %>%
  inner_join(prot_seqs[, c("gene_id", "trpt_id", "protein_id", "seq")], by = c("transcript_id" = "trpt_id"))
  
  

eogc_mutations2 <- read_tsv(file = "./data/dna/mutation/geo/eogc/WES_somaticVariants.txt") %>%
  filter(str_detect(MutationType, "splicing", negate = F)) %>%
  mutate(GenomicAnnotation = str_extract(GenomicAnnotation, "\\(.+\\)")) %>%
  mutate(GenomicAnnotation = str_replace_all(GenomicAnnotation, "\\(|\\)", "")) %>%
  mutate(GenomicAnnotation = map(.x = GenomicAnnotation, ~ as.character(str_split(.x, ",", simplify = T)) )) %>%
  unnest(cols = GenomicAnnotation) %>%
  filter(GenomicAnnotation != "") %>%
  separate(col = GenomicAnnotation, into = c("transcript_id", "exon_n", "HGVSc"), sep = ":") %>%
  select(-exon_n, -HGVSc) %>%
  mutate(HGVSp = "") %>%
  select(sample = Patient, variant_class1 = MutationType, CHROM = CHR, POS, REF, ALT, gene_symbol = Gene, transcript_id, HGVSp) %>%
  inner_join(var_class_type, by="variant_class1") %>%
  select(-variant_class1) %>%
  select(sample, variant_class, variant_type, everything()) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id", "refseq_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "refseq_id")) %>%
  select(-transcript_id) %>%
  rename(transcript_id = enst_id) %>%
  inner_join(prot_seqs[, c("gene_id", "trpt_id", "protein_id", "seq")], by = c("transcript_id" = "trpt_id"))



eogc_mutations <- bind_rows(eogc_mutations1, eogc_mutations2) %>%
  inner_join(eogc_proteomics[, c("sample1", "sample2")], by = c("sample" = "sample2")) %>%
  select(-sample) %>%
  select(sample = sample1, everything()) %>%
  group_by(sample) %>%
  mutate(id_var = 1:n()) %>%
  ungroup() %>%
  select(sample, id_var, everything()) %>%
  group_by(sample, id_var) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = gc_hcc_datasets)) %>%
  unnest(cols = data) %>%
  filter(!is.na(same)) %>%
  select(-id_var)

exclude <- eogc_mutations %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

eogc_mutations <- eogc_mutations %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)


gz1 <- gzfile("./output/files/mutations_eogc.txt.gz", "w")
write.table(eogc_mutations, gz1, sep="\t", quote=F, row.names=F)
close(gz1)


# samples
samples <- tibble(sample = unique(eogc_mutations$sample), batch = "GC", cancer = "eogc")
write.table(samples, "./output/files/mutations_eogc_samples.txt", sep="\t", quote=F, row.names=F)
