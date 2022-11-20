library(tidyverse)

source("./src/utils/mutation_prot_seq_match.R")


# assembling of mutation data from liver cancer cptac study: "Proteogenomic Characterization of HBV-Related Hepatocellular Carcinoma"


# samples with protein measurements
protein_samples <- read_tsv(file = "./output/files/proteomics_samples.txt")


# code to replace original mutation classification
var_class_type = tibble(
  variant_class1 = c("nonsynonymous SNV", "stopgain", "stoploss", "frameshift deletion", "frameshift insertion", "nonframeshift deletion", "nonframeshift insertion"),
  variant_class = c("Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins"),
  variant_type = c("SNP", "SNP", "SNP", "DEL", "INS", "DEL", "INS"))


# protein sequences
prot_seqs <- read_tsv(file = "./output/files/protein_sequences_GRCh37.txt") %>%
  filter(state == "known" & gene_biotype == "protein_coding" & trpt_biotype == "protein_coding")


# canonical uniprot transcripts
canonical_trpt <- read_tsv(file = "./data/uniprot/canonical_trpt/isoform_overrides_uniprot.txt", skip = 1) %>%
  filter(!is.na(refseq_id)) %>%
  mutate(refseq_id = str_split_fixed(refseq_id, "\\.", 2)[,1])


# mutation data
liver_mutations <- read_tsv(file = "./data/dna/mutation/liver_cancer_paper/wes_somatic_variants.txt") %>%
  select(-c(End, VariantType, VariantLocal, NormalDepth, TumorDepth)) %>%
  rename(variant_class1 = VariantClass) %>%
  mutate(GenomicAnnotation = map(.x = GenomicAnnotation, ~ as.character(str_split(.x, ",", simplify = T)) )) %>%
  unnest(cols = GenomicAnnotation) %>%
  filter(GenomicAnnotation != "") %>%
  separate(col = GenomicAnnotation, into = c("Gene", "transcript_id", "exon_n", "HGVSc", "HGVSp"), sep = ":", fill = "right") %>%
  select(-Gene, -exon_n, -HGVSc) %>%
  filter(!(is.na(transcript_id) | is.na(HGVSp))) %>%
  select(sample=Sample, variant_class1, CHROM = Chr, POS=Start, REF=Ref, ALT=Alt, gene_symbol=gene, transcript_id, HGVSp, VAF) %>%
  inner_join(var_class_type, by="variant_class1") %>%
  select(-variant_class1) %>%
  select(sample, variant_class, variant_type, everything()) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id", "refseq_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "refseq_id")) %>%
  select(-transcript_id) %>%
  rename(transcript_id = enst_id) %>%
  inner_join(prot_seqs[, c("gene_id", "trpt_id", "protein_id", "seq")], by = c("transcript_id" = "trpt_id")) %>%
  filter(sample %in% protein_samples$sample) %>%
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

exclude <- liver_mutations %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

liver_mutations <- liver_mutations %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)

liver_mutations <- liver_mutations %>%
  mutate(VAF = as.numeric(str_replace(VAF, "%", ""))) %>%
  mutate(VAF = VAF/100)


gz1 <- gzfile("./output/files/mutations_hcc.txt.gz", "w")
write.table(liver_mutations, gz1, sep="\t", quote=F, row.names=F)
close(gz1)


# samples
samples <- tibble(sample = unique(liver_mutations$sample), batch = "HBV-HCC", cancer = "HBV-HCC")
write.table(samples, "./output/files/mutations_hcc_samples.txt", sep="\t", quote=F, row.names=F)
