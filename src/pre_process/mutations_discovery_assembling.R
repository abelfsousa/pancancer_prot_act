library(tidyverse)


source("./src/utils/mutation_prot_seq_match.R")


# assembling of mutation data from discovery cptac studies
# gdc harmonized datasets


# cancer samples with proteomics measurements
protein_samples <- read_tsv(file = "./output/files/proteomics_samples.txt")


# protein sequences
prot_seqs <- read_tsv(file = "./output/files/protein_sequences_GRCh37.txt") %>%
  filter(state == "known" & gene_biotype == "protein_coding" & trpt_biotype == "protein_coding")


# amino acids code
aa_code <- read_tsv(file = "./data/protein/protein_sequences/aminoacids_code.txt", skip = c(1))


# all variant classification - consequence pairs from TCGA datasets
var_class_consq <- read_tsv(file = "./output/files/tcga_mutations_var_class_conseq.txt")


# canonical uniprot transcripts
canonical_trpt <- read_tsv(file = "./data/uniprot/canonical_trpt/isoform_overrides_uniprot.txt", skip = 1)
  #filter(!is.na(refseq_id)) %>%
  #mutate(refseq_id = str_split_fixed(refseq_id, "\\.", 2)[,1])


# ccrcc
gdc_ccrcc_meta <- read_tsv(file = "./data/dna/mutation/gdc/ccrcc/gdc_sample_sheet.2019-10-10.tsv") %>%
  select(file = `File Name`, sample = `Case ID`) %>%
  # C3L-00908 has two separate VCF calls corresponding to two diferent tumors. The protein abundance was calculated with C3L-00908-03 and not C3L-00908-01
  filter(file != "6600ee26-1bbe-4b4b-bf04-cd1e0785d97e.wxs.MuTect2.somatic_annotation.vcf.gz") %>%
  separate(col = sample, into = c("sample", "sample2"), sep = ", ", extra = "drop") %>%
  select(-sample2) %>%
  mutate(sample = str_replace(sample, "-", "."))
 
gdc_ccrcc <- data.table::fread(file = "./data/dna/mutation/gdc/ccrcc/all_samples_pass_variants.txt.gz") %>%
  as_tibble() %>%
  rename(file = sample) %>%
  inner_join(gdc_ccrcc_meta, by = "file") %>%
  select(-file) %>%
  select(sample, consequence=Consequence, variant_type=VARIANT_CLASS, CHROM, POS, REF, ALT, gene_symbol=SYMBOL, biotype=BIOTYPE, gene_id=Gene, transcript_id=Feature, protein_id=ENSP, HGVSp) %>%
  mutate(variant_type = if_else(variant_type == "SNV", "SNP", if_else(variant_type == "deletion", "DEL", "INS"))) %>%
  filter(biotype == "protein_coding") %>%
  inner_join(var_class_consq, by = c("consequence", "variant_type")) %>%
  select(sample, consequence, variant_type, variant_class, everything()) %>%
  separate(col = "HGVSp", into = c("protein_id2", "HGVSp"), sep = ":", fill = "right") %>%
  select(-protein_id2) %>%
  mutate(HGVSp = str_replace_na(HGVSp, "")) %>%
  filter(str_detect(HGVSp, "^c.", negate = T)) %>%
  group_by(sample, variant_type, variant_class, consequence, CHROM, POS, REF, ALT, gene_symbol, biotype) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "enst_id")) %>%
  ungroup() %>%
  distinct()

gdc_ccrcc2 <- gdc_ccrcc %>%
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

exclude <- gdc_ccrcc2 %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

gdc_ccrcc2 <- gdc_ccrcc2 %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)

# ccrcc samples
ccrcc_samples <- tibble(sample = unique(gdc_ccrcc$sample), batch = "discovery-ccrcc", cancer = "ccrcc")



# luad
gdc_luad_meta <- read_tsv(file = "./data/dna/mutation/gdc/luad/gdc_sample_sheet.2019-10-10.tsv") %>%
  select(file = `File Name`, sample = `Case ID`) %>%
  separate(col = sample, into = c("sample", "sample2"), sep = ", ", extra = "drop") %>%
  select(-sample2) %>%
  mutate(sample = str_replace(sample, "-", "."))

gdc_luad <- data.table::fread(file = "./data/dna/mutation/gdc/luad/all_samples_pass_variants.txt.gz") %>%
  as_tibble() %>%
  rename(file = sample) %>%
  inner_join(gdc_luad_meta, by = "file") %>%
  select(-file) %>%
  select(sample, consequence=Consequence, variant_type=VARIANT_CLASS, CHROM, POS, REF, ALT, gene_symbol=SYMBOL, biotype=BIOTYPE, gene_id=Gene, transcript_id=Feature, protein_id=ENSP, HGVSp) %>%
  mutate(variant_type = if_else(variant_type == "SNV", "SNP", if_else(variant_type == "deletion", "DEL", "INS"))) %>%
  filter(biotype == "protein_coding") %>%
  inner_join(var_class_consq, by = c("consequence", "variant_type")) %>%
  select(sample, consequence, variant_type, variant_class, everything()) %>%
  separate(col = "HGVSp", into = c("protein_id2", "HGVSp"), sep = ":", fill = "right") %>%
  select(-protein_id2) %>%
  mutate(HGVSp = str_replace_na(HGVSp, "")) %>%
  filter(str_detect(HGVSp, "^c.", negate = T)) %>%
  group_by(sample, variant_type, variant_class, consequence, CHROM, POS, REF, ALT, gene_symbol, biotype) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "enst_id")) %>%
  ungroup() %>%
  distinct()

gdc_luad2 <- gdc_luad %>%
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

exclude <- gdc_luad2 %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

gdc_luad2 <- gdc_luad2 %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)

# luad samples
luad_samples <- tibble(sample = unique(gdc_luad$sample), batch = "discovery-luad", cancer = "luad")



# ucec
gdc_ucec_meta <- read_tsv(file = "./data/dna/mutation/gdc/ucec/gdc_sample_sheet.2019-10-10.tsv") %>%
  select(file = `File Name`, sample = `Case ID`) %>%
  # C3N-01825 has two separate VCF calls corresponding to two diferent tumors. The protein abundance was calculated with C3N-01825-01 and not C3N-01825-03
  filter(file != "7e6c5e4b-861a-4268-9f51-019acaa36d87.wxs.MuTect2.somatic_annotation.vcf.gz") %>%
  separate(col = sample, into = c("sample", "sample2"), sep = ", ", extra = "drop") %>%
  select(-sample2) %>%
  mutate(sample = str_replace(sample, "-", "."))

gdc_ucec <- data.table::fread(file = "./data/dna/mutation/gdc/ucec/all_samples_pass_variants.txt.gz") %>%
  as_tibble() %>%
  rename(file = sample) %>%
  inner_join(gdc_ucec_meta, by = "file") %>%
  select(-file) %>%
  select(sample, consequence=Consequence, variant_type=VARIANT_CLASS, CHROM, POS, REF, ALT, gene_symbol=SYMBOL, biotype=BIOTYPE, gene_id=Gene, transcript_id=Feature, protein_id=ENSP, HGVSp) %>%
  mutate(variant_type = if_else(variant_type == "SNV", "SNP", if_else(variant_type == "deletion", "DEL", "INS"))) %>%
  filter(biotype == "protein_coding") %>%
  inner_join(var_class_consq, by = c("consequence", "variant_type")) %>%
  select(sample, consequence, variant_type, variant_class, everything()) %>%
  separate(col = "HGVSp", into = c("protein_id2", "HGVSp"), sep = ":", fill = "right") %>%
  select(-protein_id2) %>%
  mutate(HGVSp = str_replace_na(HGVSp, "")) %>%
  filter(str_detect(HGVSp, "^c.", negate = T)) %>%
  group_by(sample, variant_type, variant_class, consequence, CHROM, POS, REF, ALT, gene_symbol, biotype) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "enst_id")) %>%
  ungroup() %>%
  distinct()

gdc_ucec2 <- gdc_ucec %>%
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

exclude <- gdc_ucec2 %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

gdc_ucec2 <- gdc_ucec2 %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)

# ucec samples
ucec_samples <- tibble(sample = unique(gdc_ucec$sample), batch = "discovery-ucec", cancer = "ucec")



# join all samples with mutations together
discovery_samples <- bind_rows(ccrcc_samples, luad_samples) %>%
  bind_rows(ucec_samples)

write.table(discovery_samples, "./output/files/mutations_samples_discovery.txt", sep="\t", quote=F, row.names=F)



# join all mutations together (only the samples with protein abundance)

discovery_mutations <- bind_rows(
  gdc_ccrcc2, gdc_luad2) %>%
  bind_rows(gdc_ucec2)

gz1 <- gzfile("./output/files/mutations_discovery.txt.gz", "w")
write.table(discovery_mutations, gz1, sep="\t", quote=F, row.names=F)
close(gz1)



