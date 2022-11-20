library(tidyverse)


source("./src/utils/mutation_prot_seq_match.R")


# assembling of mutation data from cbttc study
# pedcbioportal datasets


# cancer samples with proteomics measurements
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt")

proteomics_cbttc <- proteomics %>%
  filter(batch == "cptac-cbttc")


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


# clinical data from CBTTC proteomics dataset
cbttc_clinical_prot <- data.table::fread("./data/protein/cptac_new/cbttc/S047_Pediatric_Brain_Cancer_Clinical_Data_r1.txt", check.names = T) %>%
  as_tibble() %>%
  filter(Kids.First.ID %in% proteomics_cbttc$sample) %>%
  select(Kids.First.ID, Clinical.Event.Id, diagnosis_type, diagnosis)


# patients with two different tumor samples from CBTTC proteomics dataset
cptac_cbttc2tumours <- read_tsv(file = "./output/files/cptac_cbttc2tumours.txt")




# Atypical Teratoid Rhabdoid Tumor (ATRT)
atrt_meta <- read_tsv("./data/dna/mutation/pedcbioportal/atrt_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, diagnosis=CANCER_TYPE)

atrt_meta_samples <- data.table::fread("./data/dna/mutation/pedcbioportal/atrt_cbttc/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  filter(Variant_Classification %in% mutations_accepted) %>%
  select(Tumor_Sample_Barcode) %>%
  distinct() %>%
  inner_join(atrt_meta[, c("Kids.First.ID", "Clinical.Event.Id", "diagnosis")], by = c("Tumor_Sample_Barcode" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, diagnosis) %>%
  tally() %>%
  ungroup()

atrt_meta_samples2 <- atrt_meta %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "diagnosis"))

atrt_mutations <- data.table::fread("./data/dna/mutation/pedcbioportal/atrt_cbttc/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  select(sample = Tumor_Sample_Barcode, consequence=Consequence, variant_class=Variant_Classification, variant_type=Variant_Type, CHROM=Chromosome, POS=Start_Position, REF=Reference_Allele, ALT=Tumor_Seq_Allele2, gene_symbol=Hugo_Symbol, biotype=BIOTYPE, gene_id=Gene, transcript_id=Feature, protein_id=ENSP, HGVSp, t_ref_count, t_alt_count) %>%
  mutate(t_alt_count = replace(t_alt_count, t_alt_count==".", NA), t_ref_count = replace(t_ref_count, t_ref_count==".", NA)) %>%
  mutate(t_alt_count = as.numeric(t_alt_count), t_ref_count = as.numeric(t_ref_count)) %>%
  mutate(VAF = t_alt_count/(t_alt_count+t_ref_count)) %>%
  select(-t_alt_count, -t_ref_count) %>%
  filter(biotype == "protein_coding") %>%
  filter(variant_class %in% mutations_accepted) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "enst_id")) %>%
  inner_join(atrt_meta_samples2[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  select(sample = Kids.First.ID, everything()) %>%
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

exclude <- atrt_mutations %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

atrt_mutations <- atrt_mutations %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)



# Craniopharyngioma
cranio_meta <- read_tsv("./data/dna/mutation/pedcbioportal/craniopharyngioma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, diagnosis=CANCER_TYPE)

cranio_meta_samples <- data.table::fread("./data/dna/mutation/pedcbioportal/craniopharyngioma_cbttc/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  filter(Variant_Classification %in% mutations_accepted) %>%
  select(Tumor_Sample_Barcode) %>%
  distinct() %>%
  inner_join(cranio_meta[, c("Kids.First.ID", "Clinical.Event.Id", "diagnosis")], by = c("Tumor_Sample_Barcode" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, diagnosis) %>%
  tally() %>%
  ungroup()

cranio_meta_samples2 <- cranio_meta %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "diagnosis"))

cranio_mutations <- data.table::fread("./data/dna/mutation/pedcbioportal/craniopharyngioma_cbttc/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  select(sample = Tumor_Sample_Barcode, consequence=Consequence, variant_class=Variant_Classification, variant_type=Variant_Type, CHROM=Chromosome, POS=Start_Position, REF=Reference_Allele, ALT=Tumor_Seq_Allele2, gene_symbol=Hugo_Symbol, biotype=BIOTYPE, gene_id=Gene, transcript_id=Feature, protein_id=ENSP, HGVSp, t_ref_count, t_alt_count) %>%
  mutate(t_alt_count = replace(t_alt_count, t_alt_count==".", NA), t_ref_count = replace(t_ref_count, t_ref_count==".", NA)) %>%
  mutate(t_alt_count = as.numeric(t_alt_count), t_ref_count = as.numeric(t_ref_count)) %>%
  mutate(VAF = t_alt_count/(t_alt_count+t_ref_count)) %>%
  select(-t_alt_count, -t_ref_count) %>%
  filter(biotype == "protein_coding") %>%
  filter(variant_class %in% mutations_accepted) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "enst_id")) %>%
  inner_join(cranio_meta_samples2[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  select(sample = Kids.First.ID, everything()) %>%
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

exclude <- cranio_mutations %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

cranio_mutations <- cranio_mutations %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)



# Ependymoma
epend_meta <- read_tsv("./data/dna/mutation/pedcbioportal/ependymoma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, diagnosis=CANCER_TYPE)

epend_meta_samples <- data.table::fread("./data/dna/mutation/pedcbioportal/ependymoma_cbttc/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  filter(Variant_Classification %in% mutations_accepted) %>%
  select(Tumor_Sample_Barcode) %>%
  distinct() %>%
  inner_join(epend_meta[, c("Kids.First.ID", "Clinical.Event.Id", "diagnosis")], by = c("Tumor_Sample_Barcode" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, diagnosis) %>%
  tally() %>%
  ungroup()

epend_meta_samples2 <- epend_meta %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "diagnosis"))

epend_mutations <- data.table::fread("./data/dna/mutation/pedcbioportal/ependymoma_cbttc/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  select(sample = Tumor_Sample_Barcode, consequence=Consequence, variant_class=Variant_Classification, variant_type=Variant_Type, CHROM=Chromosome, POS=Start_Position, REF=Reference_Allele, ALT=Tumor_Seq_Allele2, gene_symbol=Hugo_Symbol, biotype=BIOTYPE, gene_id=Gene, transcript_id=Feature, protein_id=ENSP, HGVSp, t_ref_count, t_alt_count) %>%
  mutate(t_alt_count = replace(t_alt_count, t_alt_count==".", NA), t_ref_count = replace(t_ref_count, t_ref_count==".", NA)) %>%
  mutate(t_alt_count = as.numeric(t_alt_count), t_ref_count = as.numeric(t_ref_count)) %>%
  mutate(VAF = t_alt_count/(t_alt_count+t_ref_count)) %>%
  select(-t_alt_count, -t_ref_count) %>%
  filter(biotype == "protein_coding") %>%
  filter(variant_class %in% mutations_accepted) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "enst_id")) %>%
  inner_join(epend_meta_samples2[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  select(sample = Kids.First.ID, everything()) %>%
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

exclude <- epend_mutations %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

epend_mutations <- epend_mutations %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)



# Ganglioglioma
gangli_meta <- read_tsv("./data/dna/mutation/pedcbioportal/ganglioglioma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, diagnosis=CANCER_TYPE)

gangli_meta_samples <- data.table::fread("./data/dna/mutation/pedcbioportal/ganglioglioma_cbttc/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  filter(Variant_Classification %in% mutations_accepted) %>%
  select(Tumor_Sample_Barcode) %>%
  distinct() %>%
  inner_join(gangli_meta[, c("Kids.First.ID", "Clinical.Event.Id", "diagnosis")], by = c("Tumor_Sample_Barcode" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, diagnosis) %>%
  tally() %>%
  ungroup()

gangli_meta_samples2 <- gangli_meta %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "diagnosis"))

gangli_mutations <- data.table::fread("./data/dna/mutation/pedcbioportal/ganglioglioma_cbttc/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  select(sample = Tumor_Sample_Barcode, consequence=Consequence, variant_class=Variant_Classification, variant_type=Variant_Type, CHROM=Chromosome, POS=Start_Position, REF=Reference_Allele, ALT=Tumor_Seq_Allele2, gene_symbol=Hugo_Symbol, biotype=BIOTYPE, gene_id=Gene, transcript_id=Feature, protein_id=ENSP, HGVSp, t_ref_count, t_alt_count) %>%
  mutate(t_alt_count = replace(t_alt_count, t_alt_count==".", NA), t_ref_count = replace(t_ref_count, t_ref_count==".", NA)) %>%
  mutate(t_alt_count = as.numeric(t_alt_count), t_ref_count = as.numeric(t_ref_count)) %>%
  mutate(VAF = t_alt_count/(t_alt_count+t_ref_count)) %>%
  select(-t_alt_count, -t_ref_count) %>%
  filter(biotype == "protein_coding") %>%
  filter(variant_class %in% mutations_accepted) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "enst_id")) %>%
  inner_join(gangli_meta_samples2[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  select(sample = Kids.First.ID, everything()) %>%
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

exclude <- gangli_mutations %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

gangli_mutations <- gangli_mutations %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)




# Medulloblastoma
medullo_meta <- read_tsv("./data/dna/mutation/pedcbioportal/medulloblastoma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, diagnosis=CANCER_TYPE)

medullo_meta_samples <- data.table::fread("./data/dna/mutation/pedcbioportal/medulloblastoma_cbttc/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  filter(Variant_Classification %in% mutations_accepted) %>%
  select(Tumor_Sample_Barcode) %>%
  distinct() %>%
  inner_join(medullo_meta[, c("Kids.First.ID", "Clinical.Event.Id", "diagnosis")], by = c("Tumor_Sample_Barcode" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, diagnosis) %>%
  tally() %>%
  ungroup()

medullo_meta_samples2 <- medullo_meta %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "diagnosis"))

medullo_mutations <- data.table::fread("./data/dna/mutation/pedcbioportal/medulloblastoma_cbttc/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  select(sample = Tumor_Sample_Barcode, consequence=Consequence, variant_class=Variant_Classification, variant_type=Variant_Type, CHROM=Chromosome, POS=Start_Position, REF=Reference_Allele, ALT=Tumor_Seq_Allele2, gene_symbol=Hugo_Symbol, biotype=BIOTYPE, gene_id=Gene, transcript_id=Feature, protein_id=ENSP, HGVSp, t_ref_count, t_alt_count) %>%
  mutate(t_alt_count = replace(t_alt_count, t_alt_count==".", NA), t_ref_count = replace(t_ref_count, t_ref_count==".", NA)) %>%
  mutate(t_alt_count = as.numeric(t_alt_count), t_ref_count = as.numeric(t_ref_count)) %>%
  mutate(VAF = t_alt_count/(t_alt_count+t_ref_count)) %>%
  select(-t_alt_count, -t_ref_count) %>%
  filter(biotype == "protein_coding") %>%
  filter(variant_class %in% mutations_accepted) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "enst_id")) %>%
  inner_join(medullo_meta_samples2[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  select(sample = Kids.First.ID, everything()) %>%
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

exclude <- medullo_mutations %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

medullo_mutations <- medullo_mutations %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)




# High-grade glioma
hgg_meta <- read_tsv("./data/dna/mutation/pedcbioportal/ped_high_grade_glioma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, diagnosis=CANCER_TYPE) %>%
  mutate(diagnosis = str_replace(diagnosis, " \\(WHO grade III/IV\\)", ""))

hgg_meta_samples <- data.table::fread("./data/dna/mutation/pedcbioportal/ped_high_grade_glioma_cbttc/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  filter(Variant_Classification %in% mutations_accepted) %>%
  select(Tumor_Sample_Barcode) %>%
  distinct() %>%
  inner_join(hgg_meta[, c("Kids.First.ID", "Clinical.Event.Id", "diagnosis")], by = c("Tumor_Sample_Barcode" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, diagnosis) %>%
  tally() %>%
  ungroup()

hgg_meta_samples2 <- hgg_meta %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "diagnosis"))

hgg_mutations <- data.table::fread("./data/dna/mutation/pedcbioportal/ped_high_grade_glioma_cbttc/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  select(sample = Tumor_Sample_Barcode, consequence=Consequence, variant_class=Variant_Classification, variant_type=Variant_Type, CHROM=Chromosome, POS=Start_Position, REF=Reference_Allele, ALT=Tumor_Seq_Allele2, gene_symbol=Hugo_Symbol, biotype=BIOTYPE, gene_id=Gene, transcript_id=Feature, protein_id=ENSP, HGVSp, t_ref_count, t_alt_count) %>%
  mutate(t_alt_count = replace(t_alt_count, t_alt_count==".", NA), t_ref_count = replace(t_ref_count, t_ref_count==".", NA)) %>%
  mutate(t_alt_count = as.numeric(t_alt_count), t_ref_count = as.numeric(t_ref_count)) %>%
  mutate(VAF = t_alt_count/(t_alt_count+t_ref_count)) %>%
  select(-t_alt_count, -t_ref_count) %>%
  filter(biotype == "protein_coding") %>%
  filter(variant_class %in% mutations_accepted) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "enst_id")) %>%
  inner_join(hgg_meta_samples2[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  select(sample = Kids.First.ID, everything()) %>%
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

exclude <- hgg_mutations %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

hgg_mutations <- hgg_mutations %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)




# Low-grade glioma
lgg_meta <- read_tsv("./data/dna/mutation/pedcbioportal/ped_low_grade_glioma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, diagnosis=CANCER_TYPE) %>%
  mutate(diagnosis = str_replace(diagnosis, " \\(WHO grade I/II\\)", ""))

lgg_meta_samples <- data.table::fread("./data/dna/mutation/pedcbioportal/ped_low_grade_glioma_cbttc/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  filter(Variant_Classification %in% mutations_accepted) %>%
  select(Tumor_Sample_Barcode) %>%
  distinct() %>%
  inner_join(lgg_meta[, c("Kids.First.ID", "Clinical.Event.Id", "diagnosis")], by = c("Tumor_Sample_Barcode" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, diagnosis) %>%
  tally() %>%
  ungroup()

lgg_meta_samples2 <- lgg_meta %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "diagnosis"))

lgg_mutations <- data.table::fread("./data/dna/mutation/pedcbioportal/ped_low_grade_glioma_cbttc/data_mutations_extended.txt.gz") %>%
  as_tibble() %>%
  select(sample = Tumor_Sample_Barcode, consequence=Consequence, variant_class=Variant_Classification, variant_type=Variant_Type, CHROM=Chromosome, POS=Start_Position, REF=Reference_Allele, ALT=Tumor_Seq_Allele2, gene_symbol=Hugo_Symbol, biotype=BIOTYPE, gene_id=Gene, transcript_id=Feature, protein_id=ENSP, HGVSp, t_ref_count, t_alt_count) %>%
  mutate(t_alt_count = replace(t_alt_count, t_alt_count==".", NA), t_ref_count = replace(t_ref_count, t_ref_count==".", NA)) %>%
  mutate(t_alt_count = as.numeric(t_alt_count), t_ref_count = as.numeric(t_ref_count)) %>%
  mutate(VAF = t_alt_count/(t_alt_count+t_ref_count)) %>%
  select(-t_alt_count, -t_ref_count) %>%
  filter(biotype == "protein_coding") %>%
  filter(variant_class %in% mutations_accepted) %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene_symbol" = "gene_name", "transcript_id" = "enst_id")) %>%
  inner_join(lgg_meta_samples2[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  select(sample = Kids.First.ID, everything()) %>%
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

exclude <- lgg_mutations %>%
  filter(same == 0) %>%
  group_by(gene_symbol, protein_id) %>%
  tally() %>%
  ungroup() %>%
  arrange(desc(n)) %>%
  pull(protein_id)

lgg_mutations <- lgg_mutations %>%
  filter(same == 1) %>%
  select(-same)
  #filter(!protein_id %in% exclude)




# cbttc samples
cbttc_samples <- atrt_meta_samples %>%
  select(-n) %>%
  bind_rows(cranio_meta_samples[, -c(3)]) %>%
  bind_rows(epend_meta_samples[, -c(3)]) %>%
  bind_rows(gangli_meta_samples[, -c(3)]) %>%
  bind_rows(medullo_meta_samples[, -c(3)]) %>%
  bind_rows(hgg_meta_samples[, -c(3)]) %>%
  bind_rows(lgg_meta_samples[, -c(3)]) %>%
  mutate(batch = "cbttc") %>%
  select(sample = Kids.First.ID, batch, cancer = diagnosis)

write.table(cbttc_samples, "./output/files/mutations_samples_cbttc.txt", sep="\t", quote=F, row.names=F)



# join all mutations together
cbttc_mutations <- bind_rows(
  atrt_mutations,
  cranio_mutations,
  epend_mutations,
  gangli_mutations,
  medullo_mutations,
  hgg_mutations,
  lgg_mutations)

gz1 <- gzfile("./output/files/mutations_cbttc.txt.gz", "w")
write.table(cbttc_mutations, gz1, sep="\t", quote=F, row.names=F)
close(gz1)
