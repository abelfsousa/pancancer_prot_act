library(tidyverse)


# assembling of CBTTC rna-seq data
# FPKM RSEM expression values


# cancer samples with proteomics measurements
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt")

proteomics_cbttc <- proteomics %>%
  filter(batch == "cptac-cbttc")


# clinical data from CBTTC proteomics dataset
cbttc_clinical_prot <- data.table::fread("./data/protein/cptac_new/cbttc/S047_Pediatric_Brain_Cancer_Clinical_Data_r1.txt", check.names = T) %>%
  as_tibble() %>%
  filter(Kids.First.ID %in% proteomics_cbttc$sample) %>%
  select(Kids.First.ID, Clinical.Event.Id, diagnosis_type, diagnosis)


# patients with two different tumor samples from CBTTC proteomics dataset
cptac_cbttc2tumours <- read_tsv(file = "./output/files/cptac_cbttc2tumours.txt")


# protein coding genes
protein_coding <- read_tsv("./output/files/gencode.v19.annotation.txt") %>%
  select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding")



# Atypical Teratoid Rhabdoid Tumor (ATRT)
atrt_rna_clinical <- read_tsv("./data/rna/cbttc/atrt_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, cancer=CANCER_TYPE)

atrt_rna_clinical_samples <- read_tsv("./data/rna/cbttc/atrt_cbttc/data_rna_seq_v2_mrna.txt") %>%
  rename(gene = Hugo_Symbol) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  select(sample) %>%
  distinct() %>%
  inner_join(atrt_rna_clinical[, c("Kids.First.ID", "Clinical.Event.Id", "cancer")], by = c("sample" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, cancer) %>%
  tally() %>%
  ungroup()

#atrt_rna_clinical_samples <- atrt_rna_clinical %>%
#  group_by(Kids.First.ID, cancer) %>%
#  tally() %>%
#  ungroup()

atrt_rna_clinical_prot <- atrt_rna_clinical %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "cancer" = "diagnosis")) %>%
  arrange(Kids.First.ID)

atrt_rna <- read_tsv("./data/rna/cbttc/atrt_cbttc/data_rna_seq_v2_mrna.txt") %>%
  rename(gene = Hugo_Symbol) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  inner_join(atrt_rna_clinical_prot[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  spread(key = "Kids.First.ID", value = "fpkm")



# Craniopharyngioma
cranio_rna_clinical <- read_tsv("./data/rna/cbttc/craniopharyngioma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, cancer=CANCER_TYPE)

cranio_rna_clinical_samples <- read_tsv("./data/rna/cbttc/craniopharyngioma_cbttc/data_rna_seq_v2_mrna.txt") %>%
  rename(gene = Hugo_Symbol) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  select(sample) %>%
  distinct() %>%
  inner_join(cranio_rna_clinical[, c("Kids.First.ID", "Clinical.Event.Id", "cancer")], by = c("sample" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, cancer) %>%
  tally() %>%
  ungroup()

cranio_rna_clinical_prot <- cranio_rna_clinical %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "cancer" = "diagnosis")) %>%
  arrange(Kids.First.ID)

cranio_rna <- read_tsv("./data/rna/cbttc/craniopharyngioma_cbttc/data_rna_seq_v2_mrna.txt") %>%
  rename(gene = Hugo_Symbol) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  inner_join(cranio_rna_clinical_prot[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  spread(key = "Kids.First.ID", value = "fpkm")



# Ependymoma
epend_rna_clinical <- read_tsv("./data/rna/cbttc/ependymoma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, cancer=CANCER_TYPE)

epend_rna_clinical_samples <- read_tsv("./data/rna/cbttc/ependymoma_cbttc/data_rna_seq_v2_mrna.txt") %>%
  rename(gene = Hugo_Symbol) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  select(sample) %>%
  distinct() %>%
  inner_join(epend_rna_clinical[, c("Kids.First.ID", "Clinical.Event.Id", "cancer")], by = c("sample" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, cancer) %>%
  tally() %>%
  ungroup()

epend_rna_clinical_prot <- epend_rna_clinical %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "cancer" = "diagnosis")) %>%
  arrange(Kids.First.ID)

epend_rna <- read_tsv("./data/rna/cbttc/ependymoma_cbttc/data_rna_seq_v2_mrna.txt") %>%
  rename(gene = Hugo_Symbol) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  inner_join(epend_rna_clinical_prot[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  group_by(Kids.First.ID, gene) %>%
  summarise(fpkm = mean(fpkm)) %>%
  ungroup() %>%
  spread(key = "Kids.First.ID", value = "fpkm")



# Ganglioglioma
gangli_rna_clinical <- read_tsv("./data/rna/cbttc/ganglioglioma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, cancer=CANCER_TYPE)

gangli_rna_clinical_samples <- read_tsv("./data/rna/cbttc/ganglioglioma_cbttc/data_rna_seq_v2_mrna.txt") %>%
  rename(gene = Hugo_Symbol) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  select(sample) %>%
  distinct() %>%
  inner_join(gangli_rna_clinical[, c("Kids.First.ID", "Clinical.Event.Id", "cancer")], by = c("sample" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, cancer) %>%
  tally() %>%
  ungroup()

gangli_rna_clinical_prot <- gangli_rna_clinical %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "cancer" = "diagnosis")) %>%
  arrange(Kids.First.ID)

gangli_rna <- read_tsv("./data/rna/cbttc/ganglioglioma_cbttc/data_rna_seq_v2_mrna.txt") %>%
  rename(gene = Hugo_Symbol) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  inner_join(gangli_rna_clinical_prot[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  spread(key = "Kids.First.ID", value = "fpkm")



# Medulloblastoma
medullo_rna_clinical <- read_tsv("./data/rna/cbttc/medulloblastoma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, cancer=CANCER_TYPE)

medullo_rna_clinical_samples <- read_tsv("./data/rna/cbttc/medulloblastoma_cbttc/data_rna_seq_v2_mrna.txt") %>%
  rename(gene = Hugo_Symbol) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  select(sample) %>%
  distinct() %>%
  inner_join(medullo_rna_clinical[, c("Kids.First.ID", "Clinical.Event.Id", "cancer")], by = c("sample" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, cancer) %>%
  tally() %>%
  ungroup()

medullo_rna_clinical_prot <- medullo_rna_clinical %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "cancer" = "diagnosis")) %>%
  arrange(Kids.First.ID)

medullo_rna <- read_tsv("./data/rna/cbttc/medulloblastoma_cbttc/data_rna_seq_v2_mrna.txt") %>%
  rename(gene = Hugo_Symbol) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  inner_join(medullo_rna_clinical_prot[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  group_by(Kids.First.ID, gene) %>%
  summarise(fpkm = mean(fpkm)) %>%
  ungroup() %>%
  spread(key = "Kids.First.ID", value = "fpkm")



# High-grade glioma
hgg_rna_clinical <- read_tsv("./data/rna/cbttc/ped_high_grade_glioma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, cancer=CANCER_TYPE) %>%
  mutate(cancer = str_replace(cancer, " \\(WHO grade III/IV\\)", ""))

hgg_rna_clinical_samples <- read_tsv("./data/rna/cbttc/ped_high_grade_glioma_cbttc/data_rna_seq_v2_mrna.txt") %>%
  rename(gene = Hugo_Symbol) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  select(sample) %>%
  distinct() %>%
  inner_join(hgg_rna_clinical[, c("Kids.First.ID", "Clinical.Event.Id", "cancer")], by = c("sample" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, cancer) %>%
  tally() %>%
  ungroup()

hgg_rna_clinical_prot <- hgg_rna_clinical %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "cancer" = "diagnosis")) %>%
  arrange(Kids.First.ID)

hgg_rna <- read_tsv("./data/rna/cbttc/ped_high_grade_glioma_cbttc/data_rna_seq_v2_mrna.txt") %>%
  rename(gene = Hugo_Symbol) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  inner_join(hgg_rna_clinical_prot[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  group_by(Kids.First.ID, gene) %>%
  summarise(fpkm = mean(fpkm)) %>%
  ungroup() %>%
  spread(key = "Kids.First.ID", value = "fpkm")



# Low-grade glioma
lgg_rna_clinical <- read_tsv("./data/rna/cbttc/ped_low_grade_glioma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, cancer=CANCER_TYPE) %>%
  mutate(cancer = str_replace(cancer, " \\(WHO grade I/II\\)", ""))

lgg_rna_clinical_samples <- read_tsv("./data/rna/cbttc/ped_low_grade_glioma_cbttc/data_rna_seq_v2_mrna.txt") %>%
  rename(gene = Hugo_Symbol) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  select(sample) %>%
  distinct() %>%
  inner_join(lgg_rna_clinical[, c("Kids.First.ID", "Clinical.Event.Id", "cancer")], by = c("sample" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, cancer) %>%
  tally() %>%
  ungroup()

lgg_rna_clinical_prot <- lgg_rna_clinical %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "cancer" = "diagnosis")) %>%
  arrange(Kids.First.ID)

lgg_rna <- read_tsv("./data/rna/cbttc/ped_low_grade_glioma_cbttc/data_rna_seq_v2_mrna.txt") %>%
  rename(gene = Hugo_Symbol) %>%
  gather(-gene, key = "sample", value = "fpkm") %>%
  inner_join(lgg_rna_clinical_prot[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  group_by(Kids.First.ID, gene) %>%
  summarise(fpkm = mean(fpkm)) %>%
  ungroup() %>%
  spread(key = "Kids.First.ID", value = "fpkm")




# merge all cancer types together
rna_cbttc_fpkm <- inner_join(atrt_rna, cranio_rna, by = "gene") %>%
  inner_join(epend_rna, by="gene") %>%
  inner_join(gangli_rna, by="gene") %>%
  inner_join(medullo_rna, by="gene") %>%
  inner_join(hgg_rna, by="gene") %>%
  inner_join(lgg_rna, by="gene") %>%
  filter(gene %in% protein_coding$gene_name)

gz1 <- gzfile("./output/files/transcriptomics_cbttc_fpkm.txt.gz", "w")
write.table(rna_cbttc_fpkm, gz1, sep="\t", quote=F, row.names=F)
close(gz1)


# cbttc samples
cbttc_samples <- atrt_rna_clinical_samples %>%
  select(-n) %>%
  bind_rows(cranio_rna_clinical_samples[, -c(3)]) %>%
  bind_rows(epend_rna_clinical_samples[, -c(3)]) %>%
  bind_rows(gangli_rna_clinical_samples[, -c(3)]) %>%
  bind_rows(medullo_rna_clinical_samples[, -c(3)]) %>%
  bind_rows(hgg_rna_clinical_samples[, -c(3)]) %>%
  bind_rows(lgg_rna_clinical_samples[, -c(3)]) %>%
  mutate(batch = "cbttc") %>%
  select(sample = Kids.First.ID, batch, cancer)
write.table(cbttc_samples, "./output/files/transcriptomics_samples_cbttc.txt", sep="\t", quote=F, row.names=F)

