library(tidyverse)


# assembling of cnv data from cbttc study
# GISTIC2 scores from pedcbioportal



# cancer samples with proteomics measurements
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt") %>%
  filter(str_detect(batch, "cbttc"))


# clinical data from CBTTC proteomics dataset
cbttc_clinical_prot <- data.table::fread("./data/protein/cptac_new/cbttc/S047_Pediatric_Brain_Cancer_Clinical_Data_r1.txt", check.names = T) %>%
  as_tibble() %>%
  filter(Kids.First.ID %in% proteomics$sample) %>%
  select(Kids.First.ID, Clinical.Event.Id, diagnosis_type, diagnosis)


# patients with two different tumor samples from CBTTC proteomics dataset
cptac_cbttc2tumours <- read_tsv(file = "./output/files/cptac_cbttc2tumours.txt")


# protein coding genes
protein_coding <- read_tsv("./output/files/gencode.v19.annotation.txt") %>%
  select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  select(-gene_type) %>%
  group_by(gene_name) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)


# Atypical Teratoid Rhabdoid Tumor (ATRT)
atrt_meta <- read_tsv("./data/dna/cnv/pedcbioportal/atrt_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, diagnosis=CANCER_TYPE)

atrt_meta_samples <- data.table::fread("./data/dna/cnv/pedcbioportal/atrt_cbttc/data_CNA.txt") %>%
  as_tibble() %>%
  select(gene = Hugo_Symbol, everything()) %>%
  gather(-gene, key = "sample", value = "gistic2") %>%
  select(sample) %>%
  distinct() %>%
  inner_join(atrt_meta, by = c("sample" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, diagnosis) %>%
  tally() %>%
  ungroup()

atrt_meta_samples2 <- atrt_meta %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "diagnosis"))

atrt_cnv <- data.table::fread("./data/dna/cnv/pedcbioportal/atrt_cbttc/data_CNA.txt") %>%
  as_tibble() %>%
  select(gene = Hugo_Symbol, everything()) %>%
  gather(-gene, key = "sample", value = "gistic2") %>%
  inner_join(atrt_meta_samples2[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  select(sample = Kids.First.ID, everything()) %>%
  spread(key = "sample", value = "gistic2")



# Craniopharyngioma
cranio_meta <- read_tsv("./data/dna/cnv/pedcbioportal/craniopharyngioma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, diagnosis=CANCER_TYPE)

cranio_meta_samples <- data.table::fread("./data/dna/cnv/pedcbioportal/craniopharyngioma_cbttc/data_CNA.txt") %>%
  as_tibble() %>%
  select(gene = Hugo_Symbol, everything()) %>%
  gather(-gene, key = "sample", value = "gistic2") %>%
  select(sample) %>%
  distinct() %>%
  inner_join(cranio_meta, by = c("sample" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, diagnosis) %>%
  tally() %>%
  ungroup()

cranio_meta_samples2 <- cranio_meta %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "diagnosis"))

cranio_cnv <- data.table::fread("./data/dna/cnv/pedcbioportal/craniopharyngioma_cbttc/data_CNA.txt") %>%
  as_tibble() %>%
  select(gene = Hugo_Symbol, everything()) %>%
  gather(-gene, key = "sample", value = "gistic2") %>%
  inner_join(cranio_meta_samples2[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  select(sample = Kids.First.ID, everything()) %>%
  spread(key = "sample", value = "gistic2")



# Ependymoma
epend_meta <- read_tsv("./data/dna/cnv/pedcbioportal/ependymoma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, diagnosis=CANCER_TYPE)

epend_meta_samples <- data.table::fread("./data/dna/cnv/pedcbioportal/ependymoma_cbttc/data_CNA.txt") %>%
  as_tibble() %>%
  select(gene = Hugo_Symbol, everything()) %>%
  gather(-gene, key = "sample", value = "gistic2") %>%
  select(sample) %>%
  distinct() %>%
  inner_join(epend_meta, by = c("sample" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, diagnosis) %>%
  tally() %>%
  ungroup()

epend_meta_samples2 <- epend_meta %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "diagnosis")) %>%
  filter(!Clinical.Event.Id.x %in% c("7316-322-T-233349", "7316-365-T-567972"))

epend_cnv <- data.table::fread("./data/dna/cnv/pedcbioportal/ependymoma_cbttc/data_CNA.txt") %>%
  as_tibble() %>%
  select(gene = Hugo_Symbol, everything()) %>%
  gather(-gene, key = "sample", value = "gistic2") %>%
  inner_join(epend_meta_samples2[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  select(sample = Kids.First.ID, everything()) %>%
  spread(key = "sample", value = "gistic2")



# Ganglioglioma
gangli_meta <- read_tsv("./data/dna/cnv/pedcbioportal/ganglioglioma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, diagnosis=CANCER_TYPE)

gangli_meta_samples <- data.table::fread("./data/dna/cnv/pedcbioportal/ganglioglioma_cbttc/data_CNA.txt") %>%
  as_tibble() %>%
  select(gene = Hugo_Symbol, everything()) %>%
  gather(-gene, key = "sample", value = "gistic2") %>%
  select(sample) %>%
  distinct() %>%
  inner_join(gangli_meta, by = c("sample" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, diagnosis) %>%
  tally() %>%
  ungroup()

gangli_meta_samples2 <- gangli_meta %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "diagnosis"))

gangli_cnv <- data.table::fread("./data/dna/cnv/pedcbioportal/ganglioglioma_cbttc/data_CNA.txt") %>%
  as_tibble() %>%
  select(gene = Hugo_Symbol, everything()) %>%
  gather(-gene, key = "sample", value = "gistic2") %>%
  inner_join(gangli_meta_samples2[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  select(sample = Kids.First.ID, everything()) %>%
  spread(key = "sample", value = "gistic2")



# Medulloblastoma
medullo_meta <- read_tsv("./data/dna/cnv/pedcbioportal/medulloblastoma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, diagnosis=CANCER_TYPE)

medullo_meta_samples <- data.table::fread("./data/dna/cnv/pedcbioportal/medulloblastoma_cbttc/data_CNA.txt") %>%
  as_tibble() %>%
  select(gene = Hugo_Symbol, everything()) %>%
  gather(-gene, key = "sample", value = "gistic2") %>%
  select(sample) %>%
  distinct() %>%
  inner_join(medullo_meta, by = c("sample" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, diagnosis) %>%
  tally() %>%
  ungroup()

medullo_meta_samples2 <- medullo_meta %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "diagnosis"))

medullo_cnv <- data.table::fread("./data/dna/cnv/pedcbioportal/medulloblastoma_cbttc/data_CNA.txt") %>%
  as_tibble() %>%
  select(gene = Hugo_Symbol, everything()) %>%
  gather(-gene, key = "sample", value = "gistic2") %>%
  inner_join(medullo_meta_samples2[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  select(sample = Kids.First.ID, everything()) %>%
  spread(key = "sample", value = "gistic2")



# High-grade glioma
hgg_meta <- read_tsv("./data/dna/cnv/pedcbioportal/ped_high_grade_glioma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, diagnosis=CANCER_TYPE) %>%
  mutate(diagnosis = str_replace(diagnosis, " \\(WHO grade III/IV\\)", ""))

hgg_meta_samples <- data.table::fread("./data/dna/cnv/pedcbioportal/ped_high_grade_glioma_cbttc/data_CNA.txt") %>%
  as_tibble() %>%
  select(gene = Hugo_Symbol, everything()) %>%
  gather(-gene, key = "sample", value = "gistic2") %>%
  select(sample) %>%
  distinct() %>%
  inner_join(hgg_meta, by = c("sample" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, diagnosis) %>%
  tally() %>%
  ungroup()

hgg_meta_samples2 <- hgg_meta %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "diagnosis")) %>%
  filter(str_detect(Clinical.Event.Id.x, "\\-CL\\-adh|\\-CL\\-susp", negate = T))


hgg_cnv <- data.table::fread("./data/dna/cnv/pedcbioportal/ped_high_grade_glioma_cbttc/data_CNA.txt") %>%
  as_tibble() %>%
  select(gene = Hugo_Symbol, everything()) %>%
  gather(-gene, key = "sample", value = "gistic2") %>%
  inner_join(hgg_meta_samples2[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  select(sample = Kids.First.ID, everything()) %>%
  spread(key = "sample", value = "gistic2")



# Low-grade glioma
lgg_meta <- read_tsv("./data/dna/cnv/pedcbioportal/ped_low_grade_glioma_cbttc/data_clinical_sample.txt") %>%
  slice(-c(1:4)) %>%
  select(Kids.First.ID=`#Patient Identifier`, Clinical.Event.Id=`Sample Identifier`, diagnosis_type=TUMOR_TYPE, diagnosis=CANCER_TYPE) %>%
  mutate(diagnosis = str_replace(diagnosis, " \\(WHO grade I/II\\)", ""))

lgg_meta_samples <- data.table::fread("./data/dna/cnv/pedcbioportal/ped_low_grade_glioma_cbttc/data_CNA.txt") %>%
  as_tibble() %>%
  select(gene = Hugo_Symbol, everything()) %>%
  gather(-gene, key = "sample", value = "gistic2") %>%
  select(sample) %>%
  distinct() %>%
  inner_join(lgg_meta, by = c("sample" = "Clinical.Event.Id")) %>%
  group_by(Kids.First.ID, diagnosis) %>%
  tally() %>%
  ungroup()

lgg_meta_samples2 <- lgg_meta %>%
  inner_join(cbttc_clinical_prot, by = c("Kids.First.ID", "diagnosis_type", "diagnosis")) %>%
  filter(Clinical.Event.Id.x != "7316-264-T-232408")

lgg_cnv <- data.table::fread("./data/dna/cnv/pedcbioportal/ped_low_grade_glioma_cbttc/data_CNA.txt") %>%
  as_tibble() %>%
  select(gene = Hugo_Symbol, everything()) %>%
  gather(-gene, key = "sample", value = "gistic2") %>%
  inner_join(lgg_meta_samples2[, c("Kids.First.ID", "Clinical.Event.Id.x")], by = c("sample" = "Clinical.Event.Id.x")) %>%
  select(-sample) %>%
  select(sample = Kids.First.ID, everything()) %>%
  spread(key = "sample", value = "gistic2")



# merge all cancer types together
cnv_cbttc <- inner_join(atrt_cnv, cranio_cnv, by = "gene") %>%
  inner_join(epend_cnv, by="gene") %>%
  inner_join(gangli_cnv, by="gene") %>%
  inner_join(medullo_cnv, by="gene") %>%
  inner_join(hgg_cnv, by="gene") %>%
  inner_join(lgg_cnv, by="gene") %>%
  filter(gene %in% protein_coding$gene_name)

gz1 <- gzfile("./output/files/cnv_cbttc.txt.gz", "w")
write.table(cnv_cbttc, gz1, sep="\t", quote=F, row.names=F)
close(gz1)


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
write.table(cbttc_samples, "./output/files/cnv_samples_cbttc.txt", sep="\t", quote=F, row.names=F)

