library(tidyverse)
library(data.table)



# load cptac proteomics data assembled in "Multi-omics characterization of... human protein abundance levels"
cptac1_prot <- read_tsv("./data/protein/cptac/cptac_proteomics.txt")

cptac1_prot_samples <- read_tsv("./data/protein/cptac/cptac_samples.txt") %>%
  mutate(cancer = tolower(cancer)) %>%
  mutate(batch = if_else(cancer == "brca", "cptac-tcga-brca", if_else(cancer == "coread", "cptac-tcga-coread", "cptac-tcga-hgsc"))) %>%
  select(sample, batch, cancer)



# load cell lines proteomics data assembled in "Multi-omics characterization of... human protein abundance levels"
cells_prot <- read_tsv("./data/protein/cell_lines_prot/cell_lines_proteomics.txt")

cells_prot_samples <- read_tsv("./data/protein/cell_lines_prot/cell_lines_proteomics_metadata.txt") %>%
  select(-proteomics) %>%
  select(sample = cell_line, batch, cancer = cancer_type) %>%
  mutate(batch = paste("cell-lines", batch, sep = "-"))




# load cptac proteomics data from "Colon Cancer Therapeutic Opportunities"
cptac_colon <- read_tsv("./data/protein/cptac_new/colon/Human__CPTAC_COAD__PNNL__Proteome__TMT__03_01_2017__BCM__Gene__PNNL_Tumor_TMT_UnsharedLogRatio.cct") %>%
  rename(gene = attrib_name)

cptac_colon_samples <- tibble(sample = colnames(cptac_colon)[-c(1)], batch = "cptac-colon-opportunities", cancer = "coread")


# load cptac cell renal cell carcinoma (CCRCC) from "CPTAC CCRCC Discovery Study"
cptac_ccrcc_biosp <- fread("./data/protein/cptac_new/ccrcc/S044_CPTAC_CCRCC_Discovery_Cohort_Specimens_r1_Sept2018.txt", check.names = T) %>%
  as_tibble() %>%
  select(ParticipantID, Aliquot.ID, Group) %>%
  distinct() %>%
  mutate(nchar = map_dbl(ParticipantID, nchar)) %>%
  filter(nchar == 9) %>%
  select(-nchar)

cptac_ccrcc <- read_tsv("./data/protein/cptac_new/ccrcc/CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Proteome.tmt10.tsv") %>%
  slice(-c(1:3)) %>%
  rename(gene = Gene) %>%
  select(-c("NCBIGeneID", "Authority", "Description", "Organism", "Chromosome", "Locus")) %>%
  select(c("gene", contains("Unshared"))) %>%
  select(c("gene", starts_with("CPT"))) %>%
  #rename_all(.funs = function(x) str_split_fixed(x, " ", n=2)[,1])
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  gather(key = "sample", value = "tmt", -gene) %>%
  inner_join(cptac_ccrcc_biosp, by = c("sample" = "Aliquot.ID")) %>%
  filter(Group == "Tumor") %>%
  select(-sample, -Group, sample = ParticipantID) %>%
  mutate(sample = str_replace(sample, "-", ".")) %>%
  spread(key = "sample", value = "tmt")

cptac_ccrcc_samples <- tibble(sample = colnames(cptac_ccrcc)[-c(1)], batch = "cptac-discovery-ccrcc", cancer = "ccrcc")



# load cptac lung adenocarcinoma (LUAD) from "CPTAC LUAD Discovery Study"
cptac_luad_biosp <- fread("./data/protein/cptac_new/luad/S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.txt", check.names = T) %>%
  as_tibble() %>%
  select(ParticipantID = Participant.ID..case_id., Aliquot.ID = Aliquot..Specimen.Label., Group = Type) %>%
  distinct() %>%
  filter(!ParticipantID %in% c("Normal Only IR", "Taiwanese IR", "Tumor Only IR"))

cptac_luad <- read_tsv("./data/protein/cptac_new/luad/CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv") %>%
  slice(-c(1:3)) %>%
  rename(gene = Gene) %>%
  select(-c("NCBIGeneID", "Authority", "Description", "Organism", "Chromosome", "Locus")) %>%
  select(c("gene", contains("Unshared"))) %>%
  select(c("gene", starts_with("CPT"))) %>%
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  gather(key = "sample", value = "tmt", -gene) %>%
  mutate(sample = str_split_fixed(sample, "\\.", n=2)[,1]) %>%
  group_by(sample, gene) %>%
  summarise(tmt = mean(tmt, na.rm = T)) %>%
  ungroup() %>%
  inner_join(cptac_luad_biosp, by = c("sample" = "Aliquot.ID")) %>%
  filter(Group == "Tumor") %>%
  select(-sample, -Group, sample = ParticipantID) %>%
  mutate(sample = str_replace(sample, "-", ".")) %>%
  spread(key = "sample", value = "tmt")

cptac_luad_samples <- tibble(sample = colnames(cptac_luad)[-c(1)], batch = "cptac-discovery-luad", cancer = "luad")



# load cptac uterine corpus endometrial carcinoma (UCEC) from "CPTAC UCEC Discovery Study"
cptac_ucec_biosp <- fread("./data/protein/cptac_new/ucec/S043_CPTAC_UCEC_Discovery_Cohort_Study_Specimens_r1_Sept2018.txt", check.names = T) %>%
  as_tibble() %>%
  select(ParticipantID = ParticipantID..Case_ID., Aliquot.ID = Aliquot.ID, Group) %>%
  distinct() %>%
  filter(!(ParticipantID == "" & Aliquot.ID == "" & Group == "")) %>%
  filter(!ParticipantID == "Ref") %>%
  filter(!Group == "Withdrawn") %>%
  filter(!str_detect(ParticipantID, "NX"))

cptac_ucec <- read_tsv("./data/protein/cptac_new/ucec/CPTAC3_Uterine_Corpus_Endometrial_Carcinoma_Proteome.tmt10.tsv") %>%
  slice(-c(1:3)) %>%
  rename(gene = Gene) %>%
  select(-c("NCBIGeneID", "Authority", "Description", "Organism", "Chromosome", "Locus")) %>%
  select(c("gene", contains("Unshared"))) %>%
  select(c("gene", starts_with("CPT"))) %>%
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  gather(key = "sample", value = "tmt", -gene) %>%
  inner_join(cptac_ucec_biosp, by = c("sample" = "Aliquot.ID")) %>%
  filter(Group == "Tumor") %>%
  select(-sample, -Group, sample = ParticipantID) %>%
  mutate(sample = str_replace(sample, "-", ".")) %>%
  spread(key = "sample", value = "tmt")

cptac_ucec_samples <- tibble(sample = colnames(cptac_ucec)[-c(1)], batch = "cptac-discovery-ucec", cancer = "ucec")



# load cptac early-onset gastric cancer (EOGC) from "Proteogenomic Characterization of Human Early-Onset Gastric Cancer"
cptac_eogc <- read_tsv("./data/protein/cptac_new/eogc/Gastric_Cancer_Korea_Proteome.itraq.tsv") %>%
  slice(-c(1:3)) %>%
  rename(gene = Gene) %>%
  select(-c("NCBIGeneID", "Authority", "Description", "Organism", "Chromosome", "Locus")) %>%
  select(c("gene", contains("Unshared"))) %>%
  #rename_if(.predicate = (function(x) x != "gene"), .funs = ~ str_split_fixed(.x, " ", n=3)[,3]) %>%
  #rename_if(.predicate = ~ .x != "gene", .funs = ~ str_split_fixed(.x, " ", n=3)[,3]) %>%
  rename_at(.vars = 2:length(.), .funs = ~ str_split_fixed(.x, " ", n=3)[,3]) %>%
  gather(key = "sample", value = "itraq", -gene) %>%
  separate(col = "sample", into = c("tumour", "normal"), sep = "/") %>%
  mutate(tumour = str_split_fixed(tumour, "\\.", n=2)[,1], normal = str_split_fixed(normal, "\\.", n=2)[,1]) %>%
  unite(col = "sample", tumour, normal, sep = "") %>%
  #group_by(sample) %>% count() %>% as.data.frame() %>% arrange(n)
  group_by(sample, gene)%>%
  summarise(itraq = mean(itraq, na.rm = T)) %>%
  ungroup() %>%
  spread(key = "sample", value = "itraq")

cptac_eogc_samples <- tibble(sample = colnames(cptac_eogc)[-c(1)], batch = "cptac-GC", cancer = "eogc")



# load cptac children’s brain tumor tissue consortium (CBTTC) from "Proteogenomic Analysis of Pediatric Brain Cancer Tumors Pilot Study"
cptac_cbttc_biosp <- fread("./data/protein/cptac_new/cbttc/S047_Pediatric_Brain_Cancer_Clinical_Data_r1.txt", check.names = T) %>%
  as_tibble() %>%
  select(research_id, Kids.First.ID, Clinical.Event.Id, diagnosis_type, diagnosis)

cptac_cbttc <- read_tsv("./data/protein/cptac_new/cbttc/Gygi_TCMP_HMS_Proteome.tmt11.tsv") %>%
  slice(-c(1:3)) %>%
  rename(gene = Gene) %>%
  select(-c("NCBIGeneID", "Authority", "Description", "Organism", "Chromosome", "Locus")) %>%
  select(c("gene", contains("Unshared"))) %>%
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  gather(key = "sample", value = "tmt", -gene) %>%
  mutate(sample = str_split_fixed(sample, "\\.", n=2)[,1]) %>%
  group_by(sample, gene) %>%
  summarise(tmt = mean(tmt, na.rm = T)) %>%
  ungroup() %>%
  inner_join(cptac_cbttc_biosp, by = c("sample" = "Clinical.Event.Id"))

cptac_cbttc2tumours <- cptac_cbttc_biosp %>%
  filter( research_id  %in% (cptac_cbttc %>% group_by(research_id, Kids.First.ID) %>% count() %>% arrange(n) %>% as.data.frame() %>% filter(n > 10000) %>% pull(research_id)) )

cptac_cbttc <- cptac_cbttc %>%
  #group_by(Kids.First.ID, research_id, gene) %>%
  #summarise(tmt = mean(tmt, na.rm = T)) %>%
  #ungroup() %>%
  select(Kids.First.ID, research_id, gene, tmt) %>%
  filter(!(Kids.First.ID %in% unique(cptac_cbttc2tumours$Kids.First.ID))) %>%
  select(-research_id, sample = Kids.First.ID) %>%
  spread(key = "sample", value = "tmt")


cptac_cbttc_samples <- tibble(sample = colnames(cptac_cbttc)[-c(1)], batch = "cptac-cbttc") %>%
  inner_join(cptac_cbttc_biosp %>% select(sample = Kids.First.ID, cancer = diagnosis), by = "sample")



# load cptac HBV-HCC proteomics data from "Integrated Proteogenomic Characterization of HBV-Related Hepatocellular Carcinoma"
cptac_hcc_biosp <- fread("./data/protein/cptac_new/hbv-hcc/S049_HCC_Clinical_Information_and_TMT11_Sample_Mapping_Gao2019.r1.txt", check.names = T) %>%
  as_tibble() %>%
  mutate(Tumor.ID.in.Proteome.experiment = str_c("T", Tumor.ID.in.Proteome.experiment))

cptac_hcc <- read_tsv("./data/protein/cptac_new/hbv-hcc/Zhou_Liver_Cancer_Proteome.tmt11.tsv") %>%
  slice(-c(1:3)) %>%
  rename(gene = Gene) %>%
  select(-c("NCBIGeneID", "Authority", "Description", "Organism", "Chromosome", "Locus")) %>%
  select(c("gene", contains("Unshared"))) %>%
  select(c("gene", starts_with("T"))) %>%
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  #select(c("gene", colnames(.)[colnames(.) %in% cptac_hcc_biosp$Tumor.ID.in.Proteome.experiment]))
  select_if(.predicate = colnames(.) %in% c("gene", cptac_hcc_biosp$Tumor.ID.in.Proteome.experiment))

cptac_hcc_samples <- tibble(sample = colnames(cptac_hcc)[-c(1)], batch = "cptac-hcc", cancer = "hcc")



# join sample batch and cancer together
proteomics_samples <- bind_rows(
  cptac1_prot_samples,
  cells_prot_samples,
  cptac_colon_samples,
  cptac_ccrcc_samples,
  cptac_luad_samples,
  cptac_ucec_samples,
  cptac_eogc_samples,
  cptac_cbttc_samples,
  cptac_hcc_samples
)

write.table(proteomics_samples, "./output/files/proteomics_samples.txt", sep="\t", quote=F, row.names=F)



# join all datasets together by gene
# keep all genes
# remove genes with only NAs
proteomics <- full_join(cptac1_prot, cells_prot, by = "gene") %>%
  full_join(cptac_colon, by = "gene") %>%
  full_join(cptac_ccrcc, by = "gene") %>%
  full_join(cptac_luad, by = "gene") %>%
  full_join(cptac_ucec, by = "gene") %>%
  full_join(cptac_eogc, by = "gene") %>%
  full_join(cptac_cbttc, by = "gene") %>%
  full_join(cptac_hcc, by = "gene") %>%
  gather(key = "sample", value = "log2FC", -gene) %>%
  group_by(gene) %>%
  #filter( !(sum(is.na(log2FC)) == length(unique(sample))) ) %>%
  filter( !(sum(is.na(log2FC)) == n()) ) %>%
  ungroup() %>%
  spread(key = "sample", value = "log2FC")

gz1 <- gzfile("./output/files/proteomics.txt.gz", "w")
write.table(proteomics, gz1, sep="\t", quote=F, row.names=F)
close(gz1)
