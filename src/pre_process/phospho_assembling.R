library(tidyverse)
library(data.table)


source("./src/utils/psite_prot_seq_match.R")
source("./src/utils/psites_duplicated.R")


# protein sequences
prot_seqs <- read_tsv(file = "./output/files/protein_sequences_GRCh37.txt") %>%
  filter(state == "known" & gene_biotype == "protein_coding" & trpt_biotype == "protein_coding")


# canonical uniprot transcripts
canonical_trpt <- read_tsv(file = "./data/uniprot/canonical_trpt/isoform_overrides_uniprot.txt", skip = 1)


gene_prot_seq <- canonical_trpt %>%
  select(gene_name, enst_id) %>%
  inner_join(prot_seqs[, c("trpt_id", "seq")], by = c("enst_id" = "trpt_id"))


# load cptac phospho data assembled in "Multi-omics characterization of... human protein abundance levels"
cptac1_phospho <- read_tsv("./data/protein/cptac/cptac_phospho.txt") %>%
  separate(col = phospho_site, into = c("gene", "psite"), sep = "_") %>%
  mutate(psite = toupper(psite)) %>%
  mutate(psite = paste(gene, psite, sep = "_"))

cptac1_prot_samples <- read_tsv("./data/protein/cptac/cptac_samples.txt")

cptac1_phospho_samples <- cptac1_prot_samples %>%
  filter(sample %in% colnames(cptac1_phospho)[-c(1:2)]) %>%
  mutate(cancer = tolower(cancer)) %>%
  mutate(batch = if_else(cancer == "brca", "cptac-tcga-brca", "cptac-tcga-hgsc")) %>%
  select(sample, batch, cancer)



# load cptac phospho data from "TCGA breast cancer"

# qc status samples
brca_failed_samples <- read_csv("./data/protein/cptac/CPTAC_TCGA_BreastCancer_select_clinical_data_r1_QC_status.csv") %>%
  filter(`QC Status` == "fail") %>%
  pull(`TCGA ID`) %>%
  str_replace_all(., "-", ".")

cptac_tcga_breast <- read_tsv("./data/protein/cptac/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq.tsv") %>%
  select(-c("Peptide", "Organism")) %>%
  select(gene = Gene, psite = Phosphosite, everything()) %>%
  mutate(psite = toupper(psite)) %>%
  separate(psite, into = c("prot_id", "psite"), sep = ":") %>%
  mutate(prot_id = str_replace(prot_id, "\\.[0-9]{1}", "")) %>%
  mutate(psites = str_count(psite, "[A-Z]")) %>%
  #select(gene, prot_id, psite, psites) %>%
  select(gene, prot_id, psite, psites, everything()) %>%
  inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("gene" = "gene_name")) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, everything()) %>%
  group_by(n, gene, prot_id) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = map_psite)) %>%
  unnest(cols = data) %>%
  filter(same == 1) %>%
  select(-n, -seq, -same) %>%
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  rename_at(.vars = 5:ncol(.), .funs = ~ str_c("TCGA", .x, sep="-")) %>%
  #group_by(gene, prot_id) %>% tally() %>% ungroup() %>% group_by(gene) %>% mutate(nn = n()) %>% ungroup() %>% filter(nn == 2) %>% arrange(desc(nn))
  mutate(psite = str_c(gene, psite, sep="_")) %>%
  group_by(psite) %>%
  nest() %>%
  ungroup() %>%
  mutate(n = map_dbl(data, nrow))

cptac_tcga_breast_dup <- cptac_tcga_breast %>%
  filter(n>1) %>%
  mutate(same = map_dbl(data, psites_duplicated, indx=4)) %>%
  filter(same == 1) %>%
  select(-n,-same) %>%
  unnest(cols = data) %>%
  filter(!duplicated(psite))

cptac_tcga_breast <- cptac_tcga_breast %>%
  filter(n == 1) %>%
  select(-n) %>%
  unnest(cols = data) %>%
  bind_rows(cptac_tcga_breast_dup) %>%
  select(-prot_id) %>%
  select(gene, psite, everything()) %>%
  gather(key = "sample", value = "tmt", -gene, -psite, -psites) %>%
  mutate(sample = str_split_fixed(sample, "\\.", 2)[,1]) %>%
  group_by(sample, gene, psite, psites) %>%
  summarise(tmt = mean(tmt, na.rm = T)) %>%
  ungroup() %>%
  mutate(tmt = replace(tmt, is.nan(tmt), NA)) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(sample = str_sub(sample, 1, 12)) %>%
  filter(!sample %in% brca_failed_samples) %>%
  spread(key = "sample", value = "tmt")

tcga_breast_samples <- tibble(sample = colnames(cptac_tcga_breast)[-c(1:3)], batch = "cptac-tcga-brca", cancer = "brca")



# load cptac phospho data from "TCGA ovarian cancer"
cptac_tcga_hgsc <- read_tsv("./data/protein/cptac/TCGA_Ovarian_PNNL_Phosphoproteome.phosphosite.itraq.tsv") %>%
  select(-c("Peptide", "Organism")) %>%
  select(gene = Gene, psite = Phosphosite, everything()) %>%
  mutate(psite = toupper(psite)) %>%
  separate(psite, into = c("prot_id", "psite"), sep = ":") %>%
  mutate(prot_id = str_replace(prot_id, "\\.[0-9]{1}", "")) %>%
  mutate(psites = str_count(psite, "[A-Z]")) %>%
  #select(gene, prot_id, psite, psites) %>%
  select(gene, prot_id, psite, psites, everything()) %>%
  inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("gene" = "gene_name")) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, everything()) %>%
  group_by(n, gene, prot_id) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = map_psite)) %>%
  unnest(cols = data) %>%
  filter(same == 1) %>%
  select(-n, -seq, -same) %>%
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  rename_at(.vars = 5:ncol(.), .funs = ~ str_c("TCGA", .x, sep="-")) %>%
  #group_by(gene, prot_id) %>% tally() %>% ungroup() %>% group_by(gene) %>% mutate(nn = n()) %>% ungroup() %>% filter(nn == 2) %>% arrange(desc(nn))
  mutate(psite = str_c(gene, psite, sep="_")) %>%
  group_by(psite) %>%
  nest() %>%
  ungroup() %>%
  mutate(n = map_dbl(data, nrow))

cptac_tcga_hgsc_dup <- cptac_tcga_hgsc %>%
  filter(n>1) %>%
  mutate(same = map_dbl(data, psites_duplicated, indx=4)) %>%
  filter(same == 1) %>%
  select(-n,-same) %>%
  unnest(cols = data) %>%
  filter(!duplicated(psite))

cptac_tcga_hgsc <- cptac_tcga_hgsc %>%
  filter(n == 1) %>%
  select(-n) %>%
  unnest(cols = data) %>%
  bind_rows(cptac_tcga_hgsc_dup) %>%
  select(-prot_id) %>%
  select(gene, psite, everything()) %>%
  rename_at(.vars = 4:ncol(.), .funs = ~ str_replace_all(.x, "-", ".")) %>%
  rename_at(.vars = 4:ncol(.), .funs = ~ str_sub(.x, 1, 12))

tcga_hgsc_samples <- tibble(sample = colnames(cptac_tcga_hgsc)[-c(1:3)], batch = "cptac-tcga-hgsc", cancer = "hgsc")



# load cell lines phospho data assembled in "Multi-omics characterization of... human protein abundance levels"
cells_phospho <- read_tsv("./data/protein/cell_lines_prot/cell_lines_phospho.txt") %>%
  separate(col = phospho_site, into = c("gene", "psite"), sep = "_") %>%
  filter(!(gene == "NA")) %>%
  mutate(psite = toupper(psite)) %>%
  mutate(psite = paste(gene, psite, sep = "_"))


cells_prot_samples <- read_tsv("./data/protein/cell_lines_prot/cell_lines_proteomics_metadata.txt") %>%
  select(-proteomics) %>%
  select(sample = cell_line, batch, cancer = cancer_type) %>%
  mutate(batch = paste("cell-lines", batch, sep = "-")) %>%
  filter(sample %in% colnames(cells_phospho)[-c(1:2)])



# load phospho data from Roumeliotis et al
cell_lines_rmlt <- read_csv("./data/protein/cell_lines_prot/roumeliotis_phosphoprot_50_coread_cell_lines.csv") %>%
  select(gene = `Gene name`, psite = `Protein site`, everything()) %>%
  rename_at(.vars = 3:ncol(.), .funs = ~ toupper(str_replace_all(.x, "-", ""))) %>%
  mutate(psite = toupper(psite)) %>%
  mutate(psites = str_count(psite, "[A-Z]")) %>%
  select(gene, psite, psites, everything()) %>%
  inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("gene" = "gene_name")) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, everything()) %>%
  group_by(n, gene) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = map_psite)) %>%
  unnest(cols = data) %>%
  filter(same == 1) %>%
  select(-n, -seq, -same) %>%
  mutate(psite = str_c(gene, psite, sep="_")) %>%
  group_by(psite) %>%
  nest() %>%
  ungroup() %>%
  mutate(n = map_dbl(data, nrow))

cell_lines_rmlt_dup <- cell_lines_rmlt %>%
  filter(n>1) %>%
  mutate(same = map_dbl(data, psites_duplicated, indx=3)) %>%
  filter(same == 1) %>%
  select(-n,-same) %>%
  unnest(cols = data) %>%
  filter(!duplicated(psite))

cell_lines_rmlt <- cell_lines_rmlt %>%
  filter(n == 1) %>%
  select(-n) %>%
  unnest(cols = data) %>%
  bind_rows(cell_lines_rmlt_dup) %>%
  select(gene, psite, everything()) %>%
  gather(key = "cell", value = "value", -c(gene, psite, psites)) %>%
  # divide each value by 100 (they were multiplied by 100) and transform to log2
  mutate(value = value/100) %>%
  mutate(value = log2(value)) %>%
  spread(key = "cell", value = "value")

cell_lines_rmlt_samples <- tibble(sample = colnames(cell_lines_rmlt)[-c(1:3)], batch = "cell-lines-rmlt", cancer = "coread")



# load cptac phospho data from "Colon Cancer Therapeutic Opportunities"
cptac_colon <- fread("./data/protein/cptac_new/colon/Human__CPTAC_COAD__PNNL__Phosphoproteome__TMT__03_01_2017__BCM__Site__Tumor_PNNL_TMT_LogRatio.cct.gz") %>%
  as_tibble() %>%
  rename(feature = V1) %>%
  separate(col = feature, into = c("gene", "prot_id"), sep = "__") %>%
  mutate(gene = str_split_fixed(gene, "_", 2)[,1]) %>%
  mutate(prot_id = toupper(prot_id)) %>%
  separate(col = prot_id, into = c("prot_id", "psite"), sep = ":") %>%
  mutate(psites = str_count(psite, "[A-Z]")) %>%
  select(gene, prot_id, psite, psites, everything()) %>%
  inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("gene" = "gene_name")) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, everything()) %>%
  group_by(n, gene, prot_id) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = map_psite)) %>%
  unnest(cols = data) %>%
  filter(same == 1) %>%
  select(-n, -seq, -same) %>%
  #group_by(gene, prot_id) %>% tally() %>% ungroup() %>% group_by(gene) %>% mutate(nn = n()) %>% ungroup() %>% filter(nn == 2) %>% arrange(desc(nn))
  mutate(psite = str_c(gene, psite, sep="_")) %>%
  group_by(psite) %>%
  nest() %>%
  ungroup() %>%
  mutate(n = map_dbl(data, nrow))

cptac_colon_dup <- cptac_colon %>%
  filter(n>1) %>%
  mutate(same = map_dbl(data, psites_duplicated, indx=4)) %>%
  filter(same == 1) %>%
  select(-n,-same) %>%
  unnest(cols = data) %>%
  filter(!duplicated(psite))

cptac_colon <- cptac_colon %>%
  filter(n == 1) %>%
  select(-n) %>%
  unnest(cols = data) %>%
  bind_rows(cptac_colon_dup) %>%
  select(-prot_id) %>%
  select(gene, psite, everything())

cptac_colon_samples <- tibble(sample = colnames(cptac_colon)[-c(1:3)], batch = "cptac-colon-opportunities", cancer = "colon-cancer")



# load cptac cell renal cell carcinoma (CCRCC) from "CPTAC CCRCC Discovery Study"
cptac_ccrcc_biosp <- fread("./data/protein/cptac_new/ccrcc/S044_CPTAC_CCRCC_Discovery_Cohort_Specimens_r1_Sept2018.txt", check.names = T) %>%
  as_tibble() %>%
  select(ParticipantID, Aliquot.ID, Group) %>%
  distinct() %>%
  mutate(nchar = map_dbl(ParticipantID, nchar)) %>%
  filter(nchar == 9) %>%
  select(-nchar)

cptac_ccrcc <- read_tsv("./data/protein/cptac_new/ccrcc/CPTAC3_Clear_Cell_Renal_Cell_Carcinoma_Phosphoproteome.phosphosite.tmt10.tsv") %>%
  select(-c("Peptide", "Organism")) %>%
  select(gene = Gene, psite = Phosphosite, everything()) %>%
  mutate(psite = toupper(psite)) %>%
  separate(psite, into = c("prot_id", "psite"), sep = ":") %>%
  mutate(prot_id = str_replace(prot_id, "\\.[0-9]{1}", "")) %>%
  mutate(psites = str_count(psite, "[A-Z]")) %>%
  select(gene, prot_id, psite, psites, everything()) %>%
  inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("gene" = "gene_name")) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, everything()) %>%
  group_by(n, gene, prot_id) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = map_psite)) %>%
  unnest(cols = data) %>%
  filter(same == 1) %>%
  select(-n, -seq, -same) %>%
  select(c("gene", "prot_id", "psite", "psites", starts_with("CPT"))) %>%
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  #group_by(gene, prot_id) %>% tally() %>% ungroup() %>% group_by(gene) %>% mutate(nn = n()) %>% ungroup() %>% filter(nn == 2) %>% arrange(desc(nn))
  mutate(psite = str_c(gene, psite, sep="_")) %>%
  group_by(psite) %>%
  nest() %>%
  ungroup() %>%
  mutate(n = map_dbl(data, nrow))

cptac_ccrcc_dup <- cptac_ccrcc %>%
  filter(n>1) %>%
  mutate(same = map_dbl(data, psites_duplicated, indx=4)) %>%
  filter(same == 1) %>%
  select(-n,-same) %>%
  unnest(cols = data) %>%
  filter(!duplicated(psite))

cptac_ccrcc <- cptac_ccrcc %>%
  filter(n == 1) %>%
  select(-n) %>%
  unnest(cols = data) %>%
  bind_rows(cptac_ccrcc_dup) %>%
  select(-prot_id) %>%
  select(gene, psite, everything()) %>%
  gather(key = "sample", value = "tmt", -gene, -psite, -psites) %>%
  inner_join(cptac_ccrcc_biosp, by = c("sample" = "Aliquot.ID")) %>%
  filter(Group == "Tumor") %>%
  select(-sample, -Group, sample = ParticipantID) %>%
  mutate(sample = str_replace(sample, "-", ".")) %>%
  spread(key = "sample", value = "tmt")

cptac_ccrcc_samples <- tibble(sample = colnames(cptac_ccrcc)[-c(1:3)], batch = "cptac-discovery-ccrcc", cancer = "ccrcc")



# load cptac lung adenocarcinoma (LUAD) from "CPTAC LUAD Discovery Study"
cptac_luad_biosp <- fread("./data/protein/cptac_new/luad/S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.txt", check.names = T) %>%
  as_tibble() %>%
  select(ParticipantID = Participant.ID..case_id., Aliquot.ID = Aliquot..Specimen.Label., Group = Type) %>%
  distinct() %>%
  filter(!ParticipantID %in% c("Normal Only IR", "Taiwanese IR", "Tumor Only IR"))

cptac_luad <- read_tsv("./data/protein/cptac_new/luad/CPTAC3_Lung_Adeno_Carcinoma_Phosphoproteome.phosphosite.tmt10.tsv") %>%
  select(-c("Peptide", "Organism")) %>%
  select(gene = Gene, psite = Phosphosite, everything()) %>%
  mutate(psite = toupper(psite)) %>%
  separate(psite, into = c("prot_id", "psite"), sep = ":") %>%
  mutate(prot_id = str_replace(prot_id, "\\.[0-9]{1}", "")) %>%
  mutate(psites = str_count(psite, "[A-Z]")) %>%
  select(gene, prot_id, psite, psites, everything()) %>%
  inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("gene" = "gene_name")) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, everything()) %>%
  group_by(n, gene, prot_id) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = map_psite)) %>%
  unnest(cols = data) %>%
  filter(same == 1) %>%
  select(-n, -seq, -same) %>%
  select(c("gene", "prot_id", "psite", "psites", starts_with("CPT"))) %>%
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  #group_by(gene, prot_id) %>% tally() %>% ungroup() %>% group_by(gene) %>% mutate(nn = n()) %>% ungroup() %>% filter(nn == 2) %>% arrange(desc(nn))
  mutate(psite = str_c(gene, psite, sep="_")) %>%
  group_by(psite) %>%
  nest() %>%
  ungroup() %>%
  mutate(n = map_dbl(data, nrow))

cptac_luad_dup <- cptac_luad %>%
  filter(n>1) %>%
  mutate(same = map_dbl(data, psites_duplicated, indx=4)) %>%
  filter(same == 1) %>%
  select(-n,-same) %>%
  unnest(cols = data) %>%
  filter(!duplicated(psite))

cptac_luad <- cptac_luad %>%
  filter(n == 1) %>%
  select(-n) %>%
  unnest(cols = data) %>%
  bind_rows(cptac_luad_dup) %>%
  select(-prot_id) %>%
  select(gene, psite, everything()) %>%
  gather(key = "sample", value = "tmt", -gene, -psite, -psites) %>%
  mutate(sample = str_split_fixed(sample, "\\.", n=2)[,1]) %>%
  group_by(sample, gene, psite, psites) %>%
  summarise(tmt = mean(tmt, na.rm = T)) %>%
  ungroup() %>%
  mutate(tmt = replace(tmt, is.nan(tmt), NA)) %>%
  inner_join(cptac_luad_biosp, by = c("sample" = "Aliquot.ID")) %>%
  filter(Group == "Tumor") %>%
  select(-sample, -Group, sample = ParticipantID) %>%
  mutate(sample = str_replace(sample, "-", ".")) %>%
  spread(key = "sample", value = "tmt")

cptac_luad_samples <- tibble(sample = colnames(cptac_luad)[-c(1:3)], batch = "cptac-discovery-luad", cancer = "luad")



# load cptac uterine corpus endometrial carcinoma (UCEC) from "CPTAC UCEC Discovery Study"
cptac_ucec_biosp <- fread("./data/protein/cptac_new/ucec/S043_CPTAC_UCEC_Discovery_Cohort_Study_Specimens_r1_Sept2018.txt", check.names = T) %>%
  as_tibble() %>%
  select(ParticipantID = ParticipantID..Case_ID., Aliquot.ID = Aliquot.ID, Group) %>%
  distinct() %>%
  filter(!(ParticipantID == "" & Aliquot.ID == "" & Group == "")) %>%
  filter(!ParticipantID == "Ref") %>%
  filter(!Group == "Withdrawn") %>%
  filter(!str_detect(ParticipantID, "NX"))

cptac_ucec <- read_tsv("./data/protein/cptac_new/ucec/CPTAC3_Uterine_Corpus_Endometrial_Carcinoma_Phosphoproteome.phosphosite.tmt10.tsv") %>%
  select(-c("Peptide", "Organism")) %>%
  select(gene = Gene, psite = Phosphosite, everything()) %>%
  mutate(psite = toupper(psite)) %>%
  separate(psite, into = c("prot_id", "psite"), sep = ":") %>%
  mutate(prot_id = str_replace(prot_id, "\\.[0-9]{1}", "")) %>%
  mutate(psites = str_count(psite, "[A-Z]")) %>%
  select(gene, prot_id, psite, psites, everything()) %>%
  inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("gene" = "gene_name")) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, everything()) %>%
  group_by(n, gene, prot_id) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = map_psite)) %>%
  unnest(cols = data) %>%
  filter(same == 1) %>%
  select(-n, -seq, -same) %>%
  select(c("gene", "prot_id", "psite", "psites", starts_with("CPT"))) %>%
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  #group_by(gene, prot_id) %>% tally() %>% ungroup() %>% group_by(gene) %>% mutate(nn = n()) %>% ungroup() %>% filter(nn == 2) %>% arrange(desc(nn))
  mutate(psite = str_c(gene, psite, sep="_")) %>%
  group_by(psite) %>%
  nest() %>%
  ungroup() %>%
  mutate(n = map_dbl(data, nrow))

cptac_ucec_dup <- cptac_ucec %>%
  filter(n>1) %>%
  mutate(same = map_dbl(data, psites_duplicated, indx=4)) %>%
  filter(same == 1) %>%
  select(-n,-same) %>%
  unnest(cols = data) %>%
  filter(!duplicated(psite))

cptac_ucec <- cptac_ucec %>%
  filter(n == 1) %>%
  select(-n) %>%
  unnest(cols = data) %>%
  bind_rows(cptac_ucec_dup) %>%
  select(-prot_id) %>%
  select(gene, psite, everything()) %>%
  gather(key = "sample", value = "tmt", -gene, -psite, -psites) %>%
  inner_join(cptac_ucec_biosp, by = c("sample" = "Aliquot.ID")) %>%
  filter(Group == "Tumor") %>%
  select(-sample, -Group, sample = ParticipantID) %>%
  mutate(sample = str_replace(sample, "-", ".")) %>%
  spread(key = "sample", value = "tmt")

cptac_ucec_samples <- tibble(sample = colnames(cptac_ucec)[-c(1:3)], batch = "cptac-discovery-ucec", cancer = "ucec")



# load cptac early-onset gastric cancer (EOGC) from "Proteogenomic Characterization of Human Early-Onset Gastric Cancer"
cptac_eogc <- read_tsv("./data/protein/cptac_new/eogc/Gastric_Cancer_Korea_Phosphoproteome.phosphosite.itraq.tsv") %>%
  select(-c("Peptide", "Organism")) %>%
  select(gene = Gene, psite = Phosphosite, everything()) %>%
  mutate(psite = toupper(psite)) %>%
  separate(psite, into = c("prot_id", "psite"), sep = ":") %>%
  mutate(prot_id = str_replace(prot_id, "\\.[0-9]{1}", "")) %>%
  mutate(psites = str_count(psite, "[A-Z]")) %>%
  select(gene, prot_id, psite, psites, everything()) %>%
  inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("gene" = "gene_name")) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, everything()) %>%
  group_by(n, gene, prot_id) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = map_psite)) %>%
  unnest(cols = data) %>%
  filter(same == 1) %>%
  select(-n, -seq, -same) %>%
  rename_at(.vars = 5:ncol(.), .funs = ~ str_split_fixed(.x, " ", n=2)[,2]) %>%
  #group_by(gene, prot_id) %>% tally() %>% ungroup() %>% group_by(gene) %>% mutate(nn = n()) %>% ungroup() %>% filter(nn == 2) %>% arrange(desc(nn))
  mutate(psite = str_c(gene, psite, sep="_")) %>%
  group_by(psite) %>%
  nest() %>%
  ungroup() %>%
  mutate(n = map_dbl(data, nrow))

cptac_eogc_dup <- cptac_eogc %>%
  filter(n>1) %>%
  mutate(same = map_dbl(data, psites_duplicated, indx=4)) %>%
  filter(same == 1) %>%
  select(-n,-same) %>%
  unnest(cols = data) %>%
  filter(!duplicated(psite))

cptac_eogc <- cptac_eogc %>%
  filter(n == 1) %>%
  select(-n) %>%
  unnest(cols = data) %>%
  bind_rows(cptac_eogc_dup) %>%
  select(-prot_id) %>%
  select(gene, psite, everything()) %>%
  gather(key = "sample", value = "tmt", -gene, -psite, -psites) %>%
  separate(col = "sample", into = c("tumour", "normal"), sep = "/") %>%
  mutate(tumour = str_split_fixed(tumour, "\\.", n=2)[,1], normal = str_split_fixed(normal, "\\.", n=2)[,1]) %>%
  unite(col = "sample", tumour, normal, sep = "") %>%
  #group_by(sample) %>% tally() %>% arrange(desc(n))
  group_by(sample, gene, psite, psites) %>%
  summarise(tmt = mean(tmt, na.rm = T)) %>%
  ungroup() %>%
  mutate(tmt = replace(tmt, is.nan(tmt), NA)) %>%
  spread(key = "sample", value = "tmt")

cptac_eogc_samples <- tibble(sample = colnames(cptac_eogc)[-c(1:3)], batch = "cptac-GC", cancer = "eogc")



# load cptac children’s brain tumor tissue consortium (CBTTC) from "Proteogenomic Analysis of Pediatric Brain Cancer Tumors Pilot Study"
cptac_cbttc_biosp <- fread("./data/protein/cptac_new/cbttc/S047_Pediatric_Brain_Cancer_Clinical_Data_r1.txt", check.names = T) %>%
  as_tibble() %>%
  select(research_id, Kids.First.ID, Clinical.Event.Id, diagnosis_type, diagnosis)

cptac_cbttc <- read_tsv("./data/protein/cptac_new/cbttc/Gygi_TCMP_HMS_Phosphoproteome.phosphosite.tmt11.tsv") %>%
  select(-c("Peptide", "Organism")) %>%
  select(gene = Gene, psite = Phosphosite, everything()) %>%
  mutate(psite = toupper(psite)) %>%
  separate(psite, into = c("prot_id", "psite"), sep = ":") %>%
  mutate(prot_id = str_replace(prot_id, "\\.[0-9]{1}", "")) %>%
  mutate(psites = str_count(psite, "[A-Z]")) %>%
  select(gene, prot_id, psite, psites, everything()) %>%
  inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("gene" = "gene_name")) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, everything()) %>%
  group_by(n, gene, prot_id) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = map_psite)) %>%
  unnest(cols = data) %>%
  filter(same == 1) %>%
  select(-n, -seq, -same) %>%
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  #group_by(gene, prot_id) %>% tally() %>% ungroup() %>% group_by(gene) %>% mutate(nn = n()) %>% ungroup() %>% filter(nn == 2) %>% arrange(desc(nn))
  mutate(psite = str_c(gene, psite, sep="_")) %>%
  group_by(psite) %>%
  nest() %>%
  ungroup() %>%
  mutate(n = map_dbl(data, nrow))

cptac_cbttc_dup <- cptac_cbttc %>%
  filter(n>1) %>%
  mutate(same = map_dbl(data, psites_duplicated, indx=4)) %>%
  filter(same == 1) %>%
  select(-n,-same) %>%
  unnest(cols = data) %>%
  filter(!duplicated(psite))

cptac_cbttc <- cptac_cbttc %>%
  filter(n == 1) %>%
  select(-n) %>%
  unnest(cols = data) %>%
  bind_rows(cptac_cbttc_dup) %>%
  select(-prot_id) %>%
  select(gene, psite, everything()) %>%
  gather(key = "sample", value = "tmt", -gene, -psite, -psites) %>%
  mutate(sample = str_split_fixed(sample, "\\.", n=2)[,1]) %>%
  group_by(sample, gene, psite, psites) %>%
  summarise(tmt = mean(tmt, na.rm = T)) %>%
  ungroup() %>%
  mutate(tmt = replace(tmt, is.nan(tmt), NA)) %>%
  inner_join(cptac_cbttc_biosp, by = c("sample" = "Clinical.Event.Id"))

cptac_cbttc2tumours <- cptac_cbttc_biosp %>%
  filter( Kids.First.ID  %in% c(cptac_cbttc %>% group_by(research_id, Kids.First.ID) %>% tally() %>% ungroup() %>% arrange(desc(n)) %>% slice(1:20) %>% pull(Kids.First.ID)) )
write.table(cptac_cbttc2tumours, file = "./output/files/cptac_cbttc2tumours.txt", sep = "\t", quote = F, row.names = F)

cptac_cbttc <- cptac_cbttc %>%
  select(sample = Kids.First.ID, gene, psite, psites, tmt) %>%
  filter(!(sample %in% unique(cptac_cbttc2tumours$Kids.First.ID))) %>%
  spread(key = "sample", value = "tmt")

cptac_cbttc_samples <- tibble(sample = colnames(cptac_cbttc)[-c(1:3)], batch = "cptac-cbttc") %>%
  inner_join(cptac_cbttc_biosp %>% select(sample = Kids.First.ID, cancer = diagnosis), by = "sample")



# load cptac HBV-HCC proteomics data from "Integrated Proteogenomic Characterization of HBV-Related Hepatocellular Carcinoma"
cptac_hcc_biosp <- fread("./data/protein/cptac_new/hbv-hcc/S049_HCC_Clinical_Information_and_TMT11_Sample_Mapping_Gao2019.r1.txt", check.names = T) %>%
  as_tibble() %>%
  mutate(Tumor.ID.in.Proteome.experiment = str_c("T", Tumor.ID.in.Proteome.experiment))

cptac_hcc <- read_tsv("./data/protein/cptac_new/hbv-hcc/Zhou_Liver_Cancer_Phosphoproteome.phosphosite.tmt11.tsv") %>%
  select(-c("Peptide", "Organism")) %>%
  select(gene = Gene, psite = Phosphosite, everything()) %>%
  mutate(psite = toupper(psite)) %>%
  separate(psite, into = c("prot_id", "psite"), sep = ":") %>%
  mutate(prot_id = str_replace(prot_id, "\\.[0-9]{1}", "")) %>%
  mutate(psites = str_count(psite, "[A-Z]")) %>%
  select(gene, prot_id, psite, psites, everything()) %>%
  inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("gene" = "gene_name")) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, everything()) %>%
  group_by(n, gene, prot_id) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = map_psite)) %>%
  unnest(cols = data) %>%
  filter(same == 1) %>%
  select(-n, -seq, -same) %>%
  select(c("gene", "prot_id", "psite", "psites", starts_with("T"))) %>%
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  select_if(.predicate = colnames(.) %in% c("gene", "prot_id", "psite", "psites", cptac_hcc_biosp$Tumor.ID.in.Proteome.experiment)) %>%
  #group_by(gene, prot_id) %>% tally() %>% ungroup() %>% group_by(gene) %>% mutate(nn = n()) %>% ungroup() %>% filter(nn == 2) %>% arrange(desc(nn))
  mutate(psite = str_c(gene, psite, sep="_")) %>%
  group_by(psite) %>%
  nest() %>%
  ungroup() %>%
  mutate(n = map_dbl(data, nrow))

cptac_hcc_dup <- cptac_hcc %>%
  filter(n>1) %>%
  mutate(same = map_dbl(data, psites_duplicated, indx=4)) %>%
  filter(same == 1) %>%
  select(-n,-same) %>%
  unnest(cols = data) %>%
  filter(!duplicated(psite))

cptac_hcc <- cptac_hcc %>%
  filter(n == 1) %>%
  select(-n) %>%
  unnest(cols = data) %>%
  bind_rows(cptac_hcc_dup) %>%
  select(-prot_id) %>%
  select(gene, psite, everything())

cptac_hcc_samples <- tibble(sample = colnames(cptac_hcc)[-c(1:3)], batch = "cptac-hcc", cancer = "hcc")



# join sample batch and cancer together
proteomics_samples <- bind_rows(
  tcga_breast_samples,
  tcga_hgsc_samples,
  cell_lines_rmlt_samples,
  cptac_colon_samples,
  cptac_ccrcc_samples,
  cptac_luad_samples,
  cptac_ucec_samples,
  cptac_eogc_samples,
  cptac_cbttc_samples,
  cptac_hcc_samples
)

write.table(proteomics_samples, "./output/files/phosphoproteomics_samples.txt", sep="\t", quote=F, row.names=F)



# join all datasets together by gene and psite
# keep all genes and psites
# remove genes/psites with only NAs
phospho <- full_join(cptac_tcga_breast, cptac_tcga_hgsc, by = c("gene", "psite", "psites")) %>%
  full_join(cell_lines_rmlt, by = c("gene", "psite", "psites")) %>%
  full_join(cptac_colon, by = c("gene", "psite", "psites")) %>%
  full_join(cptac_ccrcc, by = c("gene", "psite", "psites")) %>%
  full_join(cptac_luad, by = c("gene", "psite", "psites")) %>%
  full_join(cptac_ucec, by = c("gene", "psite", "psites")) %>%
  full_join(cptac_eogc, by = c("gene", "psite", "psites")) %>%
  full_join(cptac_cbttc, by = c("gene", "psite", "psites")) %>%
  full_join(cptac_hcc, by = c("gene", "psite", "psites")) %>%
  gather(key = "sample", value = "log2FC", -gene, -psite, -psites) %>%
  group_by(gene, psite, psites) %>%
  #filter( !(sum(is.na(log2FC)) == length(unique(sample))) ) %>%
  filter( !(sum(is.na(log2FC)) == n()) ) %>%
  ungroup() %>%
  spread(key = "sample", value = "log2FC")

gz1 <- gzfile("./output/files/phosphoproteomics.txt.gz", "w")
write.table(phospho, gz1, sep="\t", quote=F, row.names=F)
close(gz1)
