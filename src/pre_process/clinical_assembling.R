library(tidyverse)
library(TCGAbiolinks)


# assembling of clinical data



# cancer samples with proteomics measurements
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt")


# get TCGA cancer subtypes for breast, ovarian and colorectal tissues
tcga_subtypes <- PanCancerAtlas_subtypes() 

tcga_subtypes1 <- tcga_subtypes %>%
  filter(cancer.type %in% c("BRCA", "OVCA")) %>%
  filter(str_sub(pan.samplesID, 14, 15) == "01") %>%
  mutate(pan.samplesID = str_sub(pan.samplesID, 1, 12)) %>%
  distinct() %>%
  mutate(pan.samplesID = str_replace_all(pan.samplesID, "-", "."))

# colorectal classification comes from this publication - https://www.cell.com/cancer-cell/pdf/S1535-6108(18)30114-4.pdf
tcga_subtypes2 <- tcga_subtypes %>%
  filter(cancer.type %in% c("COAD", "READ")) %>%
  distinct() %>%
  mutate(pan.samplesID = str_replace_all(pan.samplesID, "-", "."))

# cancer subtypes from the original TCGA papers
brca <- TCGAquery_subtype("brca") %>%
  as_tibble() %>%
  select(sample = patient, cancer_subtype = BRCA_Subtype_PAM50) %>%
  mutate(sample = str_replace_all(sample, "-", "\\."))

coread <- TCGAquery_subtype("coad") %>%
  as_tibble() %>%
  select(sample = patient, methylation_subtype, expression_subtype) %>%
  mutate(sample = str_replace_all(sample, "-", "\\."))


# survival data from "An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics"
# https://www.sciencedirect.com/science/article/pii/S0092867418302290?via%3Dihub#app2
tcga_survival <- read_tsv(file = "./data/clinical/tcga_publications/survival_data/tcga_clinical_data.txt") %>%
  mutate(bcr_patient_barcode = str_replace_all(bcr_patient_barcode, "-", "\\."))



# TCGA brca
# cptac data
tcga_brca <- read_csv(file = "./data/clinical/cptac/brca/CPTAC_TCGA_BreastCancer_select_clinical_data_r1.csv") %>%
  select(sample=`Hybridization REF`, age=yearstobirth, gender, histological_type = histologicaltype, tumor_stage=neoplasm.diseasestage) %>%
  mutate(sample = toupper(str_replace_all(sample, "-", "\\."))) %>%
  distinct()

cbio_brca <- read_tsv(file = "./data/clinical/cbioportal/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", skip = 4) %>%
  select(sample=PATIENT_ID, cancer=CANCER_TYPE_ACRONYM, cancer_subtype=SUBTYPE, age=AGE, gender=GENDER, ajcc_tumor_stage=AJCC_PATHOLOGIC_TUMOR_STAGE, ethnicity=ETHNICITY, race=RACE, seq_center=CENTER) %>%
  mutate(sample = str_replace_all(sample, "-", "\\.")) %>%
  distinct()

tcga_brca_subtypes <- read_tsv(file = "./data/clinical/tcga_publications/tcga_brca/sample_metadata.txt") %>%
  select(sample = `Complete TCGA ID`, subtype = `PAM50 mRNA`) %>%
  mutate(sample = str_replace_all(sample, "-", "\\."))

tcga_brca <- tcga_brca %>%
  left_join(cbio_brca[, c("sample", "ethnicity", "race")], by = "sample") %>%
  mutate(race=tolower(race), ethnicity=tolower(ethnicity)) %>%
  mutate(batch = "tcga-brca", cancer = "brca", tissue = "breast") %>%
  select(sample, batch, tissue, cancer, everything()) %>%
  #left_join(tcga_subtypes1[, c("pan.samplesID", "Subtype_Selected")], by = c("sample" = "pan.samplesID")) %>%
  left_join(tcga_brca_subtypes, by = "sample") %>%
  rename(cancer_subtype = subtype) %>%
  left_join(tcga_survival[, c("bcr_patient_barcode", "OS.time", "OS")], by = c("sample" = "bcr_patient_barcode")) %>%
  rename(OS_time = OS.time) %>%
  mutate(OS_time = as.numeric(OS_time))



# TCGA coread
# cptac data
#primarysiteofdesease contains information about colon/rectum
tcga_coread <- read_csv(file = "./data/clinical/cptac/coread/CPTAC_TCGA_ColorectalCancer_select_clinical_data_release1_090413.csv") %>%
  select(sample=`Hybridization REF`, age=yearstobirth, gender, histological_type = histologicaltype, tumor_stage=neoplasm.diseasestage) %>%
  mutate(sample = toupper(str_replace_all(sample, "-", "\\."))) %>%
  distinct()

cbio_coread <- read_tsv(file = "./data/clinical/cbioportal/coadread_tcga_pan_can_atlas_2018/data_clinical_patient.txt", skip = 4) %>%
  select(sample=PATIENT_ID, cancer=CANCER_TYPE_ACRONYM, cancer_subtype=SUBTYPE, age=AGE, gender=SEX, ajcc_tumor_stage=AJCC_PATHOLOGIC_TUMOR_STAGE, ethnicity=ETHNICITY, race=RACE, seq_center=CENTER) %>%
  mutate(sample = str_replace_all(sample, "-", "\\.")) %>%
  distinct()

tcga_coread_subtypes <- read_tsv(file = "./data/clinical/tcga_publications/tcga_coread/sample_metadata.txt") %>%
  select(sample = `patient`, subtype = expression_subtype) %>%
  mutate(sample = str_replace_all(sample, "-", "\\."))

tcga_coread <- tcga_coread %>%
  left_join(cbio_coread[, c("sample", "ethnicity", "race")], by = "sample") %>%
  mutate(race=tolower(race), ethnicity=tolower(ethnicity)) %>%
  mutate(batch = "tcga-coread", cancer = "coread", tissue = "colorectal") %>%
  select(sample, batch, tissue, cancer, everything()) %>%
  mutate(histological_type = str_replace(histological_type, "colon |rectal ", "")) %>%
  #left_join(tcga_subtypes2[, c("pan.samplesID", "Subtype_Selected")], by = c("sample" = "pan.samplesID")) %>%
  left_join(tcga_coread_subtypes, by = "sample") %>%
  rename(cancer_subtype = subtype) %>%
  left_join(tcga_survival[, c("bcr_patient_barcode", "OS.time", "OS")], by = c("sample" = "bcr_patient_barcode")) %>%
  rename(OS_time = OS.time) %>%
  mutate(OS_time = as.numeric(OS_time))



# TCGA ov
# cptac data
tcga_ov <- read_csv(file = "./data/clinical/cptac/ov/OV_All_clinical_features_TCGAbiotab_CPTAC_S020.csv") %>%
  select(sample=bcr_patient_barcode, age=age_at_initial_pathologic_diagnosis, gender, histological_type, tumor_stage=clinical_stage) %>%
  mutate(sample = str_replace_all(sample, "-", "\\."), gender=tolower(gender), histological_type=tolower(histological_type), tumor_stage=tolower(tumor_stage)) %>%
  distinct()

cbio_ov <- read_tsv(file = "./data/clinical/cbioportal/ov_tcga_pan_can_atlas_2018/data_clinical_patient.txt", skip = 4) %>%
  select(sample=PATIENT_ID, cancer=CANCER_TYPE_ACRONYM, cancer_subtype=SUBTYPE, age=AGE, gender=GENDER, ajcc_tumor_stage=AJCC_PATHOLOGIC_TUMOR_STAGE, ethnicity=ETHNICITY, race=RACE, seq_center=CENTER) %>%
  mutate(sample = str_replace_all(sample, "-", "\\.")) %>%
  distinct()

tcga_ov <- tcga_ov %>%
  left_join(cbio_ov[, c("sample", "ethnicity", "race")], by = "sample") %>%
  mutate(race=tolower(race), ethnicity=tolower(ethnicity)) %>%
  mutate(batch = "tcga-ov", cancer = "hgsc", tissue = "ovary") %>%
  select(sample, batch, tissue, cancer, everything()) %>%
  left_join(tcga_subtypes1[, c("pan.samplesID", "Subtype_mRNA")], by = c("sample" = "pan.samplesID")) %>%
  rename(cancer_subtype = Subtype_mRNA) %>%
  left_join(tcga_survival[, c("bcr_patient_barcode", "OS.time", "OS")], by = c("sample" = "bcr_patient_barcode")) %>%
  rename(OS_time = OS.time) %>%
  mutate(OS_time = as.numeric(OS_time))



# discovery ccrcc
# cptac data

max_days <- function(x,y){
  days <- max(x,y,na.rm = T)
  if(is.infinite(days) | days < 0) days <- NA
  days
}

os_status <- function(x,y){
  os <- max(x,y,na.rm = T)
  if(is.infinite(os)) os <- NA
  os
}

discovery_ccrcc <- read_tsv(file = "./data/clinical/cptac_new/ccrcc/S044_CPTAC_CCRCC_Discovery_Cohort_Clinical_Data_r3_Mar2019.txt") %>%
  select(sample=case_id, age, gender, histological_type=histologic_type, tumor_stage=tumor_stage_pathological, ethnicity, race, height=height_in_cm, weight=weight_in_kg, bmi=BMI, OS_time_12 = days_from_initial_diagnosis_to_last_contact_at_12months_follow_up, OS_time_24 = days_from_initial_diagnosis_to_last_contact_at_24months_follow_up, OS_12 = vital_status_at_12months_follow_up, OS_24 = vital_status_at_24months_follow_up) %>%
  mutate(OS_12 = if_else(OS_12 == "Deceased", 1, 0), OS_24 = if_else(OS_24 == "Deceased", 1, 0)) %>%
  mutate(OS = map2_dbl(.x=OS_12, .y=OS_24, .f = os_status)) %>%
  mutate(OS_time = map2_dbl(.x=OS_time_12, .y=OS_time_24, .f = max_days)) %>%
  select(-OS_12, -OS_24) %>%
  select(-OS_time_12, -OS_time_24) %>%
  mutate(age = str_replace(age, ">=90", "90")) %>%
  mutate(age = as.numeric(age)) %>%
  mutate_if(.predicate = is.character, .funs = tolower) %>%
  mutate(sample = str_replace_all(toupper(sample), "-", "."), ethnicity = str_replace_all(ethnicity, "-", " ")) %>%
  mutate(batch = "discovery-ccrcc", cancer = "ccrcc", tissue = "kidney") %>%
  select(sample, batch, tissue, cancer, everything())

sampleIDs <- read_tsv(file = "./data/clinical/cptac_new/ccrcc/S044_CPTAC_CCRCC_Discovery_Cohort_Specimens_r1_Sept2018.txt") %>%
  select(`ParticipantID`, `Aliquot ID`) %>%
  na.exclude() %>%
  distinct() %>%
  set_names(nm = c("sample", "aliquot")) %>%
  mutate(sample = str_replace_all(toupper(sample), "-", ".")) %>%
  filter(str_detect(sample, "^NCI", negate = T))

subtypes <- read_tsv(file = "./data/clinical/cptac_new/ccrcc/xcell_subtypes.txt", col_names = F) %>%
  t() %>%
  as_tibble() %>%
  slice(-1) %>%
  set_names(nm = c("aliquot", "group")) %>%
  filter(str_detect(group, "NAT", negate = T)) %>%
  left_join(sampleIDs, by = "aliquot") %>%
  select(sample, group)

discovery_ccrcc <- discovery_ccrcc %>%
  left_join(subtypes, by = "sample") %>%
  rename(cancer_subtype = group)



# discovery ucec
# cptac data

ucec_gdc <- read_tsv(file = "./data/clinical/cptac_new/ucec/clinical.tsv") %>%
  select(sample = submitter_id, days_to_last_known_disease_status, days_to_last_follow_up, OS = vital_status) %>%
  mutate(OS = if_else(OS == "Dead", 1, 0)) %>%
  mutate(sample = str_replace_all(sample, "-", "\\.")) %>%
  filter(!c(days_to_last_known_disease_status == "--" | days_to_last_follow_up == "--")) %>%
  mutate(days_to_last_known_disease_status = as.numeric(days_to_last_known_disease_status)) %>%
  filter(days_to_last_known_disease_status < 19000) %>%
  mutate(days_to_last_follow_up = as.numeric(days_to_last_follow_up)) %>%
  mutate(OS_time = map2_dbl(.x=days_to_last_known_disease_status, .y=days_to_last_follow_up, .f = max_days)) %>%
  select(-days_to_last_known_disease_status, -days_to_last_follow_up)
  

discovery_ucec <- read_tsv(file = "./data/clinical/cptac_new/ucec/S043_CPTAC_UCEC_Discovery_Cohort_Clinical_Data_r1_Sept2018.txt") %>%
  select(sample=Case_ID, age=Age, gender=Gender, histological_type=Histologic_Type, ethnicity=Ethnicity_Self_Identify, race=Race, height=Height_cm, weight=Weight_kg, bmi=BMI) %>%
  mutate(age = str_replace(age, ">=90", "90")) %>%
  mutate(age = as.numeric(age)) %>%
  mutate_if(.predicate = is.character, .funs = tolower) %>%
  mutate(sample = str_replace_all(toupper(sample), "-", ".")) %>%
  mutate(batch = "discovery-ucec", cancer = "ucec", tissue = "uterus") %>%
  select(sample, batch, tissue, cancer, everything()) %>%
  left_join(ucec_gdc, by = "sample")

subtypes1 <- read_tsv(file = "./data/clinical/cptac_new/ucec/1-s2.0-S0092867420301070-mmc1.txt") %>%
  select(idx, sample = Proteomics_Participant_ID, tumor_stage = `tumor_Stage-Pathological`, cancer_subtype1 = Genomics_subtype) %>%
  mutate(sample = str_replace(sample, "-", "."), tumor_stage = tolower(tumor_stage)) %>%
  na.exclude()

subtypes2 <- read_tsv(file = "./data/clinical/cptac_new/ucec/1-s2.0-S0092867420301070-mmc7.txt") %>%
  select(idx, cancer_subtype2 = Immune_cluster)

subtypes1 <- subtypes1 %>%
  left_join(subtypes2, by = "idx") %>%
  select(-idx)

discovery_ucec <- discovery_ucec %>%
  left_join(subtypes1, by = "sample") %>%
  rename(cancer_subtype = cancer_subtype1) %>%
  select(-cancer_subtype2)



# discovery luad
# cptac data
discovery_luad <- read_tsv(file = "./data/clinical/cptac_new/luad/S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.txt") %>%
  select(sample=case_id, age, gender, histological_type=histologic_type, tumor_stage=tumor_stage_pathological, ethnicity, race, height=height_in_cm, weight=weight_in_kg, bmi=BMI, OS_time_12 = days_from_initial_diagnosis_to_last_contact_at_12months_follow_up, OS_time_24 = days_from_initial_diagnosis_to_last_contact_at_24months_follow_up, OS_12 = vital_status_at_12months_follow_up, OS_24 = vital_status_at_24months_follow_up) %>%
  filter(!is.na(sample)) %>%
  mutate(OS_12 = if_else(OS_12 == "Deceased", 1, 0), OS_24 = if_else(OS_24 == "Deceased", 1, 0)) %>%
  mutate(OS = map2_dbl(.x=OS_12, .y=OS_24, .f = os_status)) %>%
  mutate(OS_time = map2_dbl(.x=OS_time_12, .y=OS_time_24, .f = max_days)) %>%
  select(-OS_time_12, -OS_time_24) %>%
  select(-OS_12, -OS_24) %>%
  mutate_if(.predicate = is.character, .funs = tolower) %>%
  mutate(sample = str_replace_all(toupper(sample), "-", "."), ethnicity = str_replace_all(ethnicity, "-", " ")) %>%
  mutate(batch = "discovery-luad", cancer = "luad", tissue = "lung") %>%
  select(sample, batch, tissue, cancer, everything())



# cbttc
# cptac data
# let's assume the cancer subtypes are the different tumour types (medulloblastoma, ependymoma, etc)

folders <- c("atrt_cbttc", "craniopharyngioma_cbttc", "ependymoma_cbttc", "ganglioglioma_cbttc", "medulloblastoma_cbttc", "ped_high_grade_glioma_cbttc", "ped_low_grade_glioma_cbttc")

read_cbttc <- function(folder){
  
  file = str_c("./data/clinical/cbioportal/cbttc", folder, "data_clinical_patient.txt", sep = "/")
  
  dat <- read_tsv(file) %>%
    slice(-1:-4) %>%
    select(sample = `#Patient Identifier`, OS_time = OS_MONTHS, OS = OS_STATUS) %>%
    mutate(OS_time = round(as.numeric(OS_time))*30) %>%
    mutate(OS = if_else(OS == "DECEASED", 1, 0))
  
  dat
}

cbttc_survival <- map_dfr(.x = folders, .f = read_cbttc) %>%
  distinct()

cbttc <- read_tsv(file = "./data/clinical/cptac_new/cbttc/S047_Pediatric_Brain_Cancer_Clinical_Data_r1.txt") %>%
  select(sample=Kids.First.ID, age=age_at_diagnosis_days, gender, ethnicity = Ethnicity, race, cancer=diagnosis) %>%
  mutate_if(.predicate = is.character, .funs = tolower) %>%
  mutate(sample = toupper(sample), age = round(age/365)) %>%
  mutate(batch = "cbttc", tissue = "brain") %>%
  select(sample, batch, tissue, cancer, everything()) %>%
  mutate(cancer_subtype = cancer) %>%
  left_join(cbttc_survival, by = "sample")



# eogc
# paper clinical data
eogc_survival <- read_tsv(file = "./data/clinical/cptac_new/eogc/clinical_information_EOGC.txt") %>%
  select(tumor=Tumor, normal=Normal, OS_time = OS_month, OS = Death) %>%
  mutate(tumor = str_c("T",str_replace(tumor, "T", "")), normal = str_c("N",str_replace(normal, "N", ""))) %>%
  unite(col = "sample", c("tumor", "normal"), sep = "") %>%
  mutate(OS_time = round(as.numeric(OS_time))*30)

eogc <- read_tsv(file = "./data/clinical/cptac_new/eogc/clinical_data_paper.txt") %>%
  select(tumor=Tumor, normal=Normal, gender = Gender, age = `Age at the time of surgery`, histological_type = Histology, cancer_subtype = `Integrated subtype`) %>%
  mutate(tumor = str_c("T",str_replace(tumor, "T", "")), normal = str_c("N",str_replace(normal, "N", ""))) %>%
  unite(col = "sample", c("tumor", "normal"), sep = "") %>%
  mutate(gender = if_else(gender == "M", "male", "female")) %>%
  mutate(histological_type = str_replace(histological_type, "N/A", "other")) %>%
  mutate(batch = "eogc", cancer = "gc", tissue = "stomach") %>%
  select(sample, batch, tissue, cancer, everything()) %>%
  left_join(eogc_survival, by = "sample")



# colon
# linkedomics data
# only colon tumour (not rectal)
colon <- read_tsv(file = "./data/clinical/cptac_new/colon/Human__CPTAC_COAD__MS__Clinical__Clinical__03_01_2017__CPTAC__Clinical__BCM.tsi") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  as_tibble() %>%
  mutate_if(is.factor, as.character) %>%
  set_names(nm = as.character(slice(., 1))) %>%
  slice(-1) %>%
  rename(histological_type=Mucinous) %>%
  mutate(histological_type = if_else(histological_type == "Mucinous", "mucinous adenocarcinoma", "adenocarcinoma")) %>%
  select(sample = attrib_name, age = Age, gender = Gender, histological_type, tumor_stage = Stage, cancer_subtype = Integrated.Phenotype) %>%
  mutate(gender = tolower(gender), tumor_stage = tolower(tumor_stage)) %>%
  mutate(age = as.numeric(age)) %>%
  mutate(age = replace_na(age, round(mean(age, na.rm=T)))) %>%
  mutate(age=round(age/12)) %>%
  mutate(batch = "colon-oppt", cancer = "coread", tissue = "colorectal") %>%
  select(sample, batch, tissue, cancer, everything())

subtypes <- read_tsv(file = "./data/clinical/cptac_new/colon/1-s2.0-S0092867419302922-mmc1.txt") %>%
  slice(-1) %>%
  select(sample = SampleID, UMS) %>%
  mutate(UMS = replace_na(UMS, "Others"))
  # others should be NA

colon <- colon %>%
  inner_join(subtypes, by = "sample") %>%
  select(-cancer_subtype) %>%
  rename(cancer_subtype = UMS) %>%
  mutate(OS_time = NA, OS = NA)


# ccle
# broad institute annotations
ccle <- read_tsv(file = "./data/clinical/broad/ccle/Cell_lines_annotations_20181226.txt") %>%
  select(CCLE_ID, tissue=type_refined, cancer=tcga_code, age=Age, gender=Gender, race=Race, ethnicity=inferred_ethnicity, histological_type=Hist_Subtype1) %>%
  mutate(cancer = tolower(cancer)) %>%
  mutate(cancer = str_replace(cancer, "coad/read", "coread")) %>%
  mutate(histological_type = str_replace(histological_type, "NS", "other")) %>%
  mutate(histological_type = str_replace_na(histological_type, "other")) %>%
  filter(cancer %in% c("coread", "brca") & tissue %in% c("colorectal", "breast")) %>%
  separate(col = "CCLE_ID", into = c("sample", "tissue2"), sep = "_", extra = "merge") %>%
  select(-tissue2) %>%
  mutate(sample = str_replace(sample, "COLO320", "COLO320HSR")) %>%
  mutate(sample = str_replace(sample, "MDAMB175VII", "MDAMB175")) %>%
  mutate(sample = str_replace(sample, "MDAMB134VI", "MDAMB134")) %>%
  #filter(!is.na(gender)) %>%
  mutate(age = replace_na(age, round(mean(age, na.rm=T))), batch = "ccle") %>%
  select(sample, batch, tissue, cancer, everything())


subtypes <- read_tsv(file = "./data/clinical/broad/ccle/roumeliotis_subtypes.txt") %>%
  mutate(cell_line = str_replace_all(cell_line, "-", "")) %>%
  mutate(cell_line = toupper(cell_line)) %>%
  rename(sample = cell_line, cancer_subtype = cluster)

ccle <- ccle %>%
  left_join(subtypes, by = "sample") %>%
  mutate(OS_time = NA, OS = NA)



# HBV-HCC
# cptac data

hcc_survival <- read_tsv(file = "./data/clinical/cptac_new/hbv-hcc/S049_HCC_Clinical_Information_and_TMT11_Sample_Mapping_Gao2019.r1_2.txt") %>%
  select(case_ID = `Case ID`, OS_time = `Survival time (day)`, OS = `Survival Test`)

hcc_paper <- read_tsv(file = "./data/clinical/cptac_new/hbv-hcc/1-s2.0-S0092867419310037-mmc1.txt", col_names = F) %>%
  select(1, 7, 9, 25, 27) %>%
  set_names(nm = slice(., 1)) %>%
  slice(-1) %>%
  set_names(c("sample", "OS_time", "OS", "tumor_stage", "cancer_subtype")) %>%
  mutate(OS_time = round(as.numeric(OS_time))*30) %>%
  mutate(OS = as.numeric(OS)) %>%
  mutate(tumor_stage = str_c("stage", tolower(tumor_stage), sep = " "))

hcc <- read_tsv(file = "./data/clinical/cptac_new/hbv-hcc/S049_HCC_Clinical_Information_and_TMT11_Sample_Mapping_Gao2019.r1.txt") %>%
  select(sample=`Tumor ID in Proteome experiment`, age=Age, gender=Gender) %>%
  mutate(sample = str_c("T", sample)) %>%
  mutate(batch = "hbv-hcc", cancer = "hcc", tissue = "liver") %>%
  select(sample, batch, tissue, cancer, everything()) %>%
  left_join(hcc_paper, by = "sample")



all_clinical_samples <- bind_rows(
  tcga_brca,
  tcga_coread,
  tcga_ov,
  discovery_ccrcc,
  discovery_luad,
  discovery_ucec,
  cbttc,
  eogc,
  colon,
  ccle,
  hcc) %>%
  select(1:4)

write.table(all_clinical_samples, "./output/files/clinical_samples.txt", sep="\t", quote=F, row.names=F)



all_clinical <- bind_rows(
  tcga_brca,
  tcga_coread,
  tcga_ov,
  discovery_ccrcc,
  discovery_luad,
  discovery_ucec,
  cbttc,
  eogc,
  colon,
  ccle,
  hcc) %>%
  filter(sample %in% proteomics$sample)

write.table(all_clinical, "./output/files/clinical_data.txt", sep="\t", quote=F, row.names=F)
