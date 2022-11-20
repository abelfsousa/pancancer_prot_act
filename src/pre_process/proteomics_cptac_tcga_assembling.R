library(tidyverse)


# load cptac TCGA BRCA cancer study
cptac_tcga_breast <- read_tsv("/Volumes/G-DRIVE/data/cptac/tcga/TCGA Breast Cancer/TCGA_Breast_BI_Proteome_CDAP_Protein_Report.r3/Protein_data/CDAP/TCGA_Breast_BI_Proteome.itraq.tsv") %>%
  slice(-c(1:3)) %>%
  rename(gene = Gene) %>%
  select(-c("NCBIGeneID", "Authority", "Description", "Organism", "Chromosome", "Locus")) %>%
  select(c("gene", contains("Unshared"))) %>%
  select(c(1, 5:112)) %>%
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  gather(key = "sample", value = "tmt", -gene) %>%
  mutate(sample = str_split_fixed(sample, "\\.", n=2)[,1]) %>%
  group_by(sample, gene) %>%
  summarise(tmt = mean(tmt, na.rm = T)) %>%
  ungroup() %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(sample = str_sub(sample, 1, 10)) %>%
  mutate(sample = paste0("TCGA.", sample)) %>%
  spread(key = "sample", value = "tmt")



# load cptac TCGA COREAD cancer study
cptac_tcga_coread <- read_tsv("/Volumes/G-DRIVE/data/cptac/tcga/TCGA Colorectal Cancer/TCGA_Colon_VU_Proteome_CDAP_Protein_Report.r2/Protein_data/CDAP/TCGA_Colon_VU_Proteome_CDAP.r2.precursor_area.tsv") %>%
  rename(gene = Gene) %>%
  select(-c("Description", "Organism", "Chromosome", "Locus")) %>%
  select(c("gene", contains("Unshared"))) %>%
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  gather(key = "sample", value = "area", -gene) %>%
  group_by(sample) %>%
  mutate(area = area/sum(area)*1e6) %>%
  group_by(gene) %>%
  mutate(log2fc = log2((area/mean(area))+1)) %>%
  select(-area) %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(sample = str_sub(sample, 1, 10)) %>%
  mutate(sample = paste0("TCGA.", sample)) %>%
  group_by(sample, gene) %>%
  summarise(log2fc = mean(log2fc, na.rm = T)) %>%
  ungroup() %>%
  spread(key = "sample", value = "log2fc")



# load cptac TCGA OV cancer study
cptac_tcga_ov_jhu <- read_tsv("/Volumes/G-DRIVE/data/cptac/tcga/TCGA Ovarian Cancer/TCGA_Ovarian_JHU_Proteome_CDAP_Protein_Report.r3/Protein_data/CDAP/TCGA_Ovarian_JHU_Proteome.itraq.tsv") %>%
  slice(-c(1:3)) %>%
  rename(gene = Gene) %>%
  select(-c("NCBIGeneID", "Authority", "Description", "Organism", "Chromosome", "Locus")) %>%
  select(c("gene", contains("Unshared"))) %>%
  select(-contains("CONTROL")) %>%
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  gather(key = "sample", value = "tmt", -gene) %>%
  #mutate(sample = str_split_fixed(sample, "\\.", n=2)[,1]) %>%
  #group_by(sample, gene) %>%
  #summarise(tmt = mean(tmt, na.rm = T)) %>%
  #ungroup() %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(sample = str_sub(sample, 1, 10)) %>%
  mutate(sample = paste0("TCGA.", sample)) %>%
  spread(key = "sample", value = "tmt")


cptac_tcga_ov_pnnl <- read_tsv("/Volumes/G-DRIVE/data/cptac/tcga/TCGA Ovarian Cancer/TCGA_Ovarian_PNNL_Proteome_CDAP_Protein_Report.r3/Protein_data/CDAP/TCGA_Ovarian_PNNL_Proteome.itraq.tsv") %>%
  slice(-c(1:3)) %>%
  rename(gene = Gene) %>%
  select(-c("NCBIGeneID", "Authority", "Description", "Organism", "Chromosome", "Locus")) %>%
  select(c("gene", contains("Unshared"))) %>%
  rename_all(.funs = ~ str_split_fixed(.x, " ", n=2)[,1]) %>%
  gather(key = "sample", value = "tmt", -gene) %>%
  #mutate(sample = str_split_fixed(sample, "\\.", n=2)[,1]) %>%
  #group_by(sample, gene) %>%
  #summarise(tmt = mean(tmt, na.rm = T)) %>%
  #ungroup() %>%
  mutate(sample = str_replace_all(sample, "-", ".")) %>%
  mutate(sample = str_sub(sample, 1, 10)) %>%
  mutate(sample = paste0("TCGA.", sample)) %>%
  spread(key = "sample", value = "tmt")



cptac_tcga <- cptac_tcga_breast %>%
  full_join(cptac_tcga_coread, by = "gene") %>%
  full_join(cptac_tcga_ov_jhu, by = "gene") %>%
  full_join(cptac_tcga_ov_pnnl, by = "gene")

