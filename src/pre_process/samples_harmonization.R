library(tidyverse)


# load multi-omics cancer samples

# mutation
mut <- read_tsv(file = "./output/files/mutations_samples.txt")

# copy-number variation
cnv <- read_tsv(file = "./output/files/cnv_samples.txt")

# trancriptomics
rna <- read_tsv(file = "./output/files/transcriptomics_samples.txt")

# proteomics
proteomics <- read_tsv(file = "./output/files/proteomics_samples.txt")

# phopshoproteomics
phospho <- read_tsv(file = "./output/files/phosphoproteomics_samples.txt")

# clinical
clinical <- read_tsv(file = "./output/files/clinical_samples.txt")


# harmonize batch, tissue and cancer identifiers across datasets
clinical <- clinical %>%
  #mutate(batch_original = batch) %>%
  mutate(batch = str_replace(batch, "eogc", "eogc-proteogenomics")) %>%
  mutate(batch = str_replace(batch, "hbv-hcc", "hcc-proteogenomics")) %>%
  mutate(batch = str_replace(batch, "colon-oppt", "colon-opportunities")) %>%
  mutate(data = "clinical")

cancer_tissue <- clinical %>%
  group_by(batch, cancer, tissue) %>%
  tally() %>%
  ungroup() %>%
  select(-n)

mut <- mut %>%
  #mutate(batch_original = batch) %>%
  mutate(batch = tolower(batch), cancer = tolower(cancer)) %>%
  mutate(cancer = str_replace(cancer, "eogc", "gc")) %>%
  mutate(cancer = str_replace(cancer, "hbv-hcc", "hcc")) %>%
  mutate(batch = str_replace(batch, "gc", "eogc-proteogenomics")) %>%
  mutate(batch = str_replace(batch, "hbv-hcc", "hcc-proteogenomics")) %>%
  mutate(batch = str_replace(batch, "tcga-hgsc", "tcga-ov")) %>%
  inner_join(cancer_tissue, by = c("batch", "cancer")) %>%
  mutate(data = "mutation")

cnv <- cnv %>%
  #mutate(batch_original = batch) %>%
  mutate(batch = tolower(batch), cancer = tolower(cancer)) %>%
  mutate(batch = str_replace(batch, "hcc", "hcc-proteogenomics")) %>%
  mutate(batch = str_replace(batch, "tcga-hgsc", "tcga-ov")) %>%
  inner_join(cancer_tissue, by = c("batch", "cancer")) %>%
  mutate(data = "cnv")

rna <- rna %>%
  #mutate(batch_original = batch) %>%
  mutate(batch = tolower(batch), cancer = tolower(cancer)) %>%
  mutate(cancer = str_replace(cancer, "eogc", "gc")) %>%
  mutate(batch = str_replace(batch, "gc", "eogc-proteogenomics")) %>%
  mutate(batch = str_replace(batch, "hcc", "hcc-proteogenomics")) %>%
  mutate(batch = str_replace(batch, "tcga-hgsc", "tcga-ov")) %>%
  inner_join(cancer_tissue, by = c("batch", "cancer")) %>%
  mutate(data = "mRNA")

proteomics <- proteomics %>%
  #mutate(batch_original = batch) %>%
  mutate(batch = str_replace(batch, "cell-lines-law|cell-lines-lpk|cell-lines-rmlt", "ccle")) %>%
  mutate(batch = str_replace(batch, "cptac-", "")) %>%
  mutate(batch = tolower(batch), cancer = tolower(cancer)) %>%
  mutate(cancer = str_replace(cancer, "eogc", "gc")) %>%
  mutate(batch = str_replace(batch, "gc", "eogc-proteogenomics")) %>%
  mutate(batch = str_replace(batch, "hcc", "hcc-proteogenomics")) %>%
  mutate(batch = str_replace(batch, "tcga-hgsc", "tcga-ov")) %>%
  inner_join(cancer_tissue, by = c("batch", "cancer")) %>%
  mutate(data = "protein")

phospho <- phospho %>%
  #mutate(batch_original = batch) %>%
  mutate(batch = str_replace(batch, "cell-lines-rmlt", "ccle")) %>%
  mutate(batch = str_replace(batch, "cptac-", "")) %>%
  mutate(batch = tolower(batch), cancer = tolower(cancer)) %>%
  mutate(cancer = str_replace(cancer, "eogc", "gc")) %>%
  mutate(cancer = str_replace(cancer, "colon-cancer", "coread")) %>%
  mutate(batch = str_replace(batch, "gc", "eogc-proteogenomics")) %>%
  mutate(batch = str_replace(batch, "hcc", "hcc-proteogenomics")) %>%
  mutate(batch = str_replace(batch, "tcga-hgsc", "tcga-ov")) %>%
  inner_join(cancer_tissue, by = c("batch", "cancer")) %>%
  mutate(data = "phosphorylation")


# harmonized annotation for all samples in all omics
all_samples_annotation <- bind_rows(mut, cnv, rna, proteomics, phospho, clinical)
write.table(all_samples_annotation, "./output/files/all_samples_annotation.txt", sep="\t", quote=F, row.names=F)
