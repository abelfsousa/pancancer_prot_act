# regress-out mRNA from phosphosites

# load R packages
library(tidyverse)

source("./src/utils/regress_outCovs.R")
source("./src/utils/getSamples.R")

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy



# load all samples gathered in this study
all_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")

# get samples with phosphorylation and protein data
samplesPR <- getSamples(all_samples, c("phosphorylation", "mRNA", "clinical")) %>%
  select(sample, batch)


# load mRNA data (log2FC)
rna <- read_tsv(file = "./output/files/transcriptomics_log2fc.txt.gz") %>%
  pivot_longer(names_to = "sample", values_to = "rna_log2fc", -gene) %>%
  select(sample, protein=gene, rna_log2fc)

# load quantile normalized phosphorylation data (log2FC)
# remove some phosphosites with NAs in all samples
phospho <- read_tsv(file = "./output/files/phosphoproteomicsQ.txt.gz") %>%
  pivot_longer(names_to = "sample", values_to = "phos_log2fc", -c(gene, psite, psites)) %>%
  select(sample, protein=gene, psite, psites, phos_log2fc) %>%
  group_by(protein, psite, psites) %>%
  filter( !(sum(is.na(phos_log2fc)) == n()) ) %>%
  ungroup()

# load metadata
metadata <- read_tsv(file = "./output/files/clinical_data.txt") %>%
  select(sample, age, gender) %>%
  mutate(gender = replace_na(gender, "other"))



# add phosphorylation and mRNA data to respective samples
pho_rna <- samplesPR %>%
  inner_join(rna, by = "sample") %>%
  inner_join(phospho, by = c("sample", "protein")) %>%
  select(batch, sample, prot=protein, psite, psites, prot_fc=rna_log2fc, psite_fc=phos_log2fc)


# clear data
rm(phospho, rna)



# regress-out mRNA from phosphosites
reg <- pho_rna %>%
  select(prot, psite, psites, sample, psite_fc, prot_fc) %>%
  group_by(psite, psites, prot) %>%
  nest() %>%
  ungroup() %>%
  mutate(residuals = map(.x = data, .f = regress_outCovs, covs = metadata, ztransf = TRUE)) %>%
  #mutate(residuals = map(.x = data, .f = regress_outCovs, ztransf = TRUE)) %>%
  unnest()



# set up and export phosphorylation data
phospho <- reg %>%
  select(-psite_fc, -prot_fc) %>%
  group_by(prot, psite, psites) %>%
  filter( !(sum(is.na(residual)) == n()) ) %>%
  ungroup() %>%
  rename(gene = prot) %>%
  pivot_wider(names_from = "sample", values_from = "residual") %>%
  select(gene, psite, psites, everything())

gz1 <- gzfile("./output/files/phosphoproteomicsQ_RNAreg.txt.gz", "w")
write.table(phospho, gz1, sep="\t", quote=F, row.names=F)
close(gz1)
