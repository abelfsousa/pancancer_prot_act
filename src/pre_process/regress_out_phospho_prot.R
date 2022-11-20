# regress-out protein from phosphosites

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
samplesPP <- getSamples(all_samples, c("phosphorylation", "protein", "clinical")) %>%
  select(sample, batch)


# load protein data (log2FC)
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(names_to = "sample", values_to = "prot_log2fc", -gene) %>%
  select(sample, protein=gene, prot_log2fc)


# load phosphorylation data (log2FC)
# remove phosphosites with NAs in all samples
args <- commandArgs(trailingOnly = TRUE)
#args <- c("./output/files/phosphoproteomicsQ.txt.gz", "./output/files/", "phosphoproteomicsQ_Protreg", "allSamp", "withZtransf")

phospho <- read_tsv(file = args[1]) %>%
  pivot_longer(names_to = "sample", values_to = "phos_log2fc", -c(gene, psite, psites)) %>%
  select(sample, protein=gene, psite, psites, phos_log2fc) %>%
  group_by(protein, psite, psites) %>%
  filter( !(sum(is.na(phos_log2fc)) == n()) ) %>%
  ungroup()


# load metadata
metadata <- read_tsv(file = "./output/files/clinical_data.txt") %>%
  select(sample, age, gender) %>%
  mutate(gender = replace_na(gender, "other"))


# add phosphorylation and protein data to respective samples
pho_prot <- samplesPP %>%
  inner_join(protein, by = "sample") %>%
  inner_join(phospho, by = c("sample", "protein")) %>%
  select(batch, sample, prot=protein, psite, psites, prot_fc=prot_log2fc, psite_fc=phos_log2fc)


# clear data
rm(phospho, protein)


# regress-out protein from phosphosites
if(args[5] == "withZtransf"){
  zeta = TRUE
}else if(args[5] == "withoutZtransf"){
  zeta = FALSE
}else{
  stop("Flag not recognized!")
}

if(args[4] == "allSamp"){
  reg <- pho_prot %>%
    select(prot, psite, psites, sample, psite_fc, prot_fc) %>%
    group_by(psite, psites, prot) %>%
    nest() %>%
    ungroup() %>%
    mutate(residuals = map(.x = data, .f = regress_outCovs, covs = metadata, ztransf = zeta)) %>%
    #mutate(residuals = map(.x = data, .f = regress_outCovs, ztransf = TRUE)) %>%
    unnest()
} else if(args[4] == "byBatch"){
  reg <- pho_prot %>%
    select(batch, prot, psite, psites, sample, psite_fc, prot_fc) %>%
    group_by(batch, psite, psites, prot) %>%
    nest() %>%
    ungroup() %>%
    mutate(residuals = map(.x = data, .f = regress_outCovs, covs = metadata, ztransf = zeta)) %>%
    #mutate(residuals = map(.x = data, .f = regress_outCovs, ztransf = TRUE)) %>%
    unnest() %>%
    select(-batch)
} else{
  stop("Flag not recognized!")
}


# set up and export phosphorylation data
phospho <- reg %>%
  select(-psite_fc, -prot_fc) %>%
  group_by(prot, psite, psites) %>%
  filter( !(sum(is.na(residual)) == n()) ) %>%
  ungroup() %>%
  rename(gene = prot) %>%
  pivot_wider(names_from = "sample", values_from = "residual") %>%
  select(gene, psite, psites, everything())

file_name <- paste0(args[2], "/", args[3], "_", args[4], "_", args[5], ".txt.gz")
gz1 <- gzfile(file_name, "w")
write.table(phospho, gz1, sep="\t", quote=F, row.names=F)
close(gz1)
