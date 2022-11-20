# regress-out batch from phosphosites

# load R packages
library(tidyverse)

source("./src/utils/regress_outCovs2.R")
source("./src/utils/getSamples.R")

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load all samples gathered in this study
all_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt")

# get samples with phosphorylation data
samplesPP <- getSamples(all_samples, c("phosphorylation", "clinical")) %>%
  mutate(batch = map2_chr(batch, tissue, ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(sample, batch)


# load metadata
metadata <- read_tsv(file = "./output/files/clinical_data.txt") %>%
  select(sample, age, gender) %>%
  mutate(gender = replace_na(gender, "other"))


# load phosphorylation data (log2FC)
# remove phosphosites with NAs in all samples
# add batch covariate
#args <- commandArgs(trailingOnly = TRUE)
args <- c("./output/files/phosphoproteomics.txt.gz", "./output/files", "phosphoproteomics_BatchReg", "allSamp", "withZtransf")

phospho <- read_tsv(file = args[1]) %>%
  pivot_longer(names_to = "sample", values_to = "phos_log2fc", -c(gene, psite, psites)) %>%
  select(sample, gene, psite, psites, phos_log2fc) %>%
  group_by(gene, psite, psites) %>%
  filter( !(sum(is.na(phos_log2fc)) == n()) ) %>%
  ungroup() %>%
  inner_join(samplesPP, by = "sample")


# regress-out batch from phosphosites
if(args[5] == "withZtransf"){
  zeta = TRUE
}else if(args[5] == "withoutZtransf"){
  zeta = FALSE
}else{
  stop("Flag not recognized!")
}

if(args[4] == "allSamp"){
  reg <- phospho %>%
    select(psite, psites, gene, sample, phos_log2fc, batch) %>%
    group_by(psite, psites, gene) %>%
    nest() %>%
    ungroup() %>%
    mutate(residuals = map(.x = data, .f = regress_outCovs2, covs = metadata, ztransf = zeta)) %>%
    #mutate(residuals = map(.x = data, .f = regress_outCovs2, ztransf = TRUE)) %>%
    unnest()
} else {
  stop("Flag not recognized!")
}


# set up and export phosphorylation data
phospho <- reg %>%
  select(-phos_log2fc, -batch) %>%
  group_by(gene, psite, psites) %>%
  filter( !(sum(is.na(residual)) == n()) ) %>%
  ungroup() %>%
  pivot_wider(names_from = "sample", values_from = "residual") %>%
  select(gene, psite, psites, everything())

file_name <- paste0(args[2], "/", args[3], "_", args[4], "_", args[5], ".txt.gz")
gz1 <- gzfile(file_name, "w")
write.table(phospho, gz1, sep="\t", quote=F, row.names=F)
close(gz1)
