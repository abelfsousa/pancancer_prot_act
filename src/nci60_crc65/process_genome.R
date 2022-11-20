library(tidyverse)

source("./src/utils/read_hgvsp_nci60_crc65.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# -- DNA data


# load cell lines
nci60_crc65 <- read_tsv("./output/files/nci60_crc65_cell_lines.txt")


# NCI60 cell lines
nci60_dna <- read_tsv("./data/nci60_crc65/genome/nci60/nci60_dna.txt") %>%
  select(-c(3:8,10:11,14:16)) %>%
  select(identifier = `Identifier (c)`, gene=`Gene name (d)`, HGVSp = `AA Impact (h)`, mutation_type = `SNP Type (l)`, sift = `SIFT Score (k)`, everything()) %>%
  pivot_longer(-c(identifier, gene, HGVSp, mutation_type, sift), names_to = "sample", values_to = "p_alleles") %>%
  filter(p_alleles != "-") %>%
  mutate(p_alleles = as.numeric(p_alleles)) %>%
  filter(p_alleles > 50) %>%
  filter(mutation_type %in% c("Missense", "Nonsense", "Initiation_loss", "Read_through", "Frameshift", "Nonframeshift", "Silent")) %>%
  mutate(res = map2(.x=mutation_type, .y=HGVSp, .f = read_hgvsp)) %>%
  unnest() %>%
  mutate(sample = str_split_fixed(sample, ":", 2)[,2]) %>%
  mutate(sample = toupper(str_replace_all(sample, "-| ", ""))) %>%
  mutate(batch = "NCI60") %>%
  select(batch, sample, everything())


# CRC65 cell lines
depmap_info <- data.table::fread(file = "./data/nci60_crc65/genome/crc65/sample_info.csv") %>%
  as_tibble()

crc65_dna <- data.table::fread("./data/nci60_crc65/genome/crc65/CCLE_mutations.csv.gz") %>%
  as_tibble() %>%
  inner_join(depmap_info[,1:2], by = c("DepMap_ID")) %>%
  select(stripped_cell_line_name, everything()) %>%
  select(-DepMap_ID, sample = stripped_cell_line_name) %>%
  filter(sample %in% nci60_crc65[nci60_crc65$batch == "CRC65", "cell_line", drop=T]) %>%
  select(sample, mutation_type=Variant_Classification, identifier=Genome_Change, gene=Hugo_Symbol, HGVSp=Protein_Change) %>%
  filter(mutation_type %in%
           c("Missense_Mutation",
             "Nonsense_Mutation",
             "Start_Codon_Del",
             "Start_Codon_Ins",
             "Start_Codon_SNP",
             "Nonstop_Mutation",
             "Stop_Codon_Del",
             "Stop_Codon_Ins",
             "Frame_Shift_Del",
             "Frame_Shift_Ins",
             "In_Frame_Del",
             "In_Frame_Ins",
             "Silent")) %>%
  mutate(HGVSp = str_replace(HGVSp, "^p\\.", "")) %>%
  mutate(res = map2(.x=mutation_type, .y=HGVSp, .f = read_hgvsp)) %>%
  unnest() %>%
  mutate(mutation_type = str_replace(mutation_type, "Missense_Mutation", "Missense")) %>%
  mutate(mutation_type = str_replace(mutation_type, "Nonsense_Mutation", "Nonsense")) %>%
  mutate(mutation_type = str_replace(mutation_type, "Frame_Shift_Ins|Frame_Shift_Del", "Frameshift")) %>%
  mutate(mutation_type = str_replace(mutation_type, "In_Frame_Del|In_Frame_Ins", "Nonframeshift")) %>%
  mutate(mutation_type = str_replace(mutation_type, "Nonstop_Mutation|Stop_Codon_Ins|Stop_Codon_Del", "Read_through")) %>%
  mutate(mutation_type = str_replace(mutation_type, "Start_Codon_SNP|Start_Codon_Del|Start_Codon_Ins", "Initiation_loss")) %>%
  mutate(sift = NA, p_alleles = NA, batch = "CRC65") %>%
  select(batch, sample, everything())


mutations <- bind_rows(nci60_dna, crc65_dna)
write_tsv(mutations, "./output/files/nci60_crc65_mutations.txt.gz")
