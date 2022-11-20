# load R packages
library(tidyverse)
source("./src/utils/psite_prot_seq_match.R")


# load kinase-substrate list
ks <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz")


# in vitro KS
in_vt_substrates <- ks %>%
  filter(source_type == "in_vitro") %>%
  select(substrate) %>%
  distinct()
write.table(in_vt_substrates, "./output/files/in_vt_substrates.txt", sep="\t", quote=F, row.names=F, col.names = F)

in_vt_map <- read_tsv(file = "./data/kinase_substrate/NKIN/in_vt_substrates_uni_mapping.txt") %>%
  select(substrate = `yourlist:M20191122216DA2B77BFBD2E6699CA9B6D1C41EB25641D2L`, uniprotID = `Entry name`) %>%
  group_by(substrate) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)

in_vt_substrates_toKIN <- ks %>%
  filter(source_type == "in_vitro") %>%
  select(substrate, position, residue) %>%
  distinct() %>%
  inner_join(in_vt_map, by = "substrate") %>%
  select(uniprotID, position, residue) %>%
  filter(!(uniprotID == "KRI1_HUMAN" | uniprotID == "TIGD5_HUMAN"))
write.table(in_vt_substrates_toKIN, "./output/files/in_vt_substrates_toKIN.txt", sep="\t", quote=F, row.names=F, col.names = F)




# in vivo KS
in_vv_substrates <- ks %>%
  filter(source_type == "in_vivo") %>%
  select(substrate) %>%
  distinct()
write.table(in_vv_substrates, "./output/files/in_vv_substrates.txt", sep="\t", quote=F, row.names=F, col.names = F)

in_vv_map <- read_tsv(file = "./data/kinase_substrate/NKIN/in_vv_substrates_uni_mapping.txt") %>%
  select(substrate = `yourlist:M201911225C475328CEF75220C360D524E9D456CE564CBDF`, uniprotID = `Entry name`) %>%
  group_by(substrate) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)

remove <- c("OTUD4_HUMAN", "ZO3_HUMAN", "MIA2_HUMAN", "ZF64B_HUMAN", "ARH38_HUMAN", "CCD61_HUMAN")
in_vv_substrates_toKIN <- ks %>%
  filter(source_type == "in_vivo") %>%
  select(substrate, position, residue) %>%
  distinct() %>%
  inner_join(in_vv_map, by = "substrate") %>%
  select(uniprotID, position, residue) %>%
  filter(!(uniprotID %in% remove))
write.table(in_vv_substrates_toKIN, "./output/files/in_vv_substrates_toKIN.txt", sep="\t", quote=F, row.names=F, col.names = F)

# import NetworKIN predictions
in_vvNKIN <- read_tsv(file = "./data/kinase_substrate/NKIN/networkin_predictions_invivo.tsv") %>%
  filter(networkin_score > 2) %>%
  select(kinase = id, substrate = `#substrate`, position) %>%
  #mutate(kinase = toupper(kinase)) %>%
  filter(!str_detect(kinase, "[a-z]")) %>%
  filter(!str_detect(kinase, "\\_")) %>%
  inner_join(in_vv_map, by = c("substrate" = "uniprotID")) %>%
  select(kinase, substrate.y, position) %>%
  rename(substrate = substrate.y)


in_vv_KS <- ks %>%
  filter(source_type == "in_vivo") %>%
  select(kinase, substrate, position, residue) %>%
  semi_join(in_vvNKIN, by = c("kinase", "substrate", "position"))
  


