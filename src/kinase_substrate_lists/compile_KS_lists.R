#' ---
#' title: "Compile kinase-substrate lists"
#' author: "Abel Sousa"
#' ---


#+ setup, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)



#' load R packages
library(tidyverse)
source("./src/utils/psite_prot_seq_match.R")


#' load protein sequences for the canonical uniprot transcripts
gene_prot_seq <- read_tsv(file = "./output/files/uniprot_canTrpt_protSeqGRCh37.txt")


#' in vitro list from "large-scale discovery of substrates of the human kinome"
in_vitro_mapping <- read_tsv(file = "./data/kinase_substrate/Vivo_vitrodata/in_vitro_entrynames_genename.tab")

duplicates <- in_vitro_mapping %>%
  filter(duplicated(From)) %>%
  pull(From) %>%
  unique()

in_vitro_mapping <- in_vitro_mapping %>%
  filter(!From %in% duplicates)

in_vitro <- read_tsv(file = "./data/kinase_substrate/Vivo_vitrodata/41598_2019_46385_MOESM3_ESM.txt") %>%
  select(-type, -protein_description) %>%
  mutate(residue = str_extract(position, "[SYT]{1}")) %>%
  mutate(position = as.numeric(str_extract(position, "[0-9]+"))) %>%
  inner_join(in_vitro_mapping, by = c("substrate_uniprot_ID" = "From")) %>%
  select(-substrate_uniprot_ID, substrate = To) %>%
  mutate(kinase = str_split_fixed(kinase, "\\[", 2)[,1]) %>%
  mutate(kinase = str_split_fixed(kinase, "/", 2)[,1]) %>%
  mutate(kinase = toupper(kinase)) %>%
  mutate(pair = str_c(kinase, substrate, sep = "_")) %>%
  select(kinase, substrate, pair, everything()) %>%
  distinct() %>%
  mutate(source = "in_vitro") %>%
  inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("substrate" = "gene_name")) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, everything()) %>%
  group_by(n, kinase, substrate) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = map_psite2)) %>%
  unnest(cols = data) %>%
  filter(same == 1) %>%
  select(-n, -seq, -same)



#' in vivo list from Pedro's collaborator
in_vivo <- read_tsv(file = "./data/kinase_substrate/Vivo_vitrodata/41587_2019_391_MOESM5_ESM-3.txt") %>%
  distinct() %>%
  select(-X4, -`Kinase Family`) %>%
  rename(substrate = `Putative Downtream Target`) %>%
  mutate(position = str_extract(substrate, "\\([A-Za-z]+[0-9]+\\)")) %>%
  mutate(position = str_replace_all(position, "\\(|\\)", "")) %>%
  mutate(residue = str_extract(position, "^[A-Z]+")) %>%
  mutate(position = str_extract(position, "[0-9]+")) %>%
  mutate(substrate = str_replace(substrate, "\\([A-Za-z]+[0-9]+\\)", "")) %>%
  filter(!(is.na(residue) & is.na(position))) %>%
  mutate(position = as.numeric(position)) %>%
  mutate(pair = str_c(kinase, substrate, sep = "_")) %>%
  select(kinase, substrate, pair, everything()) %>%
  mutate(source = "in_vivo") %>%
  inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("substrate" = "gene_name")) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, everything()) %>%
  group_by(n, kinase, substrate) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = map_psite2)) %>%
  unnest(cols = data) %>%
  filter(same == 1) %>%
  select(-n, -seq, -same)


# in_vivo <- read_tsv(file = "./data/kinase_substrate/Vivo_vitrodata/55349_0_supp_454633_pjdby1.txt") %>%
#   distinct() %>%
#   mutate(substrate = str_replace(substrate, "\\(", "_")) %>%
#   mutate(substrate = str_replace(substrate, "\\)", "")) %>%
#   separate(col = substrate, into = c("substrate", "position"), sep = "_") %>%
#   mutate(residue = str_extract(position, "[SYT]{1}")) %>%
#   mutate(position = as.numeric(str_extract(position, "[0-9]+"))) %>%
#   filter(!is.na(position)) %>%
#   mutate(n = str_count(kinase, "\\."))
# 
# in_vivo2 <- in_vivo %>%
#   filter(n>1) %>%
#   mutate(kinase = str_replace(kinase, "\\.[0-9B]{1}$", "")) %>%
#   mutate(kinase = str_replace(kinase, "\\.[0-9B]{1}", "")) %>%
#   select(-n)
# 
# in_vivo <- in_vivo %>%
#   filter(n == 1) %>%
#   select(-n) %>%
#   bind_rows(in_vivo2) %>%
#   mutate(kinase = str_split(kinase, "\\.", 2)) %>%
#   unnest(cols = "kinase") %>%
#   mutate(pair = str_c(kinase, substrate, sep = "_")) %>%
#   select(kinase, substrate, pair, everything()) %>%
#   distinct() %>%
#   mutate(source = "in_vivo") %>%
#   inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("substrate" = "gene_name")) %>%
#   mutate(n = 1:nrow(.)) %>%
#   select(n, everything()) %>%
#   group_by(n, kinase, substrate) %>%
#   nest() %>%
#   ungroup() %>%
#   mutate(same = map_dbl(.x = data, .f = map_psite2)) %>%
#   unnest(cols = data) %>%
#   filter(same == 1) %>%
#   select(-n, -seq, -same)


#' join in vivo and in vitro pairs  
in_vitro_vivo <- bind_rows(
  in_vitro,
  in_vivo) %>%
  group_by(kinase, substrate, pair, position, residue) %>%
  summarise(source = str_c(source, collapse = "|"), source_n = n()) %>%
  ungroup()

in_vv_pairs <- intersect(in_vitro$pair,in_vivo$pair)



#' phosphosite plus kinase-substrate dataset
pplus <- data.table::fread(file = "./data/kinase_substrate/phosphositeplus/Kinase_Substrate_Dataset") %>%
  as_tibble() %>%
  filter(KIN_ORGANISM == "human" & SUB_ORGANISM == "human") %>%
  select(kinase = GENE, substrate = SUB_GENE, position = SUB_MOD_RSD) %>%
  distinct() %>%
  mutate(residue = str_extract(position, "[SYT]{1}")) %>%
  mutate(position = as.numeric(str_extract(position, "[0-9]+"))) %>%
  mutate(pair = str_c(kinase, substrate, sep = "_")) %>%
  select(kinase, substrate, pair, everything()) %>%
  mutate(source = "pplus") %>%
  inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("substrate" = "gene_name")) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, everything()) %>%
  group_by(n, kinase, substrate) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = map_psite2)) %>%
  unnest(cols = data) %>%
  filter(same == 1) %>%
  select(-n, -seq, -same)


#' join in vivo/vitro and pplus pairs
in_vv_pplus <- bind_rows(
  in_vitro,
  in_vivo,
  pplus) %>%
  group_by(kinase, substrate, pair, position, residue) %>%
  summarise(source = str_c(source, collapse = "|"), source_n = n()) %>%
  ungroup()

#dplyr::intersect(in_vitro %>% filter(pair %in% in_vv_pairs) %>% select(-source), pplus %>% select(-source))
#dplyr::intersect(in_vivo %>% filter(pair %in% in_vv_pairs) %>% select(-source), pplus %>% select(-source))



#' protmapper phosphorylation data
protmapper <- read_csv(file = "./data/kinase_substrate/protmapper/export.csv") %>%
  filter(CTRL_IS_KINASE) %>%
  select(kinase=CTRL_GENE_NAME, substrate=TARGET_GENE_NAME, position=TARGET_POS, residue=TARGET_RES, source=SOURCES) %>%
  filter(!is.na(kinase)) %>%
  mutate(source = str_split(source, ",")) %>%
  unnest(cols = source) %>%
  distinct() %>%
  mutate(pair = str_c(kinase, substrate, sep = "_")) %>%
  select(kinase, substrate, pair, everything()) %>%
  inner_join(gene_prot_seq[, c("gene_name", "seq")], by = c("substrate" = "gene_name")) %>%
  mutate(n = 1:nrow(.)) %>%
  select(n, everything()) %>%
  group_by(n, kinase, substrate) %>%
  nest() %>%
  ungroup() %>%
  mutate(same = map_dbl(.x = data, .f = map_psite2)) %>%
  unnest(cols = data) %>%
  filter(same == 1) %>%
  select(-n, -seq, -same)




#' join in vivo/vitro and protmapper pairs
in_vv_pmapper <- bind_rows(
  in_vitro,
  in_vivo,
  protmapper)

source_type <-
  in_vv_pmapper %>%
  group_by(source) %>%
  tally() %>%
  arrange(desc(n)) %>%
  select(-n) %>%
  bind_cols(
    source_type = c(
      "in_vitro",
      "in_vivo",
      "database",
      "database",
      "text-mining",
      "text-mining",
      "text-mining",
      "database",
      "database",
      "database"))

in_vv_pmapper <- in_vv_pmapper %>%
  inner_join(source_type, by = "source")

head(in_vv_pmapper)
dim(in_vv_pmapper)


gz1 <- gzfile("./output/files/kinase_substrate_list.txt.gz", "w")
write.table(in_vv_pmapper, gz1, sep="\t", quote=F, row.names=F)
close(gz1)
