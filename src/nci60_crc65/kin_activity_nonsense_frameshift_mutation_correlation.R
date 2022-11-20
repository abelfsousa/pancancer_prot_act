#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)

#' load R packages
library(tidyverse)

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load CPTAC/CCLE samples
cptac_samples <- read_tsv(file = "./output/files/all_samples_annotation.txt") %>%
  filter(data == "phosphorylation")


#' load NCI60/CRC65 cell lines
nci60_crc65 <- read_tsv(file = "./output/files/nci60_crc65_cell_lines.txt")

shared_cell_lines <- nci60_crc65 %>%
  group_by(batch) %>%
  summarise(cell_line = list(cell_line)) %>%
  pull(cell_line) %>%
  reduce(intersect)

alternative_names <- nci60_crc65 %>%
  filter(!is.na(alternative_names)) %>%
  mutate(alternative_names = str_split(alternative_names, ";")) %>%
  unnest() %>%
  distinct()


#' load kinase activity inference data from NCI60/CRC65 dataset\
#' select quantifications with more than 3 substrates
nci60_crc65_KA <- read_tsv(file = "./output/files/nci60_crc65_kinase_activities.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= 3) %>%
  select(-source_type, -n) %>%
  filter(!(batch == "CRC65" & sample %in% shared_cell_lines)) %>%
  select(sample, kinase, kin_activity=log10P)


#' load kinase activity inference data without imputations from CPATC/CCLE dataset\
#' (quantile-normalized protein regressed-out phosphorylation data)\
#' select quantifications with more than 3 substrates
ccle_KA <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= 3) %>%
  select(-source_type, -n) %>%
  select(sample, kinase, kin_activity=log10P) %>%
  inner_join(cptac_samples[, c("sample", "batch")], by = "sample") %>%
  filter(batch == "ccle") %>%
  select(-batch) %>%
  anti_join(nci60_crc65_KA, by = "sample") %>%
  anti_join(alternative_names, by = c("sample" = "alternative_names"))


#' bind the kinase activities of both datasets
kin_activity <- bind_rows(ccle_KA, nci60_crc65_KA)


#' load NCI60/CRC65 mutations
nci60_crc65_mutation <- read_tsv(file = "./output/files/nci60_crc65_mutations.txt.gz")


#' load Roumeliotis mutations
roum_mutation <- read_tsv(file = "./output/files/mutations_protpos.txt.gz") %>%
  inner_join(cptac_samples[, c("sample", "batch")], by = "sample") %>%
  filter(batch == "ccle") %>%
  select(-batch) %>%
  anti_join(nci60_crc65_KA, by = "sample") %>%
  anti_join(alternative_names, by = c("sample" = "alternative_names"))


#' load protein sequences
prot_seqs <- read_tsv(file = "./output/files/protein_sequences_GRCh37.txt") %>%
  filter(state == "known" & gene_biotype == "protein_coding" & trpt_biotype == "protein_coding") %>%
  select(trpt_id, protein_id, seq)


#' load canonical uniprot transcripts
canonical_trpt <- read_tsv(file = "./data/uniprot/canonical_trpt/isoform_overrides_uniprot.txt", skip = 1)


#' # -- assess the impact of Nonsense mutations on kinase activities
nonsense_nci60_crc65 <- nci60_crc65_mutation %>%
  filter(mutation_type == "Nonsense") %>%
  select(sample, gene, identifier, HGVSp, prot_pos_i)

nonsense_roum <- roum_mutation %>%
  filter(variant_class == "Nonsense_Mutation") %>%
  select(sample, gene=gene_symbol, CHROM, POS, REF, ALT, HGVSp, prot_pos_i = prot_pos) %>%
  unite(col = "identifier","CHROM", "POS", "REF", "ALT")

#' ## using only NCI60/CRC65 datasets
nonsense_impact <- nonsense_nci60_crc65 %>%
  inner_join(kin_activity, by = c("sample", "gene" = "kinase")) %>%
  #group_by(sample, gene) %>%
  #filter(n() == 1) %>%
  #ungroup() %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene" = "gene_name")) %>%
  inner_join(prot_seqs[, -2], by = c("enst_id" = "trpt_id")) %>%
  mutate(protei_len = nchar(seq), pos_p = (prot_pos_i/protei_len)*100) %>%
  pivot_longer(cols = c("pos_p", "prot_pos_i"), names_to = "position_type", values_to = "position") %>%
  ggplot(mapping = aes(x = position, y = kin_activity)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggpubr::stat_cor(method = "pearson") +
  facet_wrap(~ position_type, ncol = 2, scales = "free", labeller = labeller(position_type = c("pos_p"="%", "prot_pos_i"="aa")))

#+ fig.width=8, fig.height=4
nonsense_impact

#' ## using NCI60/CRC65 datasets + Roumeliotis
nonsense_impact <- nonsense_nci60_crc65 %>%
  bind_rows(nonsense_roum) %>%
  inner_join(kin_activity, by = c("sample", "gene" = "kinase")) %>%
  #group_by(sample, gene) %>%
  #filter(n() == 1) %>%
  #ungroup() %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene" = "gene_name")) %>%
  inner_join(prot_seqs[, -2], by = c("enst_id" = "trpt_id")) %>%
  mutate(protei_len = nchar(seq), pos_p = (prot_pos_i/protei_len)*100) %>%
  pivot_longer(cols = c("pos_p", "prot_pos_i"), names_to = "position_type", values_to = "position") %>%
  ggplot(mapping = aes(x = position, y = kin_activity)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggpubr::stat_cor(method = "pearson") +
  facet_wrap(~ position_type, ncol = 2, scales = "free", labeller = labeller(position_type = c("pos_p"="%", "prot_pos_i"="aa")))

#+ fig.width=8, fig.height=4
nonsense_impact


#' # -- assess the impact of Frameshift mutations on kinase activities
frameshift_nci60_crc65 <- nci60_crc65_mutation %>%
  filter(mutation_type == "Frameshift") %>%
  select(sample, gene, identifier, HGVSp, prot_pos_i)

frameshift_roum <- roum_mutation %>%
  filter(variant_class == "Frame_Shift_Del" | variant_class == "Frame_Shift_Ins") %>%
  select(sample, gene=gene_symbol, CHROM, POS, REF, ALT, HGVSp, prot_pos_i = prot_pos) %>%
  unite(col = "identifier","CHROM", "POS", "REF", "ALT")

#' ## using only NCI60/CRC65 datasets
frameshift_impact <- frameshift_nci60_crc65 %>%
  inner_join(kin_activity, by = c("sample", "gene" = "kinase")) %>%
  group_by(sample, gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene" = "gene_name")) %>%
  inner_join(prot_seqs[, -2], by = c("enst_id" = "trpt_id")) %>%
  mutate(protei_len = nchar(seq), pos_p = (prot_pos_i/protei_len)*100) %>%
  pivot_longer(cols = c("pos_p", "prot_pos_i"), names_to = "position_type", values_to = "position") %>%
  ggplot(mapping = aes(x = position, y = kin_activity)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggpubr::stat_cor(method = "pearson") +
  facet_wrap(~ position_type, ncol = 2, scales = "free", labeller = labeller(position_type = c("pos_p"="%", "prot_pos_i"="aa")))

#+ fig.width=8, fig.height=4
frameshift_impact


#' ## using NCI60/CRC65 datasets + Roumeliotis
frameshift_impact <- frameshift_nci60_crc65 %>%
  bind_rows(frameshift_roum) %>%
  inner_join(kin_activity, by = c("sample", "gene" = "kinase")) %>%
  group_by(sample, gene) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  inner_join(canonical_trpt[, c("gene_name", "enst_id")], by = c("gene" = "gene_name")) %>%
  inner_join(prot_seqs[, -2], by = c("enst_id" = "trpt_id")) %>%
  mutate(protei_len = nchar(seq), pos_p = (prot_pos_i/protei_len)*100) %>%
  pivot_longer(cols = c("pos_p", "prot_pos_i"), names_to = "position_type", values_to = "position") %>%
  ggplot(mapping = aes(x = position, y = kin_activity)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggpubr::stat_cor(method = "pearson") +
  facet_wrap(~ position_type, ncol = 2, scales = "free", labeller = labeller(position_type = c("pos_p"="%", "prot_pos_i"="aa")))

#+ fig.width=8, fig.height=4
frameshift_impact

