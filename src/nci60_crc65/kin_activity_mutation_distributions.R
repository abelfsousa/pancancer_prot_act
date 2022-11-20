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


#' # -- assess the kinase activity distributions between different mutation types
mut1 <- nci60_crc65_mutation %>%
  select(sample, gene, mutation_type) %>%
  distinct()

mut2 <- roum_mutation %>%
  select(sample, gene=gene_symbol, mutation_type=variant_class) %>%
  filter(!mutation_type == "Splice_Site") %>%
  mutate(mutation_type = str_replace(mutation_type, "Missense_Mutation", "Missense")) %>%
  mutate(mutation_type = str_replace(mutation_type, "Nonsense_Mutation", "Nonsense")) %>%
  mutate(mutation_type = str_replace(mutation_type, "Frame_Shift_Ins|Frame_Shift_Del", "Frameshift")) %>%
  mutate(mutation_type = str_replace(mutation_type, "In_Frame_Del|In_Frame_Ins", "Nonframeshift")) %>%
  mutate(mutation_type = str_replace(mutation_type, "Nonstop_Mutation", "Read_through")) %>%
  distinct()


#' ## using only NCI60/CRC65 datasets
mutations <- mut1 %>%
  semi_join(kin_activity, by = c("gene" = "kinase", "sample")) %>%
  filter(!mutation_type %in% c("Nonframeshift", "Read_through")) %>%
  group_by(mutation_type) %>%
  nest() %>%
  ungroup() #%>%
  #mutate(data2 = (function(x) map2(.x=1:length(x),.y=x,.f=function(a,b)reduce(.x=c(list(b),x[-a]),.f=setdiff))) (data))

mutations2 <- mutations %>%
  mutate(data2 = map2(.x = mutation_type, .y = data, .f = function(x,y) reduce(.x=c(list(y), mutations[mutations$mutation_type != x, "data", drop = T]), .f=setdiff))) %>%
  select(-data) %>%
  unnest()

mutations <- mutations %>%
  unnest()

distribution <- kin_activity %>%
  inner_join(mutations2, by = c("sample", "kinase" = "gene")) %>%
  ggplot(mapping = aes(x=mutation_type, y = kin_activity, fill = mutation_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, show.legend = F) +
  geom_jitter(size = 1, width = 0.1, alpha = 0.1, show.legend = F) +
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(3,4), c(2,4),c(1,4))) +
  coord_flip()

#+ fig.width=7, fig.height=3
distribution

distribution <- kin_activity %>%
  inner_join(mutations2, by = c("sample", "kinase" = "gene")) %>%
  ggplot(mapping = aes(x=kin_activity, fill = mutation_type, color = mutation_type)) +
  geom_density(alpha = 0.5)

#+ fig.width=7, fig.height=2
distribution

frameshift <- kin_activity %>%
  inner_join(mutations2, by = c("sample", "kinase" = "gene")) %>%
  filter(mutation_type == "Frameshift") %>%
  arrange((kin_activity))


#' ## using NCI60/CRC65 datasets + Roumeliotis
mutations <- bind_rows(mut1, mut2) %>%
  semi_join(kin_activity, by = c("gene" = "kinase", "sample")) %>%
  filter(!mutation_type %in% c("Nonframeshift", "Read_through")) %>%
  group_by(mutation_type) %>%
  nest() %>%
  ungroup()

mutations2 <- mutations %>%
  mutate(data2 = map2(.x = mutation_type, .y = data, .f = function(x,y) reduce(.x=c(list(y), mutations[mutations$mutation_type != x, "data", drop = T]), .f=setdiff))) %>%
  select(-data) %>%
  unnest()

mutations <- mutations %>%
  unnest()

distribution <- kin_activity %>%
  inner_join(mutations2, by = c("sample", "kinase" = "gene")) %>%
  ggplot(mapping = aes(x=mutation_type, y = kin_activity, fill = mutation_type)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, show.legend = F) +
  geom_jitter(size = 2, width = 0.1, alpha = 0.1, show.legend = F) +
  ggpubr::stat_compare_means(comparisons = list(c(1,2),c(3,4), c(2,4),c(1,4))) +
  coord_flip()

#+ fig.width=7, fig.height=3
distribution

distribution <- kin_activity %>%
  inner_join(mutations2, by = c("sample", "kinase" = "gene")) %>%
  ggplot(mapping = aes(x=kin_activity, fill = mutation_type, color = mutation_type)) +
  geom_density(alpha = 0.5)

#+ fig.width=7, fig.height=2
distribution


#' load phosphorylation data
phospho <- read_tsv(file = "./output/files/nci60_crc65_phospho_log2fc_ProtRegOut.txt.gz") %>%
  filter(!(batch == "CRC65" & sample %in% shared_cell_lines))

#' load kinase-substrate lists
ks_lists <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz") %>%
  select(-pair, -source) %>%
  distinct() %>%
  mutate(source_type = str_replace(source_type, "database|text-mining", "DB_text-mining")) %>%
  distinct() %>%
  filter(source_type %in% "DB_text-mining") %>%
  select(-source_type)


#' # -- Compare PRKCB phosphorylation targets between HCT116 and other cell lines\
plot <- phospho %>%
  filter(!is.na(log2fc)) %>%
  inner_join(ks_lists[ks_lists$kinase == "PRKCB", -c(1)], by = c("gene" = "substrate", "position", "residue")) %>%
  mutate(is_sample = if_else(sample == "HCT116", "yes", "no")) %>%
  mutate(sample = fct_reorder(sample, log2fc, median)) %>%
  ggplot(mapping = aes(x = sample, y = log2fc, fill = is_sample)) +
  geom_boxplot(outlier.size = 1) +
  scale_fill_manual(name = "Sample", values = c("grey", "red"), labels = c("others", "HCT116")) +
  scale_x_discrete(name = "") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#+ fig.width=10, fig.height=4
plot

plot <- phospho %>%
  filter(!is.na(log2fc)) %>%
  inner_join(ks_lists[ks_lists$kinase == "PRKCB", -c(1)], by = c("gene" = "substrate", "position", "residue")) %>%
  mutate(is_sample = if_else(sample == "HCT116", "yes", "no")) %>%
  #group_by(is_sample, gene, position, residue) %>%
  #summarise(log2fc = median(log2fc, na.rm = T)) %>%
  #ungroup() %>%
  mutate(is_sample = fct_reorder(is_sample, log2fc, median)) %>%
  ggplot(mapping = aes(x = is_sample, y = log2fc, fill = is_sample)) +
  geom_boxplot() +
  ggpubr::stat_compare_means() +
  scale_fill_manual(name = "Sample", values = c("red", "grey"), labels = c("HCT116", "others")) +
  scale_x_discrete(name = "") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#+ fig.width=4, fig.height=4
plot


#' # -- Compare ARAF phosphorylation targets between SW48 and other cell lines\
plot <- phospho %>%
  filter(!is.na(log2fc)) %>%
  inner_join(ks_lists[ks_lists$kinase == "ARAF", -c(1)], by = c("gene" = "substrate", "position", "residue")) %>%
  mutate(is_sample = if_else(sample == "SW48", "yes", "no")) %>%
  mutate(sample = fct_reorder(sample, log2fc, median)) %>%
  ggplot(mapping = aes(x = sample, y = log2fc, fill = is_sample)) +
  geom_boxplot(outlier.size = 1) +
  scale_fill_manual(name = "Sample", values = c("grey", "red"), labels = c("others", "SW48")) +
  scale_x_discrete(name = "") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#+ fig.width=10, fig.height=4
plot

plot <- phospho %>%
  filter(!is.na(log2fc)) %>%
  inner_join(ks_lists[ks_lists$kinase == "ARAF", -c(1)], by = c("gene" = "substrate", "position", "residue")) %>%
  mutate(is_sample = if_else(sample == "SW48", "yes", "no")) %>%
  #group_by(is_sample, gene, position, residue) %>%
  #summarise(log2fc = median(log2fc, na.rm = T)) %>%
  #ungroup() %>%
  mutate(is_sample = fct_reorder(is_sample, log2fc, median)) %>%
  ggplot(mapping = aes(x = is_sample, y = log2fc, fill = is_sample)) +
  geom_boxplot() +
  ggpubr::stat_compare_means() +
  scale_fill_manual(name = "Sample", values = c("red", "grey"), labels = c("SW48", "others")) +
  scale_x_discrete(name = "") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#+ fig.width=4, fig.height=4
plot

