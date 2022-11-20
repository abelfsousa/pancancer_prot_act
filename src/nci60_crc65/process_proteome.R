library(tidyverse)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# -- phosphorylation data


# load processed phosphosites
psites <- read_tsv("./output/files/nci60_crc65_processed_psites.txt")


# NCI60 cell lines
nci60_phospho <- read_tsv("./data/nci60_crc65/proteome/phospho/nci60_phospho.txt") %>%
  inner_join(psites, by = "ID") %>%
  select(-ID) %>%
  select(gene, position, residue, everything()) %>%
  pivot_longer(-c(gene, position, residue), names_to = "sample", values_to = "log10AU") %>%
  select(sample, gene, position, residue, log10AU) %>%
  mutate(AU = 10^log10AU) %>%
  select(-log10AU) %>%
  # take the mean of repeated psites
  group_by(gene, position, residue, sample) %>%
  summarise(AU = mean(AU, na.rm = T)) %>%
  ungroup() %>%
  mutate(AU = replace_na(AU, NA)) %>%
  # remove the psites with NAs in all samples
  group_by(gene, position, residue) %>%
  filter(!sum(is.na(AU)) == n()) %>%
  ungroup() %>%
  # calculate log2 fold-changes
  group_by(gene, position, residue) %>%
  mutate(log2fc = log2(AU/median(AU, na.rm = T))) %>%
  ungroup() %>%
  select(-AU) %>%
  mutate(batch = "NCI60") %>%
  select(batch, sample, everything())


# CRC65 cell lines
crc65_phospho <- read_tsv("./data/nci60_crc65/proteome/phospho/crc65_phospho.txt") %>%
  inner_join(psites, by = "ID") %>%
  select(-ID) %>%
  select(gene, position, residue, everything()) %>%
  pivot_longer(-c(gene, position, residue), names_to = "sample", values_to = "log10AU") %>%
  select(sample, gene, position, residue, log10AU) %>%
  # process cell line IDs
  mutate(sample = toupper(str_replace_all(sample, "-|/| ", ""))) %>%
  mutate(AU = 10^log10AU) %>%
  select(-log10AU) %>%
  # take the mean of repeated psites
  group_by(gene, position, residue, sample) %>%
  summarise(AU = mean(AU, na.rm = T)) %>%
  ungroup() %>%
  mutate(AU = replace_na(AU, NA)) %>%
  # remove the psites with NAs in all samples
  group_by(gene, position, residue) %>%
  filter(!sum(is.na(AU)) == n()) %>%
  ungroup() %>%
  # calculate log2 fold-changes
  group_by(gene, position, residue) %>%
  mutate(log2fc = log2(AU/median(AU, na.rm = T))) %>%
  ungroup() %>%
  select(-AU) %>%
  mutate(batch = "CRC65") %>%
  select(batch, sample, everything())

phospho <- bind_rows(nci60_phospho, crc65_phospho)
write_tsv(phospho, "./output/files/nci60_crc65_phospho_log2fc.txt.gz")

phospho_boxplot <- phospho %>%
  mutate(batch = fct_relevel(batch, "NCI60", "CRC65")) %>%
  ggplot(mapping = aes(x = sample, y = log2fc, fill = batch)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_wrap(~ batch, ncol = 2, scales = "free_x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    strip.text = element_text(size=16),
    legend.text = element_text(size=14),
    legend.title = element_text(size=16))

ggsave(filename = "nci60_crc65_phospho_boxplot.png", plot = phospho_boxplot, path = "./output/plots/nci60_crc65/", width = 12, height = 8)
unlink("nci60_crc65_phospho_boxplot.png")


# -- protein data


# load protein annotation
nci60_prot_annot <- read_tsv("./output/files/nci60_protein_annotation.txt")
crc65_prot_annot <- read_tsv("./output/files/crc65_protein_annotation.txt")


# NCI60 cell lines
nci60_protein <- read_tsv("./data/nci60_crc65/proteome/protein/nci60_protein.txt") %>%
  inner_join(nci60_prot_annot, by = "ID") %>%
  select(-ID) %>%
  select(gene, everything()) %>%
  pivot_longer(-gene, names_to = "sample", values_to = "log10iBAQ") %>%
  select(sample, gene, log10iBAQ) %>%
  mutate(iBAQ = 10^log10iBAQ) %>%
  select(-log10iBAQ) %>%
  # take the mean of repeated proteins
  group_by(gene, sample) %>%
  summarise(iBAQ = mean(iBAQ, na.rm = T)) %>%
  ungroup() %>%
  mutate(iBAQ = replace_na(iBAQ, NA)) %>%
  # remove the proteins with NAs in all samples
  group_by(gene) %>%
  filter(!sum(is.na(iBAQ)) == n()) %>%
  ungroup() %>%
  # calculate log2 fold-changes
  group_by(gene) %>%
  mutate(log2fc = log2(iBAQ/median(iBAQ, na.rm = T))) %>%
  ungroup() %>%
  select(-iBAQ) %>%
  mutate(batch = "NCI60") %>%
  select(batch, sample, everything())


# CRC65 cell lines
crc65_protein <- read_tsv("./data/nci60_crc65/proteome/protein/crc65_protein.txt") %>%
  inner_join(crc65_prot_annot, by = "ID") %>%
  select(-ID) %>%
  select(gene, everything()) %>%
  pivot_longer(-gene, names_to = "sample", values_to = "log10iBAQ") %>%
  select(sample, gene, log10iBAQ) %>%
  # process cell line IDs
  mutate(sample = toupper(str_replace_all(sample, "-|/| ", ""))) %>%
  mutate(iBAQ = 10^log10iBAQ) %>%
  select(-log10iBAQ) %>%
  # take the mean of repeated proteins
  group_by(gene, sample) %>%
  summarise(iBAQ = mean(iBAQ, na.rm = T)) %>%
  ungroup() %>%
  mutate(iBAQ = replace_na(iBAQ, NA)) %>%
  # remove the proteins with NAs in all samples
  group_by(gene) %>%
  filter(!sum(is.na(iBAQ)) == n()) %>%
  ungroup() %>%
  # calculate log2 fold-changes
  group_by(gene) %>%
  mutate(log2fc = log2(iBAQ/median(iBAQ, na.rm = T))) %>%
  ungroup() %>%
  select(-iBAQ) %>%
  mutate(batch = "CRC65") %>%
  select(batch, sample, everything())

protein <- bind_rows(nci60_protein, crc65_protein)
write_tsv(protein, "./output/files/nci60_crc65_protein_log2fc.txt.gz")

protein_boxplot <- protein %>%
  mutate(batch = fct_relevel(batch, "NCI60", "CRC65")) %>%
  ggplot(mapping = aes(x = sample, y = log2fc, fill = batch)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_wrap(~ batch, ncol = 2, scales = "free_x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    strip.text = element_text(size=16),
    legend.text = element_text(size=14),
    legend.title = element_text(size=16))

ggsave(filename = "nci60_crc65_protein_boxplot.png", plot = protein_boxplot, path = "./output/plots/nci60_crc65/", width = 12, height = 8)
unlink("nci60_crc65_protein_boxplot.png")


# -- correlate protein abundance and phosphorylation data

system.time(correlation <- phospho %>%
  rename(log2fc_phospho = log2fc) %>%
  inner_join(protein, by = c("batch", "sample", "gene")) %>%
  rename(log2fc_protein = log2fc) %>%
  filter(!(is.na(log2fc_phospho) | is.na(log2fc_protein))) %>%
  group_by(batch, gene, position, residue) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(cor = map_dbl(data, ~ cor.test(.x$log2fc_protein, .x$log2fc_phospho)$estimate)) %>%
  select(-data))

# considerably faster this way
system.time(correlation2 <- phospho %>%
  rename(log2fc_phospho = log2fc) %>%
  inner_join(protein, by = c("batch", "sample", "gene")) %>%
  rename(log2fc_protein = log2fc) %>%
  filter(!(is.na(log2fc_phospho) | is.na(log2fc_protein))) %>%
  group_by(batch, gene, position, residue) %>%
  filter(n() >= 10) %>%
  summarise(cor = cor(log2fc_protein, log2fc_phospho)) %>%
  ungroup())

correlation_boxplot <- correlation %>%
  mutate(batch = fct_relevel(batch, "NCI60", "CRC65")) %>%
  ggplot(mapping = aes(x = batch, y = cor, fill = batch)) +
  geom_boxplot(outlier.size = 0.5)

ggsave(filename = "nci60_crc65_protein_phospho_cor_boxplot.png", plot = correlation_boxplot, path = "./output/plots/nci60_crc65/", width = 3, height = 4)
unlink("nci60_crc65_protein_phospho_cor_boxplot.png")


# -- regress-out protein abundance from phosphorylation

source("./src/utils/regress_outCovs.R")

phospho_regout <- phospho %>%
  rename(log2fc_phospho = log2fc) %>%
  inner_join(protein, by = c("batch", "sample", "gene")) %>%
  rename(log2fc_protein = log2fc) %>%
  group_by(batch, gene, position, residue) %>%
  nest() %>%
  ungroup() %>%
  mutate(residual = map(data, regress_outCovs)) %>%
  unnest() %>%
  select(-starts_with("log2")) %>%
  select(batch, sample, gene, position, residue, log2fc = residual)

write_tsv(phospho_regout, "./output/files/nci60_crc65_phospho_log2fc_ProtRegOut.txt.gz")

phospho_boxplot <- phospho_regout %>%
  mutate(batch = fct_relevel(batch, "NCI60", "CRC65")) %>%
  ggplot(mapping = aes(x = sample, y = log2fc, fill = batch)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_wrap(~ batch, ncol = 2, scales = "free_x") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=14),
    axis.title = element_text(size=16),
    strip.text = element_text(size=16),
    legend.text = element_text(size=14),
    legend.title = element_text(size=16))

ggsave(filename = "nci60_crc65_phospho_boxplot_ProtRegOut.png", plot = phospho_boxplot, path = "./output/plots/nci60_crc65/", width = 12, height = 8)
unlink("nci60_crc65_phospho_boxplot_ProtRegOut.png")


# -- correlate protein abundance and protein regressed-out phosphorylation data

correlation <- phospho_regout %>%
  rename(log2fc_phospho = log2fc) %>%
  inner_join(protein, by = c("batch", "sample", "gene")) %>%
  rename(log2fc_protein = log2fc) %>%
  filter(!(is.na(log2fc_phospho) | is.na(log2fc_protein))) %>%
  group_by(batch, gene, position, residue) %>%
  filter(n() >= 10) %>%
  summarise(cor = cor(log2fc_protein, log2fc_phospho)) %>%
  ungroup()

correlation_boxplot <- correlation %>%
  mutate(batch = fct_relevel(batch, "NCI60", "CRC65")) %>%
  ggplot(mapping = aes(x = batch, y = cor, fill = batch)) +
  geom_boxplot(outlier.size = 0.5)

ggsave(filename = "nci60_crc65_protein_phospho_cor_boxplot_ProtRegOut.png", plot = correlation_boxplot, path = "./output/plots/nci60_crc65/", width = 4, height = 4)
unlink("nci60_crc65_protein_phospho_cor_boxplot_ProtRegOut.png")


# -- cell lines

nci60 <- read_tsv("./data/nci60_crc65/metadata/nci60_cell_lines.txt") %>%
  select(cell_line = Name, alternative_names = `DTP name`)

crc65 <- read_tsv("./data/nci60_crc65/metadata/crc65_cell_lines.txt") %>%
  select(cell_line = Name, alternative_names = `Alternative name`) %>%
  mutate(cell_line = toupper(str_replace_all(cell_line, "-|/| ", ""))) %>%
  mutate(alternative_names = str_split(alternative_names, "; ")) %>%
  unnest() %>%
  mutate(alternative_names = toupper(str_replace_all(alternative_names, "-|/| ", ""))) %>%
  distinct() %>%
  group_by(cell_line) %>%
  summarise(alternative_names = str_c(alternative_names, collapse = ";")) %>%
  ungroup()

cell_lines <- phospho_regout %>%
  select(batch, sample) %>%
  distinct() %>%
  rename(cell_line = sample) %>%
  left_join(crc65, by = "cell_line")

write_tsv(cell_lines, "./output/files/nci60_crc65_cell_lines.txt")
