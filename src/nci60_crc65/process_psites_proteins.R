library(tidyverse)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# -- phosphosites


# load UniProtKB/Swiss-Prot canonical protein sequences
seqs <- Biostrings::readAAStringSet(filepath = "./data/nci60_crc65/proteome/phospho/uniprot_sprot.fasta.gz")
seqs <- tibble(ids = names(seqs), seq = BiocGenerics::paste(seqs))
seqs <- seqs %>%
  filter(str_detect(seqs$ids, "OS=Homo sapiens OX=9606")) %>%
  mutate(protein = str_extract(ids, "sp\\|.*\\|")) %>%
  mutate(protein = str_replace_all(protein, "sp\\||\\|", "")) %>%
  select(protein, seq)


# load UniProtKB/Swiss-Prot non-canonical protein sequences (isoforms)
iso <- Biostrings::readAAStringSet(filepath = "./data/nci60_crc65/proteome/phospho/uniprot_sprot_varsplic.fasta.gz")
iso <- tibble(ids = names(iso), seq = BiocGenerics::paste(iso))
iso <- iso %>%
  filter(str_detect(iso$ids, "OS=Homo sapiens OX=9606")) %>%
  mutate(protein = str_extract(ids, "sp\\|.*\\|")) %>%
  mutate(protein = str_replace_all(protein, "sp\\||\\|", "")) %>%
  select(protein, seq)

seqs <- bind_rows(seqs, iso)


# load phosphosites annotation
psites1 <- data.table::fread("./data/nci60_crc65/proteome/phospho/phosphosites_annotation.txt") %>%
  as_tibble() %>%
  select(ID, gene = `Gene names`, protein = Protein, position = Position) %>%
  filter(str_detect(gene, "-Sep", negate = F))

map <- psites1 %>%
  select(protein) %>%
  distinct() %>%
  mutate(gene = c("SEPTIN4", "SEPTIN6", "SEPTIN2", "SEPTIN7", "SEPTIN12", "SEPTIN1", "SEPTIN8", "SEPTIN5", "SEPTIN5", "SEPTIN11", "SEPTIN10", "SEPTIN9", "SEPTIN9"))

psites1 <- psites1 %>%
  inner_join(map, by = "protein") %>%
  select(-gene.x, gene = gene.y)

psites2 <- data.table::fread("./data/nci60_crc65/proteome/phospho/phosphosites_annotation.txt") %>%
  as_tibble() %>%
  select(ID, gene = `Gene names`, protein = Protein, position = Position) %>%
  filter(str_detect(gene, "-Mar", negate = F))

map <- psites2 %>%
  select(protein) %>%
  distinct() %>%
  mutate(gene = c("MARCHF8", "MARCHF10", "MARCHF7", "MTARC2"))

psites2 <- psites2 %>%
  inner_join(map, by = "protein") %>%
  select(-gene.x, gene = gene.y)

psites3 <- data.table::fread("./data/nci60_crc65/proteome/phospho/phosphosites_annotation.txt") %>%
  as_tibble() %>%
  select(ID, gene = `Gene names`, protein = Protein, position = Position) %>%
  filter(str_detect(gene, "-Sep|-Mar", negate = T))

psites <- bind_rows(psites1, psites2, psites3) %>%
  left_join(seqs, by = "protein") %>%
  mutate(residue = map2_chr(.x = seq, .y = position, .f = ~ if(is.na(.x)){NA}else{str_sub(.x, .y, .y)})) %>%
  filter(residue %in% c("S", "T", "Y")) %>%
  select(-protein, -seq) %>%
  mutate(gene = str_split(gene, ";")) %>%
  unnest() %>%
  filter(!gene == "")

write_tsv(psites, "./output/files/nci60_crc65_processed_psites.txt")



# -- proteins

# load NCI60 protein annotation
nci60_proteins <- data.table::fread("./data/nci60_crc65/proteome/protein/nci60_proteins_annotation.txt") %>%
  as_tibble() %>%
  select(ID, gene = `Gene names`) %>%
  filter(!gene == "") %>%
  mutate(gene = str_split(gene, ";")) %>%
  unnest()

write_tsv(nci60_proteins, "./output/files/nci60_protein_annotation.txt")

# load CRC65 protein annotation
crc65_proteins <- data.table::fread("./data/nci60_crc65/proteome/protein/crc65_proteins_annotation.txt") %>%
  as_tibble() %>%
  select(ID, gene = `Gene names`) %>%
  filter(!gene == "") %>%
  mutate(gene = str_split(gene, ";")) %>%
  unnest()

write_tsv(crc65_proteins, "./output/files/crc65_protein_annotation.txt")

