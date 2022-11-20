library(tidyverse)


# load paper's functional score prediction
fun_score <- read_tsv("./data/protein/fun_score/media-3.txt")

# get uniprot to gene name mapping
mapping <- read_tsv("./data/protein/fun_score/uniprot_gene_name_mappings.tab")

#up <- UniProt.ws(taxId=9606)
#mapping <- UniProt.ws::select(up, unique(fun_score$uniprot), keytype = "UNIPROTKB", columns = c("GENECARDS","UNIPROTKB")) %>%
#  as_tibble() %>%
#  select(From=UNIPROTKB, To=GENECARDS) %>%
#  filter(!(is.na(From) | is.na(To)))

genes_remove <- mapping %>% group_by(To) %>% summarise(n = n()) %>% ungroup() %>% filter(n > 1)
proteins_remove <- mapping %>% group_by(From) %>% summarise(n = n()) %>% ungroup() %>% filter(n > 1)

mapping <- mapping %>%
  filter(!(From %in% proteins_remove$From | To %in% genes_remove$To))

# map uniprot accession IDs to gene name
fun_score <- fun_score %>%
  inner_join(mapping, by = c("uniprot" = "From")) %>%
  select(uniprot, gene=To, position, functional_score)

write_tsv(x = fun_score, path = "./output/files/fun_score_David_Paper.txt")

