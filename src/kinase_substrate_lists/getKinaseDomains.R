# get the kinase domain position for each kinase in the kinase-substrate list


# load R packages
library(tidyverse)

source("./src/utils/matchTool.R")

set.seed(123)

# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load AA code
aa <- read_tsv("./data/protein/protein_sequences/aminoacids_code.txt", skip = 1)


# load protein sequences for the canonical uniprot transcripts
prot_seq <- read_tsv(file = "./output/files/uniprot_canTrpt_protSeqGRCh37.txt")


# load KinBase kinase domain sequences
kin_domains <- Biostrings::readAAStringSet(filepath = "./data/protein/protein_sequences/KinBase_kinase_domains.txt")
kin_domains <- tibble(kinase = names(kin_domains), domain = BiocGenerics::paste(kin_domains)) %>%
  mutate(kinase = str_extract(kinase, "gene=[A-Za-z0-9]+"), kinase = str_replace(kinase, "gene=", "")) %>%
  group_by(kinase) %>%
  summarise(domain = list(domain)) %>%
  ungroup()



align_domain <- function(protSeq, kinDomain){
  s <- protSeq
  p <- kinDomain
  
  pos <- as.numeric(str_locate(s, p))
  
  if(sum(!is.na(pos)) > 0){
    res <- tibble(start = pos[1], end = pos[2])
  } else {
    res <- tibble(start = NA, end = NA)
    
    sl <- nchar(s)
    pl <- nchar(p)
    N = (sl - pl) + 1
    
    # alPos <- tibble(pos = 1:N)
    #
    # bestAl <- alPos %>%
    #   mutate(
    #     score = map_dbl(
    #       .x = pos,
    #       .f = ~ Biostrings::pairwiseAlignment(s, str_sub(s, .x, .x+pl-1), scoreOnly = TRUE))) %>%
    #   pull(score) %>%
    #   mean()
    # 
    # worstAl <- alPos %>%
    #   mutate(
    #     score = map_dbl(
    #       .x = pos,
    #       .f = ~ Biostrings::pairwiseAlignment(s, str_c(sample(aa[[1]], pl, replace = T), collapse = ""), scoreOnly = TRUE))) %>%
    #   pull(score) %>%
    #   mean()
    
    bestSeq <- str_sub(s, 1, 1+pl-1)
    worstSeq <- str_c(sample(aa[[1]], pl, replace = T), collapse = "")
    
    bestAl <- Biostrings::pairwiseAlignment(s, bestSeq, scoreOnly = TRUE)
    worstAl <- Biostrings::pairwiseAlignment(s, worstSeq, scoreOnly = TRUE)
    
    al <- Biostrings::pairwiseAlignment(p, s)
    
    norm_score <- (al@score-worstAl)/(bestAl-worstAl)
    
    if(norm_score > 0.5){
      al_seq <- as.character(al@subject)
      pos <- as.numeric(str_locate(s, al_seq))
      
      res <- tibble(start = pos[1], end = pos[2])
    }
  }
  res
}



# load kinase substrate list
# get respective protein sequence
# get respective kinase domain
kin_list <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz") %>%
  select(kinase) %>%
  distinct() %>%
  inner_join(prot_seq[, c("gene_name", "seq")], by = c("kinase" = "gene_name")) %>%
  rename(sequence = seq) %>%
  inner_join(kin_domains, by = "kinase") %>%
  arrange(desc(map_dbl(domain, ~ length(.x)))) %>%
  unnest() %>%
  mutate(pos = map2(.x = sequence, .y = domain, .f = align_domain)) %>%
  unnest() %>%
  filter(!is.na(start)) %>%
  select(-sequence, -domain)
  
write_tsv(x = kin_list, path = "./output/files/kinase_domain_positions.txt")
