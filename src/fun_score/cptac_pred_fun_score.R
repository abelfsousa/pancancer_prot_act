library(UniProt.ws)
library(tidyverse)
library(funscoR)
library(knitr)


# load paper's functional score prediction
fun_score <- read_tsv("./data/protein/fun_score/media-3.txt") %>%
  mutate(psite = str_c(uniprot, position, sep = "_"))


# load CPTAC phosphosites
cptac <- read_tsv("./output/files/phosphoproteomicsQ.txt.gz")

cptac_pho <- cptac %>%
  dplyr::select(c(1:3)) %>%
  filter(psites == 1) %>%
  dplyr::select(-psites) %>%
  mutate(residue = str_extract(psite, "_[A-Z]{1}"), position = str_extract(psite, "_[A-Z]{1}[0-9]+")) %>%
  mutate(residue = str_replace(residue, "_", ""), position = str_replace(position, "_[A-Z]{1}", "")) %>%
  dplyr::select(-psite)



# get gene name - uniprot accession mappings
up <- UniProt.ws(taxId=9606)
mappings <- UniProt.ws::select(up, unique(cptac_pho$gene),
       keytype = "GENECARDS",
       columns = c("GENECARDS","UNIPROTKB"))

mappings2 <- mappings %>%
  as_tibble() %>%
  filter(!(is.na(UNIPROTKB) | is.na(GENECARDS))) %>%
  group_by(GENECARDS) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  dplyr::select(-n)


# map gene name to uniprot accession IDs
cptac_pho <- cptac_pho %>%
  inner_join(mappings2, by = c("gene" = "GENECARDS")) %>%
  dplyr::select(-gene) %>%
  dplyr::select(acc = UNIPROTKB, position, residue) %>%
  mutate(psite = str_c(acc, position, sep = "_"))


# annotate phosphoproteome with features
annotated_phos <- annotate_sites(cptac_pho[,-c(4)])


# preprocess features for training
ST_features <- preprocess_features(annotated_phos, "ST")
Y_features <- preprocess_features(annotated_phos, "Y")


# train new model
ST_model <- train_funscore(ST_features, "ST", psp, ncores = 1)
Y_model <- train_funscore(Y_features, "Y", psp, ncores = 1)


# predict funcscoR for all sites
ST_scores <- predict_funscore(ST_features, ST_model, ncores = 1)
Y_scores <- predict_funscore(Y_features, Y_model, ncores = 1)


# gather all predictions
all_scores <- bind_rows(ST_scores, Y_scores) %>%
  mutate(probabilities = log_scaling(probabilities)) %>%
  mutate(in_paper = if_else(sites %in% fun_score$psite, "1", "0"))

all_scores2 <- all_scores %>%
  inner_join(fun_score[, c("psite", "functional_score")], by = c("sites" = "psite"))

plot(
  all_scores2$functional_score,
  all_scores2$probabilities,
  xlab = "paper functional score",
  ylab = "model functional score")


# plot density distribution between in paper and not in paper
density_plot <- ggplot(data = all_scores) +
  geom_density(mapping = aes(x = probabilities, fill = in_paper, color = in_paper), alpha = 0.5)
density_plot

