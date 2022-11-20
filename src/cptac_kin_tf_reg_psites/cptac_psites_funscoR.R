# get CPTAC phosphosites functional score

library(funscoR)
library(knitr)
library(tidyverse)
library(viridis)

# load kinase-substrate lists
# select kinases
kin_sub <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz")
kinases <- kin_sub %>%
  select(kinase) %>%
  distinct() %>%
  pull(kinase)


# load CPTAC phosphosites
cptac_psites <- read_tsv("./output/files/phosphoproteomicsQ.txt.gz") %>%
  select(c(1:3)) %>%
  filter(psites == 1) %>%
  select(-psites) %>%
  mutate(residue = str_extract(psite, "_[A-Z]{1}"), position = str_extract(psite, "_[A-Z]{1}[0-9]+")) %>%
  mutate(residue = str_replace(residue, "_", ""), position = str_replace(position, "_[A-Z]{1}", "")) %>%
  select(-psite)

# get uniprot accession to gene name mappings
up <- UniProt.ws::UniProt.ws(taxId=9606)
mappings <- UniProt.ws::select(up, unique(cptac_psites$gene), keytype = "GENECARDS", columns = c("GENECARDS","UNIPROTKB"))
mappings2 <- mappings %>%
  as_tibble() %>%
  na.exclude() %>%
  group_by(GENECARDS) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 1) %>%
  select(-n)

cptac_psites <- cptac_psites %>%
  inner_join(mappings2, by = c("gene" = "GENECARDS")) %>%
  select(gene, acc=UNIPROTKB, position, residue)


# select kinases from CPTAC phosphosites
cptac_kin <- cptac_psites %>%
  filter(gene %in% kinases)


# load PhosphoSitePlus (PSP) regulatory sites
psp_reg <- read_tsv(file = "./data/PhosphoSitePlus/Regulatory_sites.gz", skip = 3) %>%
  filter(ORGANISM == "human") %>%
  select(gene=GENE, protein=PROTEIN, acc=ACC_ID, psite=MOD_RSD) %>%
  na.exclude() %>%
  separate(col = psite, into = c("psite", "modification")) %>%
  filter(modification == "p") %>%
  select(-modification) %>%
  mutate(residue = str_extract(psite, "^[A-Z]{1}"), position = str_extract(psite, "[0-9]+")) %>%
  select(-psite)


# overlap and select the cptac kinase psites with the PSP regulatory psites
# the psites that do not overlap are unknown regarding their regulatory function and will be predicted
cptac_kin_reg <- cptac_kin %>%
  semi_join(psp_reg, by = c("gene", "position", "residue"))

cptac_kin_unknown <- cptac_kin %>%
  setdiff(cptac_kin_reg)


# set up the funscoR model (Functional scoring of human phosphosites)
# https://evocellnet.github.io/funscoR/

## training the model

# annotate CPTAC psites with features
psites <- annotate_sites(cptac_psites[-c(1)])

# preprocess features for training
ST_features <- preprocess_features(psites, "ST")
Y_features <- preprocess_features(psites, "Y")

# train new model
# Gradient Boosting Machine by default
# use the most recent PSP regulatory sites as Gold Standard
psp2 <- psp_reg %>%
  select(acc, position) %>%
  as.data.frame()

ST_model <- train_funscore(ST_features, "ST", psp2, ncores = 3)
Y_model <- train_funscore(Y_features, "Y", psp2, ncores = 3)

## predicting CPTAC psites funscoR
ST_scores <- predict_funscore(ST_features, ST_model, ncores = 3)
Y_scores <- predict_funscore(Y_features, Y_model, ncores = 3)

## gather all predictions
all_scores <- bind_rows(ST_scores, Y_scores) %>%
  as_tibble() %>%
  mutate(probabilities = log_scaling(probabilities)) %>%
  separate(sites, c("acc", "position"), "_") %>%
  inner_join(cptac_psites, by = c("acc", "position")) %>%
  select(gene, acc, position, residue, probabilities)

# export all scores
write_tsv(x = all_scores, path = "./output/files/cptac_psites_funscoR.txt.gz")

cptac_kin_unknown <- all_scores %>%
  inner_join(cptac_kin_unknown) %>%
  select(gene, position, residue, probabilities) %>%
  mutate(PSP_reg_status = "unknown")

cptac_kin_reg <- all_scores %>%
  inner_join(cptac_kin_reg) %>%
  select(gene, position, residue, probabilities) %>%
  mutate(PSP_reg_status = "known")

cptac_kin <- bind_rows(cptac_kin_unknown, cptac_kin_reg)

# export kinases scores
write_tsv(x = cptac_kin, path = "./output/files/cptac_kin_psites_funscoR.txt.gz")

# plot probability distribution of known and unknown PSP kinase psites
cptac_kin_plot <- cptac_kin %>%
  ggplot(mapping = aes(x = probabilities, fill = PSP_reg_status, color = PSP_reg_status)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  scale_fill_viridis(discrete = T, option = "C") +
  scale_color_viridis(discrete = T, option = "C") +
  theme(
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    legend.title = element_text(colour = "black", size = 10),
    legend.text = element_text(colour = "black", size = 8))

ggsave(filename = "cptac_kin_reg_unknown_sites.png", plot = cptac_kin_plot, path = "./output/plots/cptac_kin_regulatory_sites/", height = 2, width = 6)
unlink("cptac_kin_reg_unknown_sites.png")
