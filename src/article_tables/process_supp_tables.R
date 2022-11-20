library(tidyverse)


# load samples metadata
source("./src/utils/getSamples.R")
samples_metadata <- read_tsv("./output/files/all_samples_annotation.txt")
samples_metadata <- getSamples(samples_metadata, data_types = c("protein")) %>%
  mutate(batch = map2_chr(.x = batch, .y = tissue, .f = ~ if(.x == "ccle"){str_c(.x, .y, sep = "-")}else{.x})) %>%
  select(sample, dataset=batch, tissue) 

write_tsv(samples_metadata, "./output/files/samples_metadata_article.txt")


# load kinase activity data
k_subN <- 3
kin_act <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  #kin_act <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf_notKinSites.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-n, -source_type) %>%
  rename(kin_activity = log10P) %>%
  select(kinase, everything()) %>%
  pivot_wider(names_from = "sample", values_from = "kin_activity")

# kin_act <- kin_act %>%
#   column_to_rownames(var = "kinase") %>%
#   t() %>%
#   as.data.frame() %>%
#   rownames_to_column(var = "sample") %>%
#   as_tibble()

write_tsv(kin_act, "./output/files/kin_act_matrix_article.txt")


# load TF activity data
tf_act <- data.table::fread("./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  as_tibble() %>%
  rename(TF=V1)

write_tsv(tf_act, "./output/files/tf_act_matrix_article.txt")


# load string network (vertices names sorted for each edge; duplicates removed)
string_net <- read_tsv("./output/files/string_network_sorted_pairs.txt.gz")


# function to sort the protein pairs alphabetically
sortPrt <- function(x, y){
  sorted <- sort(c(x,y))
  tibble(a=sorted[1],b=sorted[2])
}

# load genetic associations with kinase activities
kin_genetic <- read_tsv(file = "./output/files/ka_mutStatus_allGenes.txt.gz") %>%
  filter(kinase != gene) %>%
  mutate(p.adjust = p.adjust(p.value, "BH")) %>%
  filter(p.adjust < 0.05) %>%
  mutate(x = map2(kinase, gene, sortPrt)) %>%
  unnest(cols = x) %>%
  left_join(string_net, by = c("a", "b")) %>%
  select(-a, -b) %>%
  mutate(protein_type = "kinase") %>%
  select(protein_type, protein = kinase, mutated_gene = gene, mut_effect_size = estimate, pval = p.value, adjusted_pval = p.adjust, string_net_weight = weight)

kin_genetic2 <- read_tsv(file = "output/files/ka_mutStatus_inKinase.txt.gz") %>%
  filter(p.adjust < 0.05)

# load genetic associations with TF activities
tf_genetic <- read_tsv(file = "./output/files/tf_mutStatus_allGenes.txt.gz") %>%
  filter(tf != gene) %>%
  mutate(p.adjust = p.adjust(p.value, "BH")) %>%
  filter(p.adjust < 0.05) %>%
  mutate(x = map2(tf, gene, sortPrt)) %>%
  unnest(cols = x) %>%
  left_join(string_net, by = c("a", "b")) %>%
  select(-a, -b) %>%
  mutate(protein_type = "TF") %>%
  select(protein_type, protein = tf, mutated_gene = gene, mut_effect_size = estimate, pval = p.value, adjusted_pval = p.adjust, string_net_weight = weight)

tf_genetic2 <- read_tsv(file = "output/files/tf_mutStatus_inTF.txt.gz") %>%
  filter(p.adjust < 0.05) %>%
  mutate(protein_type = "TF", mutated_gene = tf) %>%
  select(protein_type, protein = tf, mutated_gene, mut_effect_size = estimate, pval = p.value, adjusted_pval = p.adjust)

tf_genetic <- bind_rows(tf_genetic, tf_genetic2)


# merge kinase and TF genetic associations
genetic_associations <- bind_rows(kin_genetic, tf_genetic) %>%
  arrange(protein_type, adjusted_pval)

write_tsv(genetic_associations, "./output/files/genetic_associations_article.txt")


# load associations between kinases and TFs
kin_tf_associations <- read_tsv(file = "./output/files/kinaseImputed_TF_activity_associations.txt.gz") %>%
  filter(padj < 0.05) %>%
  mutate(x = map2(kinase, tf, sortPrt)) %>%
  unnest(cols = x) %>%
  left_join(string_net, by = c("a", "b")) %>%
  select(-a, -b) %>%
  select(kinase, TF = tf, kin_effect_size = kin_beta, pval = kin_pval, adjusted_pval = padj, string_net_weight = weight) %>%
  arrange(adjusted_pval)

write_tsv(kin_tf_associations, "./output/files/kin_tf_associations_article.txt")


# load survival analysis data with kinase activities
kinases_km_curves_logtest <- read_tsv("./output/files/kinases_survival_curves_logrank_test.txt") %>%
  filter(padj < 0.2) %>%
  rename(logrank_test_pval = pval, logrank_test_padj = padj)

kinases_cox_regression <- read_tsv("./output/files/kinases_cox_regression.txt") %>%
  select(kinase, cancer, cox_act_pval = act_pval, cox_act_adj_pval = padj1, cox_act_effect_size = act_beta, cox_act_hr = act_hr, cox_hrCI95_lower = hrCI95_lower, cox_hrCI95_upper = hrCI95_upper)

kinases_survival <- kinases_km_curves_logtest %>%
  inner_join(kinases_cox_regression, by = c("kinase", "cancer")) %>%
  mutate(protein_type = "kinase") %>%
  select(protein_type, protein = kinase, tissue_type = cancer, everything())


# load survival analysis data with TF activities
tfs_km_curves_logtest <- read_tsv("./output/files/TFs_survival_curves_logrank_test.txt") %>%
  filter(padj < 0.05) %>%
  rename(logrank_test_pval = pval, logrank_test_padj = padj)

tfs_cox_regression <- read_tsv("./output/files/tfs_cox_regression.txt") %>%
  select(tf, cancer, cox_act_pval = act_pval, cox_act_adj_pval = padj1, cox_act_effect_size = act_beta, cox_act_hr = act_hr, cox_hrCI95_lower = hrCI95_lower, cox_hrCI95_upper = hrCI95_upper)

tfs_survival <- tfs_km_curves_logtest %>%
  inner_join(tfs_cox_regression, by = c("tf", "cancer")) %>%
  mutate(protein_type = "TF") %>%
  select(protein_type, protein = tf, tissue_type = cancer, everything())

# merge kinase and TF survival analysis
survival_analysis <- bind_rows(kinases_survival, tfs_survival) %>%
  arrange(protein_type, logrank_test_pval)

write_tsv(survival_analysis, "./output/files/survival_analysis_article.txt")
