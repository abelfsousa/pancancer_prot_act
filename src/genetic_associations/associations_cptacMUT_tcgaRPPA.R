# perform genetic associations using CPTAC mutations and TCGA RPPA data


library(tidyverse)
source("./src/utils/KA_mut.R")
source("./src/utils/getSamples.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load Danish's RPPA phosphorylation data
load("./data/protein/tcga/rppa/pan_can_tcga_rppa_mat_updated.Rdata")
pan_can_tcga_rppa_mat <- pan_can_tcga_rppa_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "feature") %>%
  as_tibble() %>%
  select_at(.vars = c(1, which(str_sub(colnames(.),14,15) == "01"))) %>%
  pivot_longer(-feature, names_to = "sample", values_to = "expr") %>%
  mutate(sample = str_replace_all(str_sub(sample, 1, 12), "-", "."))

load("./data/protein/tcga/rppa/matched_prot_phosphoprot_list_updated.Rdata")
matched_prot_phosphoprot_list <- matched_prot_phosphoprot_list %>%
  as_tibble()

rppa_pho <- pan_can_tcga_rppa_mat %>%
  inner_join(matched_prot_phosphoprot_list, by = c("feature" = "fil_rppa_phospho_ab")) %>%
  select(prot_gene = fil_gene_prot_list, prot = fil_rppa_prot_ab, psite = feature, sample, psite_expr = expr) %>%
  mutate_if(is.factor, as.character)


# load samples metadata
metadata <- read_tsv(file = "./output/files/all_samples_annotation.txt")
metadata <- getSamples(metadata, c("mutation", "protein"))


# load mutation matrix
mut <- read_tsv(file = "./output/files/mutations_matrix.txt.gz")


# prepare batch covariate
batch <- metadata %>%
  mutate(batch = str_c(batch, tissue, sep = "_")) %>%
  select(sample, batch)


# -- perform the associations


# filter mutation matrix
mutMat <- mut %>%
  filter(n >= 20) %>%
  select(-n) %>%
  pivot_longer(cols = -c(gene), names_to = "sample", values_to = "mutated")

# prepare data to model
dataMod <- rppa_pho %>%
  inner_join(mutMat, by = c("sample")) %>%
  select(prot_gene, prot, psite, mut_gene = gene, sample, psite_expr, mutated) %>%
  group_by(prot_gene, prot, psite, mut_gene) %>%
  nest() %>%
  ungroup()

# perform the linear modelling
cores = detectCores()-1
cl <- makeCluster(cores)
clusterEvalQ(cl, library(tidyverse))

models <- parLapply(
  cl = cl,
  X = dataMod$data,
  fun = linear_model,
  covs = batch,
  mutN = 5,
  NonMutN = 0,
  kaV = "psite_expr",
  mutV = "mutated")

stopCluster(cl)

dataMod <- dataMod %>%
  mutate(model = models) %>%
  select(-data) %>%
  unnest() %>%
  filter(!is.na(estimate)) %>%
  mutate(p.adjust = p.adjust(p.value, "BH")) %>%
  rename(rppa_beta = estimate, rppa_pvalue = p.value, rppa_padj = p.adjust)


# load kinase associations from CPTAC
assoc_cptac <- read_tsv(file = "./output/files/ka_mutStatus_allGenes.txt.gz") %>%
  rename(mut_gene = gene, prot_gene = kinase) %>%
  rename(cptac_beta = estimate, cptac_pvalue = p.value, cptac_padj = p.adjust)


# overlap between significant associations
assoc_cptac %>%
  filter(cptac_padj < 0.05) %>%
  inner_join(
    dataMod %>%
      filter(rppa_padj < 0.05),
    by = c("prot_gene", "mut_gene"))


# compare associations by correlating P-values and effect sizes (betas)
assoc_cptac_rppa2 <- assoc_cptac %>%
  inner_join(dataMod, by = c("prot_gene", "mut_gene"))

assoc_cptac_rppa_pval_plot <- assoc_cptac_rppa2 %>%
  ggplot(mapping = aes(x = -log10(cptac_pvalue), y = -log10(rppa_pvalue))) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "CPTAC P-value (-log10)", y = "RPPA P-value (-log10)", title = "Genetic associations (CPTAC vs RPPA)")
ggsave("associations_cptac_rppa_pvalue_corr.png", plot = assoc_cptac_rppa_pval_plot, path = "./output/plots/genetic_associations/", width = 3, height = 3)

assoc_cptac_rppa_beta_plot <- assoc_cptac_rppa2 %>%
  ggplot(mapping = aes(x = cptac_beta, y = rppa_beta)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "CPTAC beta (effect size)", y = "RPPA beta (effect size)", title = "Genetic associations (CPTAC vs RPPA)")
ggsave("associations_cptac_rppa_beta_corr.png", plot = assoc_cptac_rppa_beta_plot, path = "./output/plots/genetic_associations/", width = 3, height = 3)


# example of association between FLT3 and AKT1/Akt_pS473 using CPTAC data

# load mutation matrix
samp_mut = 20
mut <- read_tsv(file = "./output/files/mutations_matrix.txt.gz") %>%
  filter(n >= samp_mut) %>%
  select(-n) %>%
  pivot_longer(cols = -c(gene), names_to = "sample", values_to = "mutated")

# load kinase-activity inference data
k_subN = 3
k_sampN = 10
ka <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  group_by(kinase) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n > k_sampN) %>%
  select(-n)

example <- assoc_cptac %>%
  filter(prot_gene == "AKT1", mut_gene == "FLT3") %>%
  mutate(pair = str_c(mut_gene, prot_gene, sep = "--")) %>%
  inner_join(ka, by = c("prot_gene" = "kinase")) %>%
  inner_join(mut, by = c("mut_gene" = "gene", "sample")) %>%
  mutate(mutated = fct_relevel(as.character(mutated), "0", "1"))

helper <- function(vect, lims){
  x <- vect
  x <- x[x >= min(lims) & x <= max(lims)]
  
  return(x)
}

example2 <- example %>%
  group_by(pair, mutated) %>%
  summarise(log10P = list(log10P)) %>%
  ungroup() %>%
  mutate(limits = map(.x = log10P, .f = ~ c(quantile(.x, 0.25)-(1.5*IQR(.x)), quantile(.x, 0.75)+(1.5*IQR(.x))))) %>%
  mutate(log10P = map2(.x = log10P, .y = limits, .f = helper)) %>%
  select(-limits) %>%
  unnest()

helper <- function(vect, value){
  x <- vect
  #x[x == min(x)] <- min(x)-value
  x[x == max(x)] <- max(x)+value
  
  return(x)
}

dummy <- example2 %>%
  group_by(pair, mutated) %>%
  summarise(log10P = list(range(log10P))) %>%
  ungroup() %>%
  mutate(increment = rep(0.1,2)) %>%
  mutate(log10P = map2(.x = log10P, .y = increment, .f = helper)) %>%
  unnest()

p_values <- example %>%
  select(pair, mutated, log10P) %>%
  group_by(pair) %>%
  nest() %>%
  ungroup() %>%
  mutate(p_value = map_dbl(.x = data, .f = ~ t.test(log10P ~ mutated, data = .x)$p.value)) %>%
  select(-data) %>%
  inner_join(
    dummy %>% select(-increment) %>% group_by(pair) %>% summarise(log10P = max(log10P)) %>% ungroup(),
    by = "pair") %>%
  rename(y.position = log10P) %>%
  mutate(group1 = "0", group2 = "1") %>%
  mutate(p_value = signif(p_value,2))

N <- example %>%
  #filter(pair %in% c("STK11--PRKACA", "TP53--CDK1")) %>%
  group_by(pair, mutated) %>%
  tally() %>%
  ungroup()

plot1 <- example2 %>%
  #filter(pair %in% c("STK11--PRKACA", "TP53--CDK1")) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = pair, y = log10P, fill = mutated), notch = FALSE, outlier.shape = NA, alpha = 0.8, size = 1, fatten = 1.5) +
  geom_point(mapping = aes(x = pair, y = log10P, fill = mutated), position = position_jitterdodge(jitter.width = 0.2), size = 1.5, show.legend = F, pch = 21, color = "white", alpha = 0.5) +
  #geom_violin(mapping = aes(x = pair, y = log10P, fill = mutated), alpha = 0.8, size = 1) +
  #geom_jitter(mapping = aes(x = pair, y = log10P, alpha = mutated, group = mutated), width = 0.1, size = 0.5, show.legend = F, color = "black") +
  #geom_blank(data = dummy %>% filter(pair %in% c("STK11--PRKACA", "TP53--CDK1")), mapping = aes(x = mutated, y = log10P)) +
  #facet_wrap(~ pair, scales = "fixed", nrow = 1) +
  #stat_compare_means(mapping = aes(x = mutated, y = log10P), method = "t.test", size = 3.5, label = "p.format") +
  #stat_pvalue_manual(data = p_values %>% filter(pair %in% c("STK11--PRKACA", "TP53--CDK1")), label = "p = {p_value}", xmax = NULL, size = 3.5) +
  geom_text(data = p_values, mapping = aes(x = pair, y = Inf, label = paste0("p = ", p_value)), vjust = 1, size = 4) +
  geom_text(data = N, mapping = aes(x = pair, y = -3, label = n, group = mutated), size = 4, position = position_dodge(width = 0.8)) +
  #scale_alpha_manual(values = c(0.5, 0.5)) +
  #scale_fill_manual(values = c("#1b9e77", "#d7191c")) +
  #scale_fill_manual(values = c("#bdbdbd", "#67a9cf"), labels = c("0" = "Non-mutated", "1" = "Mutated")) +
  scale_fill_manual(values = c("#d8b365", "#5ab4ac"), labels = c("0" = "Non-mutated", "1" = "Mutated")) +
  scale_x_discrete(labels = c("FLT3--AKT1" = "FLT3 - AKT1")) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(colour = "black", size = 18),
    axis.text = element_text(colour = "black", size = 16),
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", size = 16),
    strip.text = element_text(colour = "black", size = 12),
    strip.background = element_blank(),
    legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(size=0.5))) +
  labs(x = "Association", y = "Kinase activity score")
ggsave(filename = "kinase_association_example_cptac_FLT3_AKT1.png", plot = plot1, path = "./output/plots/genetic_associations/", width = 5, height = 5)
ggsave(filename = "kinase_association_example_cptac_FLT3_AKT1.pdf", plot = plot1, path = "./output/plots/genetic_associations/", width = 5, height = 5)


# example of association between FLT3 and AKT1/Akt_pS473 using RPPA data
example <- dataMod %>%
  filter(prot_gene == "AKT1", mut_gene == "FLT3", prot == "Akt", psite == "Akt_pS473") %>%
  mutate(pair = str_c(mut_gene, psite, sep = "--")) %>%
  inner_join(rppa_pho, by = c("psite", "prot_gene", "prot")) %>%
  inner_join(mut, by = c("mut_gene" = "gene", "sample")) %>%
  mutate(mutated = fct_relevel(as.character(mutated), "0", "1"))

helper <- function(vect, lims){
  x <- vect
  x <- x[x >= min(lims) & x <= max(lims)]
  
  return(x)
}

example2 <- example %>%
  group_by(pair, mutated) %>%
  summarise(psite_expr = list(psite_expr)) %>%
  ungroup() %>%
  mutate(limits = map(.x = psite_expr, .f = ~ c(quantile(.x, 0.25)-(1.5*IQR(.x)), quantile(.x, 0.75)+(1.5*IQR(.x))))) %>%
  mutate(psite_expr = map2(.x = psite_expr, .y = limits, .f = helper)) %>%
  select(-limits) %>%
  unnest()

helper <- function(vect, value){
  x <- vect
  #x[x == min(x)] <- min(x)-value
  x[x == max(x)] <- max(x)+value
  
  return(x)
}

dummy <- example2 %>%
  group_by(pair, mutated) %>%
  summarise(psite_expr = list(range(psite_expr))) %>%
  ungroup() %>%
  mutate(increment = rep(0.1,2)) %>%
  mutate(psite_expr = map2(.x = psite_expr, .y = increment, .f = helper)) %>%
  unnest()

p_values <- example %>%
  select(pair, mutated, psite_expr) %>%
  group_by(pair) %>%
  nest() %>%
  ungroup() %>%
  mutate(p_value = map_dbl(.x = data, .f = ~ t.test(psite_expr ~ mutated, data = .x)$p.value)) %>%
  select(-data) %>%
  inner_join(
    dummy %>% select(-increment) %>% group_by(pair) %>% summarise(psite_expr = max(psite_expr)) %>% ungroup(),
    by = "pair") %>%
  rename(y.position = psite_expr) %>%
  mutate(group1 = "0", group2 = "1") %>%
  mutate(p_value = signif(p_value,2))

N <- example %>%
  #filter(pair %in% c("STK11--PRKACA", "TP53--CDK1")) %>%
  group_by(pair, mutated) %>%
  tally() %>%
  ungroup()

plot1 <- example2 %>%
  #filter(pair %in% c("STK11--PRKACA", "TP53--CDK1")) %>%
  ggplot() +
  geom_boxplot(mapping = aes(x = pair, y = psite_expr, fill = mutated), notch = FALSE, outlier.shape = NA, alpha = 0.8, size = 1, fatten = 1.5) +
  geom_point(mapping = aes(x = pair, y = psite_expr, fill = mutated), position = position_jitterdodge(jitter.width = 0.2), size = 1.5, show.legend = F, pch = 21, color = "white", alpha = 0.5) +
  #geom_violin(mapping = aes(x = pair, y = psite_expr, fill = mutated), alpha = 0.8, size = 1) +
  #geom_jitter(mapping = aes(x = pair, y = psite_expr, alpha = mutated, group = mutated), width = 0.1, size = 0.5, show.legend = F, color = "black") +
  #geom_blank(data = dummy %>% filter(pair %in% c("STK11--PRKACA", "TP53--CDK1")), mapping = aes(x = mutated, y = psite_expr)) +
  #facet_wrap(~ pair, scales = "fixed", nrow = 1) +
  #stat_compare_means(mapping = aes(x = mutated, y = psite_expr), method = "t.test", size = 3.5, label = "p.format") +
  #stat_pvalue_manual(data = p_values %>% filter(pair %in% c("STK11--PRKACA", "TP53--CDK1")), label = "p = {p_value}", xmax = NULL, size = 3.5) +
  geom_text(data = p_values, mapping = aes(x = pair, y = Inf, label = paste0("p = ", p_value)), vjust = 1, size = 4) +
  geom_text(data = N, mapping = aes(x = pair, y = -2, label = n, group = mutated), size = 4, position = position_dodge(width = 0.8)) +
  #scale_alpha_manual(values = c(0.5, 0.5)) +
  #scale_fill_manual(values = c("#1b9e77", "#d7191c")) +
  #scale_fill_manual(values = c("#bdbdbd", "#67a9cf"), labels = c("0" = "Non-mutated", "1" = "Mutated")) +
  scale_fill_manual(values = c("#d8b365", "#5ab4ac"), labels = c("0" = "Non-mutated", "1" = "Mutated")) +
  scale_x_discrete(labels = c("FLT3--Akt_pS473" = "FLT3 - Akt_pS473")) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(colour = "black", size = 18),
    axis.text = element_text(colour = "black", size = 16),
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", size = 16),
    strip.text = element_text(colour = "black", size = 12),
    strip.background = element_blank(),
    legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(size=0.5))) +
  labs(x = "Association", y = "Kinase activity score")
ggsave(filename = "kinase_association_example_rppa_FLT3_Akt_pS473.png", plot = plot1, path = "./output/plots/genetic_associations/", width = 5, height = 5)
ggsave(filename = "kinase_association_example_rppa_FLT3_Akt_pS473.pdf", plot = plot1, path = "./output/plots/genetic_associations/", width = 5, height = 5)

