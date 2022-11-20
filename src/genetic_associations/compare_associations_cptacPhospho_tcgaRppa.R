#' # compare genetic associations obtained with TCGA RPPA data (from Danish) and CPTAC phosphorylation data


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T)


#' load R packages and external R functions
library(tidyverse)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load kinase associations from CPTAC
assoc_cptac <- read_tsv(file = "./output/files/ka_mutStatus_allGenes.txt.gz") %>%
  rename(mutation = gene) %>%
  rename(cptac_beta = estimate, cptac_pvalue = p.value, cptac_padj = p.adjust)


#' load Danish's phosphosite associations from TCGA (RPPA)
load("./data/protein/tcga/rppa/matched_prot_phosphoprot_list_updated.Rdata")
matched_prot_phosphoprot_list <- matched_prot_phosphoprot_list %>%
  as_tibble() %>%
  rename_with(.fn = ~ c("gene", "protein", "psite"))

load("./data/Danish/cor_tcga_mut_rppa_sig.Rdata")
assoc_rppa_pvalues <- cor_tcga_mut_rppa_sig %>%
  as.data.frame() %>%
  rownames_to_column(var = "mut_gene") %>%
  as_tibble() %>%
  pivot_longer(-mut_gene, names_to = "psite", values_to = "rppa_pvalue") %>%
  mutate(protein = str_extract(psite, "^[^_]+")) %>%
  #mutate(protein = str_extract(psite, "^[A-Za-z0-9.]+")) %>%
  #mutate(protein = str_split_fixed(feature, "_", n = 2)[,1], psite = str_split_fixed(feature, "_", n = 2)[,2]) %>%
  #separate(col = "feature", into = c("protein", "psite"), sep = "_", extra = "merge") %>%
  inner_join(distinct(matched_prot_phosphoprot_list[, c("gene", "protein")]), by = "protein") %>%
  select(kinase = gene, psite, mutation = mut_gene, rppa_pvalue) %>%
  mutate(rppa_padj = p.adjust(rppa_pvalue, method = "BH"))

load("./data/Danish/cor_tcga_mut_rppa_beta.Rdata")
assoc_rppa_betas <- cor_tcga_mut_rppa_beta %>%
  as.data.frame() %>%
  rownames_to_column(var = "mut_gene") %>%
  as_tibble() %>%
  pivot_longer(-mut_gene, names_to = "psite", values_to = "rppa_beta") %>%
  mutate(protein = str_extract(psite, "^[^_]+")) %>%
  #mutate(protein = str_extract(psite, "^[A-Za-z0-9.]+")) %>%
  #mutate(protein = str_split_fixed(feature, "_", n = 2)[,1], psite = str_split_fixed(feature, "_", n = 2)[,2]) %>%
  #separate(col = "feature", into = c("protein", "psite"), sep = "_", extra = "merge") %>%
  inner_join(distinct(matched_prot_phosphoprot_list[, c("gene", "protein")]), by = "protein") %>%
  select(kinase = gene, psite, mutation = mut_gene, rppa_beta)

assoc_rppa <- inner_join(assoc_rppa_pvalues, assoc_rppa_betas, by = c("kinase", "psite", "mutation")) %>%
  select(kinase, psite, mutation, rppa_beta, everything())


#' load kinase activity - RPPA correlations
kinase_rppa_cor <- read_tsv(file = "./output/files/corr_kin_activity_rppa.txt")


#' overlap between significant associations
assoc_cptac %>%
  filter(cptac_padj < 0.05) %>%
  inner_join(
    assoc_rppa %>%
      filter(rppa_padj < 0.05),
    by = c("kinase", "mutation")) %>%
  nrow()

assoc_cptac %>%
  filter(cptac_padj < 0.05) %>%
  select(kinase, mutation) %>%
  intersect(
    assoc_rppa %>%
      filter(rppa_padj < 0.05) %>%
      select(kinase, mutation)) %>%
  nrow()


#' join CPTAC and RPPA associations
assoc_cptac_rppa <- assoc_cptac %>%
  inner_join(assoc_rppa, by = c("kinase", "mutation")) %>%
  select(kinase, psite, mutation, everything())


#' # -- compare associations by correlating P-values and effect sizes (betas)\

#' ## perform correlations between betas and P-values across all kinases - all associations
assoc_cptac_rppa_pval_plot <- assoc_cptac_rppa %>%
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
ggsave("associations_cptac_rppa_danish_pvalue_corr.png", plot = assoc_cptac_rppa_pval_plot, path = "./output/plots/genetic_associations/", width = 3, height = 3)

#+ fig.width = 3, fig.height = 3
assoc_cptac_rppa_pval_plot

assoc_cptac_rppa_beta_plot <- assoc_cptac_rppa %>%
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
ggsave("associations_cptac_rppa_danish_beta_corr.png", plot = assoc_cptac_rppa_beta_plot, path = "./output/plots/genetic_associations/", width = 3, height = 3)

#+ fig.width = 3, fig.height = 3
assoc_cptac_rppa_beta_plot


#' example of association between PTEN and AKT1 using CPTAC data (ask Danish for the example using RPPA data)

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
  filter(kinase == "AKT1", mutation == "PTEN") %>%
  mutate(pair = str_c(mutation, kinase, sep = "--")) %>%
  inner_join(ka, by = c("kinase")) %>%
  inner_join(mut, by = c("mutation" = "gene", "sample")) %>%
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
  #ggpubr::stat_compare_means(mapping = aes(x = pair, y = log10P, group = mutated), method = "t.test", size = 3.5, label = "p.format") +
  #stat_pvalue_manual(data = p_values %>% filter(pair %in% c("STK11--PRKACA", "TP53--CDK1")), label = "p = {p_value}", xmax = NULL, size = 3.5) +
  geom_text(data = p_values, mapping = aes(x = pair, y = Inf, label = paste0("p = ", p_value)), vjust = 1, size = 4) +
  geom_text(data = N, mapping = aes(x = pair, y = -4.2, label = n, group = mutated), size = 4, position = position_dodge(width = 0.8)) +
  #scale_alpha_manual(values = c(0.5, 0.5)) +
  #scale_fill_manual(values = c("#1b9e77", "#d7191c")) +
  #scale_fill_manual(values = c("#bdbdbd", "#67a9cf"), labels = c("0" = "Non-mutated", "1" = "Mutated")) +
  scale_fill_manual(values = c("#d8b365", "#5ab4ac"), labels = c("0" = "Non-mutated", "1" = "Mutated")) +
  scale_x_discrete(labels = c("PTEN--AKT1" = "PTEN - AKT1")) +
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
ggsave(filename = "kinase_association_example_cptac_PTEN_AKT1.png", plot = plot1, path = "./output/plots/genetic_associations/", width = 5, height = 5)
ggsave(filename = "kinase_association_example_cptac_PTEN_AKT1.pdf", plot = plot1, path = "./output/plots/genetic_associations/", width = 5, height = 5)

#+ fig.width = 5, fig.height = 5
plot1


#' ## perform correlations between betas and P-values across all kinases - select CPTAC associations that are significant
assoc_cptac_rppa_pval_plot <- assoc_cptac_rppa %>%
  filter(cptac_padj < 0.05) %>%
  ggplot(mapping = aes(x = -log10(cptac_pvalue), y = -log10(rppa_pvalue))) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "CPTAC P-value (-log10)", y = "RPPA P-value (-log10)", title = "Genetic associations (CPTAC vs RPPA)")

#+ fig.width = 5, fig.height = 5
assoc_cptac_rppa_pval_plot

assoc_cptac_rppa_beta_plot <- assoc_cptac_rppa %>%
  filter(cptac_padj < 0.05) %>%
  ggplot(mapping = aes(x = cptac_beta, y = rppa_beta)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "CPTAC beta (effect size)", y = "RPPA beta (effect size)", title = "Genetic associations (CPTAC vs RPPA)")

#+ fig.width = 5, fig.height = 5
assoc_cptac_rppa_beta_plot


#' ## perform correlations between betas and P-values across all kinases - select RPPA associations that are significant
assoc_cptac_rppa_pval_plot <- assoc_cptac_rppa %>%
  filter(rppa_padj < 0.05) %>%
  ggplot(mapping = aes(x = -log10(cptac_pvalue), y = -log10(rppa_pvalue))) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "CPTAC P-value (-log10)", y = "RPPA P-value (-log10)", title = "Genetic associations (CPTAC vs RPPA)")

#+ fig.width = 5, fig.height = 5
assoc_cptac_rppa_pval_plot

assoc_cptac_rppa_beta_plot <- assoc_cptac_rppa %>%
  filter(rppa_padj < 0.05) %>%
  ggplot(mapping = aes(x = cptac_beta, y = rppa_beta)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "CPTAC beta (effect size)", y = "RPPA beta (effect size)", title = "Genetic associations (CPTAC vs RPPA)")
ggsave(filename = "cptac_rppa_beta_corr_rppa_signf.png", plot = assoc_cptac_rppa_beta_plot, path = "./output/plots/genetic_associations/", width = 4, height = 4)
ggsave(filename = "cptac_rppa_beta_corr_rppa_signf.pdf", plot = assoc_cptac_rppa_beta_plot, path = "./output/plots/genetic_associations/", width = 4, height = 4)

#+ fig.width = 4, fig.height = 4
assoc_cptac_rppa_beta_plot


#' ## perform correlations between betas and P-values across all kinases - select CPTAC and RPPA associations that are significant
assoc_cptac_rppa %>%
  filter(cptac_padj < 0.05 & rppa_padj < 0.05) %>%
  nrow()


#' ## perform correlations between betas and P-values by kinase
assoc_cptac_rppa_pval_plot <- assoc_cptac_rppa %>%
  ggplot(mapping = aes(x = -log10(cptac_pvalue), y = -log10(rppa_pvalue))) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(kinase), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "CPTAC P-value (-log10)", y = "RPPA P-value (-log10)", title = "Genetic associations (CPTAC vs RPPA)")

#+ fig.width = 8, fig.height = 8
assoc_cptac_rppa_pval_plot

assoc_cptac_rppa_beta_plot <- assoc_cptac_rppa %>%
  ggplot(mapping = aes(x = cptac_beta, y = rppa_beta)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(kinase), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "CPTAC beta (effect size)", y = "RPPA beta (effect size)", title = "Genetic associations (CPTAC vs RPPA)")

#+ fig.width = 8, fig.height = 8
assoc_cptac_rppa_beta_plot


#' ## perform correlations between betas and P-values by kinase and RPPA psite
assoc_cptac_rppa_pval_plot <- assoc_cptac_rppa %>%
  unite(kinase, psite, col = "kin_psite", sep = " ~ ") %>%
  ggplot(mapping = aes(x = -log10(cptac_pvalue), y = -log10(rppa_pvalue))) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(kin_psite), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "CPTAC P-value (-log10)", y = "RPPA P-value (-log10)", title = "Genetic associations (CPTAC vs RPPA)")

#+ fig.width = 8, fig.height = 8
assoc_cptac_rppa_pval_plot

assoc_cptac_rppa_beta_plot <- assoc_cptac_rppa %>%
  unite(kinase, psite, col = "kin_psite", sep = " ~ ") %>%
  ggplot(mapping = aes(x = cptac_beta, y = rppa_beta)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(kin_psite), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "CPTAC beta (effect size)", y = "RPPA beta (effect size)", title = "Genetic associations (CPTAC vs RPPA)")

#+ fig.width = 8, fig.height = 8
assoc_cptac_rppa_beta_plot


#' ## perform correlations between betas and P-values by kinase and RPPA psite - select CPTAC associations that are significant
assoc_cptac_rppa_pval_plot <- assoc_cptac_rppa %>%
  filter(cptac_padj < 0.05) %>%
  unite(kinase, psite, col = "kin_psite", sep = " ~ ") %>%
  filter(kin_psite == "RAF1 ~ C.Raf_pS338") %>%
  ggplot(mapping = aes(x = -log10(cptac_pvalue), y = -log10(rppa_pvalue))) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(kin_psite), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "CPTAC P-value (-log10)", y = "RPPA P-value (-log10)", title = "Genetic associations (CPTAC vs RPPA)")

#+ fig.width = 5, fig.height = 5
assoc_cptac_rppa_pval_plot

assoc_cptac_rppa_beta_plot <- assoc_cptac_rppa %>%
  filter(cptac_padj < 0.05) %>%
  unite(kinase, psite, col = "kin_psite", sep = " ~ ") %>%
  filter(kin_psite == "RAF1 ~ C.Raf_pS338") %>%
  ggplot(mapping = aes(x = cptac_beta, y = rppa_beta)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(kin_psite), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "CPTAC beta (effect size)", y = "RPPA beta (effect size)", title = "Genetic associations (CPTAC vs RPPA)")

#+ fig.width = 5, fig.height = 5
assoc_cptac_rppa_beta_plot


#' ## perform correlations between betas and P-values by kinase and RPPA psite - select RPPA associations that are significant
assoc_cptac_rppa_pval_plot <- assoc_cptac_rppa %>%
  filter(rppa_padj < 0.05) %>%
  unite(kinase, psite, col = "kin_psite", sep = " ~ ") %>%
  filter(kin_psite == "SRC ~ Src_pY527") %>%
  ggplot(mapping = aes(x = -log10(cptac_pvalue), y = -log10(rppa_pvalue))) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(kin_psite), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "CPTAC P-value (-log10)", y = "RPPA P-value (-log10)", title = "Genetic associations (CPTAC vs RPPA)")

#+ fig.width = 5, fig.height = 5
assoc_cptac_rppa_pval_plot

assoc_cptac_rppa_beta_plot <- assoc_cptac_rppa %>%
  filter(rppa_padj < 0.05) %>%
  unite(kinase, psite, col = "kin_psite", sep = " ~ ") %>%
  filter(kin_psite == "SRC ~ Src_pY527") %>%
  ggplot(mapping = aes(x = cptac_beta, y = rppa_beta)) +
  geom_point() +
  geom_smooth(method = "lm", se = F, size = 0.5) +
  ggpubr::stat_cor(label.y.npc = 1, size = 2.5) +
  facet_wrap(facets = vars(kin_psite), scales = "free") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "CPTAC beta (effect size)", y = "RPPA beta (effect size)", title = "Genetic associations (CPTAC vs RPPA)")

#+ fig.width = 5, fig.height = 5
assoc_cptac_rppa_beta_plot


#' ## perform correlations between betas and P-values by kinase - correlate with RPPA vs kinase activity correlations
kin_act_rppa_corr <- kinase_rppa_cor %>%
  filter(class1 == "Kinase-RPPA", class2 == "Substrate-based", class3 == "Protein reg-out") %>%
  select(kinase = kinaseB, psite, kin_act_rppa_corr = pearson_r)

kin_act_rppa_beta_cor_plot <- assoc_cptac_rppa %>%
  group_by(kinase, psite) %>%
  summarise(beta_corr = cor(cptac_beta, rppa_beta)) %>%
  ungroup() %>%
  inner_join(kin_act_rppa_corr, by = c("kinase", "psite")) %>%
  ggplot(mapping = aes(x = kin_act_rppa_corr, y = beta_corr)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor() +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "Kinase activity vs RPPA correlation", y = "Beta correlation of genetic associations")

#+ fig.width = 5, fig.height = 5
kin_act_rppa_beta_cor_plot

kin_act_rppa_pval_cor_plot <- assoc_cptac_rppa %>%
  group_by(kinase, psite) %>%
  summarise(pval_corr = cor(-log10(cptac_pvalue), -log10(rppa_pvalue))) %>%
  ungroup() %>%
  inner_join(kin_act_rppa_corr, by = c("kinase", "psite")) %>%
  ggplot(mapping = aes(x = kin_act_rppa_corr, y = pval_corr)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor() +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "Kinase activity vs RPPA correlation", y = "FDR correlation of genetic associations")

#+ fig.width = 5, fig.height = 5
kin_act_rppa_pval_cor_plot

kin_act_rppa_pval_cor_plot <- assoc_cptac_rppa %>%
  group_by(kinase, psite) %>%
  summarise(beta_corr = cor(cptac_beta, rppa_beta), pval_corr = cor(-log10(cptac_pvalue), -log10(rppa_pvalue))) %>%
  ungroup() %>%
  inner_join(kin_act_rppa_corr, by = c("kinase", "psite")) %>%
  pivot_longer(cols = c("beta_corr", "pval_corr"), names_to = "statistic", values_to = "value") %>%
  ggplot(mapping = aes(x = kin_act_rppa_corr, y = value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(
    facets = vars(statistic),
    scales = "free",
    labeller = labeller(statistic = c("beta_corr" = "Beta", "pval_corr" = "FDR"))) +
  ggpubr::stat_cor() +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    plot.title = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10)) +
  labs(x = "Kinase activity vs RPPA correlation", y = "Beta/FDR correlation of genetic associations")

#+ fig.width = 10, fig.height = 5
kin_act_rppa_pval_cor_plot
