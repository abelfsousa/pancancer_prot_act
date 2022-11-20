library(ggpubr)
library(tidyverse)
library(ggbeeswarm)
library(RColorBrewer)
library(viridis)


# -- BRAF V600E mutation


## CPTAC data


# load CPTAC mutation data
cptac_mutations <- read_tsv(file = "./output/files/mutations_protpos.txt.gz")


# load CPTAC tumour purity scores
cptac_purity <- read_tsv(file = "./output/files/samples_purity.txt")
pure_samples <- cptac_purity %>%
  filter(score > 0.7)


# load CPTAC kinase activity inference data
# (quantile-normalized protein regressed-out phosphorylation data)
k_subN <- 3
cptac_kin_activity <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ_Protregout_allSamp_withZtransf.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_subN) %>%
  select(-source_type, -n) %>%
  select(sample, gene=kinase, activity=log10P)


# load CPTAC TF activities
cptac_tf_activity <- read_csv(file = "./data/dorothea/CPTAC/PKG_TF_activity_log2FC_allsamples.csv") %>%
  rename(X1=`...1`) %>%
  pivot_longer(-X1, names_to = "sample", values_to = "activity") %>%
  rename(tf=X1) %>%
  select(sample, gene = tf, activity)


# brafMut <- cptac_mutations %>%
#   filter(variant_class == "Missense_Mutation") %>%
#   select(sample, gene=gene_symbol, prot_pos, aa_wt, aa_mut) %>%
#   distinct() %>%
#   #filter(!(gene == "BRAF" & !(prot_pos == 600 & aa_wt == "V" & aa_mut == "E"))) %>%
#   mutate(sel = pmap_dbl(.l = ., .f = ~ if(..2 != "BRAF"){res=1} else if(..2 == "BRAF" & ..3 == "600" & ..4 == "V" & ..5 == "E"){res=1}else{res=0})) %>%
#   filter(sel == 1) %>%
#   select(-sel) %>%
#   select(sample, gene) %>%
#   distinct() %>%
#   mutate(mutated = 1) %>%
#   group_by(gene) %>%
#   mutate(n = n()) %>%
#   ungroup() %>%
#   pivot_wider(names_from = "sample", values_from = "mutated", values_fill = list(mutated = 0)) %>%
#   select(n, gene, everything()) %>%
#   pivot_longer(-c("gene", "n"), names_to = "sample", values_to = "mutated") %>%
#   filter(gene == "BRAF") %>%
#   rename(V600E = mutated) %>%
#   select(-n)

brafMut <- cptac_mutations %>%
  filter(gene_symbol == "BRAF", prot_pos == 600, aa_wt == "V", aa_mut == "E") %>%
  select(gene=gene_symbol, sample) %>%
  distinct() %>%
  mutate(V600E = 1)

brafNotMut <- tibble(gene = "BRAF", V600E = 0) %>%
  mutate(sample = list(setdiff(unique(cptac_mutations$sample), brafMut$sample))) %>%
  unnest(cols = "sample")

brafMut <- brafMut %>%
  bind_rows(brafNotMut)


cptac_BRAF_act_BRAFV600E_plot <- brafMut %>%
  #filter(sample %in% pure_samples$sample) %>%
  inner_join(cptac_kin_activity, by = c("gene", "sample")) %>%
  ggplot(mapping = aes(x = V600E, y = activity)) +
  geom_boxplot(mapping = aes(group = as.character(V600E), fill = as.character(V600E)), outlier.shape = NA, show.legend = F) +
  geom_jitter(mapping = aes(alpha = as.character(V600E)), show.legend = F, width = 0.1, size = 1) +
  stat_cor(label.x = 0.1, size = 4.5, label.sep = "\n") +
  theme_classic() +
  scale_x_continuous(breaks = c(0,1)) +
  scale_alpha_manual(values = c(0.1, 0.5)) +
  theme(
    plot.title = element_text(color = "black", size = 13),
    axis.title = element_text(color = "black", size = 14),
    axis.text = element_text(color = "black", size = 13)) +
  labs(x = "BRAF V600E\nmutation", y = "BRAF activity", title = "CPTAC tumours")

ggsave("cptac_BRAF_activity_BRAFV600E_mutation_plot.pdf", plot = cptac_BRAF_act_BRAFV600E_plot, path = "./output/plots/genetic_associations/", height = 5, width = 2)


brafMutKin <- brafMut %>%
  inner_join(cptac_kin_activity, by = "sample") %>%
  rename(kinase1 = gene.x, kinase2 = gene.y)

brafMutKin_model <- brafMutKin %>%
  group_by(kinase1, kinase2) %>%
  nest() %>%
  ungroup() %>%
  mutate(model = map(.x=data, .f = ~ broom::tidy(lm(formula = activity ~ V600E, data = .x)) %>% filter(term == "V600E"))) %>%
  #select(-data) %>%
  unnest(cols = model) %>%
  filter(!is.na(estimate)) %>%
  mutate(padj = p.adjust(p = p.value, method = "BH")) %>%
  arrange(padj)

volcano <- brafMutKin_model %>%
  #mutate(l1 = pmap_chr(.l = ., .f = ~ if(..9 < 0.05){str_c(..1, ..2, sep = "--")}else{NA})) %>%
  mutate(l1 = pmap_chr(.l = ., .f = ~ if(..2 %in% c("BRAF", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2")){..2}else{if(..9 < 0.05){..2}else{NA}})) %>%
  mutate(l2 = if_else(padj < 0.05, "signif", "non_signif")) %>%
  mutate(padj = -log10(padj)) %>%
  ggplot(mapping = aes(x = estimate, y = padj, color = l2, label = l1)) +
  geom_point(size = 3, alpha = 0.5) +
  ggrepel::geom_text_repel(size = 6, color = "black", box.padding = 1.5, seed = 123) +
  scale_colour_manual(values = c("#525252", "red"), guide = F) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, colour = "black", hjust = 0.5),
    axis.text = element_text(size = 18, colour = "black"),
    axis.title = element_text(size = 20, colour = "black")) +
  labs(x = expression("Effect size (" ~ beta ~ ")"), y = "Adjusted P-value (-log10)", title = "BRAF V600E associations")

ggsave(filename = "cptac_kinase_activity_BRAFV600E_mutation_volcano.pdf", plot = volcano, path = "./output/plots/genetic_associations/", width = 5, height = 5)

examples <- brafMutKin_model %>%
  filter(kinase2 %in% c("CDK7", "CDK1")) %>%
  unnest(cols = data) %>%
  mutate(V600E = as.character(V600E))

helper <- function(vect, lims){
  x <- vect
  x <- x[x >= min(lims) & x <= max(lims)]
  
  return(x)
}

examples2 <- examples %>%
  group_by(kinase2, V600E) %>%
  summarise(activity = list(activity)) %>%
  ungroup() %>%
  mutate(limits = map(.x = activity, .f = ~ c(quantile(.x, 0.25)-(1.5*IQR(.x)), quantile(.x, 0.75)+(1.5*IQR(.x))))) %>%
  mutate(activity = map2(.x = activity, .y = limits, .f = helper)) %>%
  select(-limits) %>%
  unnest(cols = c(activity))

p_values <- examples %>%
  select(kinase2, V600E, activity) %>%
  group_by(kinase2) %>%
  nest() %>%
  ungroup() %>%
  mutate(p_value = map_dbl(.x = data, .f = ~ wilcox.test(activity ~ V600E, data = .x)$p.value)) %>%
  select(-data) %>%
  mutate(p_value = signif(p_value,2))

N <- examples %>%
  group_by(kinase2, V600E) %>%
  tally() %>%
  ungroup()

cptac_CDK1_act_BRAFV600E_plot <- examples2 %>%
  ggplot(mapping = aes(x = V600E, y = activity)) +
  #geom_violin(mapping = aes(fill = V600E), show.legend = F, alpha = 0.9, size = 1) +
  geom_boxplot(mapping = aes(fill = V600E), outlier.shape = NA, show.legend = F, alpha = 0.8, size = 1, fatten = 1.5) +
  geom_jitter(mapping = aes(alpha = V600E, fill = V600E), show.legend = F, width = 0.2, size = 2.5, pch = 21, color = "white") +
  facet_wrap(~ kinase2, scales = "free") +
  #stat_compare_means(mapping = aes(label = paste0("p = ", ..p.format..)), method = "wilcox.test", size = 6, label.y.npc = 1) +
  geom_text(data = p_values, mapping = aes(x = 1.5, y = Inf, label = paste0("p = ", p_value)), vjust = 1, size = 6) +
  geom_text(data = N, mapping = aes(x = V600E, y = -Inf, label = n), size = 6, vjust = -0.2) +
  theme_classic() +
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_fill_manual(values = c("#d8b365", "#5ab4ac")) +
  #scale_x_discrete(labels = c("0" = "Non-mutated", "1" = "Mutated")) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 20),
    axis.title = element_text(color = "black", size = 20),
    axis.text = element_text(color = "black", size = 18)) +
  labs(x = "BRAF V600E", y = "Activity score")

ggsave("cptac_CDK1_activity_BRAFV600E_mutation_plot.pdf", plot = cptac_CDK1_act_BRAFV600E_plot, path = "./output/plots/genetic_associations/", height = 5, width = 5)


# extension of the previous model with the KRAS mutations G12C, G12D and G12R
krasMutG12C <- cptac_mutations %>%
  filter(gene_symbol == "KRAS", prot_pos == 12, aa_wt == "G", aa_mut == "C") %>%
  select(gene=gene_symbol, sample) %>%
  distinct() %>%
  mutate(KRAS_G12C = 1)

krasNotMutG12C <- tibble(gene = "KRAS", KRAS_G12C = 0) %>%
  mutate(sample = list(setdiff(unique(cptac_mutations$sample), krasMutG12C$sample))) %>%
  unnest(cols = "sample")

krasMutG12C <- krasMutG12C %>%
  bind_rows(krasNotMutG12C)


krasMutG12D <- cptac_mutations %>%
  filter(gene_symbol == "KRAS", prot_pos == 12, aa_wt == "G", aa_mut == "D") %>%
  select(gene=gene_symbol, sample) %>%
  distinct() %>%
  mutate(KRAS_G12D = 1)

krasNotMutG12D <- tibble(gene = "KRAS", KRAS_G12D = 0) %>%
  mutate(sample = list(setdiff(unique(cptac_mutations$sample), krasMutG12D$sample))) %>%
  unnest(cols = "sample")

krasMutG12D <- krasMutG12D %>%
  bind_rows(krasNotMutG12D)


BrafKrasMutKin <- brafMut %>%
  rename(BRAF_V600E = V600E) %>%
  select(-gene) %>%
  inner_join(krasMutG12C[, -which(colnames(krasMutG12C) == "gene")], by = "sample") %>%
  inner_join(krasMutG12D[, -which(colnames(krasMutG12D) == "gene")], by = "sample") %>%
  inner_join(cptac_kin_activity, by = "sample") %>%
  rename(kinase = gene)

BrafKrasMutKin_model <- BrafKrasMutKin %>%
  filter(kinase %in% c("BRAF", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2")) %>%
  group_by(kinase) %>%
  nest() %>%
  ungroup() %>%
  mutate(model = map(.x=data, .f = ~ broom::tidy(lm(formula = activity ~ BRAF_V600E + KRAS_G12C + KRAS_G12D, data = .x)))) %>%
  select(-data) %>%
  unnest(cols = model)

# effect size barplot
BrafKrasMutKin_model %>%
  filter(term != "(Intercept)") %>%
  ggplot(mapping = aes(x = term, y = estimate, color = kinase, fill = p.value)) +
  geom_col(position = "dodge", size = 1.5) +
  scale_fill_continuous(limits = c(0,1)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text =  element_text(size = 12)) +
  scale_x_discrete(labels = ~ str_replace(.x, "_", " ")) +
  labs(x = "Mutation", y = "Effect size", color = "Kinase", fill = "P-value")
ggsave(filename = "braf_kras_mut_kin_activity_barplot.png", path = "./output/plots/genetic_associations/", width = 6, height = 4)
ggsave(filename = "braf_kras_mut_kin_activity_barplot.pdf", path = "./output/plots/genetic_associations/", width = 6, height = 4)

# violin plot KRAS_G12C
cptac_kin_activity %>%
  filter(gene %in% c("BRAF", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2")) %>%
  inner_join(krasMutG12C[, -which(colnames(krasMutG12C) == "gene")], by = "sample") %>%
  ggplot(mapping = aes(x = as.character(KRAS_G12C), y = activity, fill = gene)) +
  geom_violin(position = position_dodge(0.8), width = 0.8, trim = FALSE, alpha = 0.6) +
  geom_boxplot(position = position_dodge(0.8), width = 0.2, outlier.shape = NA, show.legend = F, size = 0.3, color = "black") +
  scale_fill_viridis(discrete = T) +
  theme_bw() +
  labs(x = "KRAS_G12C", y = "Kinase activity", fill = "Kinase")

# box plot KRAS_G12C
cptac_kin_activity %>%
  filter(gene %in% c("BRAF", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2")) %>%
  inner_join(krasMutG12C[, -which(colnames(krasMutG12C) == "gene")], by = "sample") %>%
  ggplot(mapping = aes(x = as.character(KRAS_G12C), y = activity, fill = gene)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = T) +
  theme_bw() +
  labs(x = "KRAS_G12C", y = "Kinase activity", fill = "Kinase")


example <- cptac_kin_activity %>%
  filter(gene %in% c("BRAF", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2")) %>%
  inner_join(krasMutG12C[, -which(colnames(krasMutG12C) == "gene")], by = "sample")

helper <- function(vect, lims){
  x <- vect
  x <- x[x >= min(lims) & x <= max(lims)]
  
  return(x)
}

example2 <- example %>%
  group_by(gene, KRAS_G12C) %>%
  summarise(activity = list(activity)) %>%
  ungroup() %>%
  mutate(limits = map(.x = activity, .f = ~ c(quantile(.x, 0.25)-(1.5*IQR(.x)), quantile(.x, 0.75)+(1.5*IQR(.x))))) %>%
  mutate(activity = map2(.x = activity, .y = limits, .f = helper)) %>%
  select(-limits) %>%
  unnest(cols = c(activity))

p_values <- example %>%
  select(gene, KRAS_G12C, activity) %>%
  group_by(gene) %>%
  nest() %>%
  ungroup() %>%
  mutate(p_value = map_dbl(.x = data, .f = ~ wilcox.test(activity ~ KRAS_G12C, data = .x)$p.value)) %>%
  select(-data) %>%
  mutate(p_value = signif(p_value,2))

N <- example %>%
  group_by(gene, KRAS_G12C) %>%
  tally() %>%
  ungroup()

# violin plot KRAS_G12C after removing outliers
example2 %>%
  ggplot(mapping = aes(x = as.character(KRAS_G12C), y = activity, fill = gene)) +
  geom_violin(position = position_dodge(0.8), width = 0.8, trim = FALSE, alpha = 0.6) +
  geom_boxplot(position = position_dodge(0.8), width = 0.2, outlier.shape = NA, show.legend = F, size = 0.3, color = "black") +
  scale_fill_viridis(discrete = T) +
  theme_bw() +
  labs(x = "KRAS_G12C", y = "Kinase activity", fill = "Kinase")

# box plot KRAS_G12C after removing outliers
set.seed(123)
example2 %>%
  mutate(KRAS_G12C = as.character(KRAS_G12C)) %>%
  ggplot(mapping = aes(x = KRAS_G12C, y = activity, fill = gene)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(mapping = aes(alpha = KRAS_G12C), position = position_jitterdodge(jitter.width = 0.05), pch = 21, show.legend = F) +
  #scale_fill_viridis(discrete = T, direction = -1) +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  scale_alpha_manual(values = c(0.2, 1)) +
  theme_bw() +
  scale_x_discrete(labels = c("0" = "KRAS WT", "1" = "KRAS G12C")) +
  labs(x = "KRAS mutation status", y = "Kinase activity", fill = "Kinase")
ggsave(filename = "KRAS_G12C_kin_activity_boxplot.png", path = "./output/plots/genetic_associations/", width = 6, height = 4)
ggsave(filename = "KRAS_G12C_kin_activity_boxplot.pdf", path = "./output/plots/genetic_associations/", width = 6, height = 4)

# box plot KRAS_G12C after removing outliers (without jitter)
set.seed(123)
example2 %>%
  mutate(KRAS_G12C = as.character(KRAS_G12C)) %>%
  ggplot(mapping = aes(x = KRAS_G12C, y = activity, fill = gene)) +
  geom_boxplot(outlier.shape = NA) +
  #scale_fill_viridis(discrete = T, direction = -1) +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  scale_alpha_manual(values = c(0.5, 1)) +
  theme_bw() +
  scale_x_discrete(labels = c("0" = "KRAS WT", "1" = "KRAS G12C")) +
  labs(x = "KRAS mutation status", y = "Kinase activity", fill = "Kinase")

# box plot (horizontal) KRAS_G12C after removing outliers
set.seed(123)
example2 %>%
  mutate(KRAS_G12C = as.character(KRAS_G12C)) %>%
  mutate(KRAS_G12C = fct_relevel(KRAS_G12C, "1", "0"), gene = fct_rev(gene)) %>%
  ggplot(mapping = aes(x = KRAS_G12C, y = activity, fill = gene)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(mapping = aes(alpha = KRAS_G12C), position = position_jitterdodge(jitter.width = 0.2), pch = 21, show.legend = F) +
  #scale_fill_viridis(discrete = T, direction = -1) +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  scale_alpha_manual(values = c(1, 0.5)) +
  theme_bw() +
  coord_flip() +
  scale_x_discrete(labels = c("0" = "KRAS WT", "1" = "KRAS G12C")) +
  guides(fill = guide_legend(reverse = T)) +
  labs(x = "KRAS mutation status", y = "Kinase activity", fill = "Kinase")

# box plot (horizontal) KRAS_G12C after removing outliers (without jitter)
set.seed(123)
example2 %>%
  mutate(KRAS_G12C = as.character(KRAS_G12C)) %>%
  mutate(KRAS_G12C = fct_relevel(KRAS_G12C, "1", "0"), gene = fct_rev(gene)) %>%
  ggplot(mapping = aes(x = KRAS_G12C, y = activity, fill = gene)) +
  geom_boxplot(outlier.shape = NA) +
  #scale_fill_viridis(discrete = T, direction = -1) +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  scale_alpha_manual(values = c(1, 0.5)) +
  theme_bw() +
  coord_flip() +
  scale_x_discrete(labels = c("0" = "KRAS WT", "1" = "KRAS G12C")) +
  guides(fill = guide_legend(reverse = T)) +
  labs(x = "KRAS mutation status", y = "Kinase activity", fill = "Kinase")


# phosphatases differentially expressed (at mRNA and protein level) between BRAFV600E samples with low and high BRAF activity
brafMut %>%
  filter(V600E == 1) %>%
  inner_join(cptac_kin_activity, by = c("sample", "gene")) %>%
  mutate(sample = fct_reorder(sample, activity, .fun = function(x) x), sample = fct_rev(sample)) %>%
  ggplot(mapping = aes(x = sample, y = activity, fill = activity)) +
  geom_col(show.legend = F) +
  theme_bw() +
  coord_flip() +
  labs(x = "Sample", y = "BRAF activity")

brafMut_LowHigh <- brafMut %>%
  filter(V600E == 1) %>%
  inner_join(cptac_kin_activity, by = c("sample", "gene")) %>%
  filter(!sample %in% c("05CO041", "01CO022", "MDST8", "01CO014")) %>%
  mutate(class = if_else(activity < 0, "low-BRAF", "high-BRAF"))


# load RNA data
rna <- read_tsv(file = "./output/files/transcriptomics_log2fpkm.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "rna_log2fpkm") %>%
  filter(!is.na(rna_log2fpkm))


# load genes related with phosphatase activity (GOMF_PHOSPHATASE_ACTIVITY)
phosphatases <- read_tsv(file = "data/gene_collections/GOMF_PHOSPHATASE_ACTIVITY.txt", col_names = c("gene")) 

wt_f <- function(tab, x, y){
  pval <- wilcox.test(formula = as.formula(str_c(y, "~", x)), data = tab) %>%
    pluck("p.value")
  lfc <- mean(tab[tab[[x]] == "low-BRAF", y, drop = T]) - mean(tab[tab[[x]] == "high-BRAF", y, drop = T])
  
  res <- tibble(pval, lfc)
  
  return(res)
}

diff_phospha_rna <- rna %>%
  semi_join(phosphatases, by = "gene") %>%
  inner_join(brafMut_LowHigh[, c("sample", "class")], by = "sample") %>%
  group_by(gene) %>%
  filter(n() == 14) %>%
  nest() %>%
  ungroup() %>%
  mutate(res = map(.x = data, .f = wt_f, x = "class", y = "rna_log2fpkm")) %>%
  unnest(res) %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "BH"))


diff_phospha_rna %>%
  mutate(flag = if_else(padj < 0.15 & lfc > 0, "1", "0")) %>%
  mutate(annot = ifelse(padj < 0.15 & lfc > 0, gene, NA)) %>%
  ggplot(mapping = aes(x = lfc, y = -log10(padj))) +
  geom_point(mapping = aes(color = flag), show.legend = F, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.15), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(mapping = aes(label = annot), box.padding = 0.5, size = 4, segment.size = 0.3, max.overlaps = 100, seed = 123) +
  scale_color_manual(values = c("1" = "red", "0" = "black")) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text =  element_text(size = 12)) +
  scale_y_continuous(limits = c(NA, 1.2)) +
  labs(x = "Fold-change (log2)", y = "FDR (-log10)")
ggsave(filename = "phosphatase_diff_expr_volcano.png", path = "./output/plots/genetic_associations/", width = 6, height = 4)
ggsave(filename = "phosphatase_diff_expr_volcano.pdf", path = "./output/plots/genetic_associations/", width = 6, height = 4)


diff_phospha_rna %>%
  slice(1) %>%
  unnest(data) %>%
  ggplot(mapping = aes(x = class, y = rna_log2fpkm)) +
  geom_boxplot() +
  geom_jitter(width = 0.1)


# load protein data
protein <- read_tsv(file = "./output/files/proteomics.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "protein_log2fc") %>%
  filter(!is.na(protein_log2fc))

# protein %>%
#   semi_join(phosphatases, by = "gene") %>%
#   inner_join(brafMut_LowHigh[, c("sample", "class")], by = "sample") %>%
#   group_by(gene, class) %>%
#   summarise(n = n()) %>%
#   ungroup() %>%
#   group_by(gene) %>%
#   filter(n() == 2) %>%
#   filter(n[which(class == "high-BRAF")] >= 4 & n[which(class == "low-BRAF")] >= 4) %>%
#   ungroup()

diff_phospha_prot <- protein %>%
  semi_join(phosphatases, by = "gene") %>%
  inner_join(brafMut_LowHigh[, c("sample", "class")], by = "sample") %>%
  group_by(gene) %>%
  filter(length(unique(class)) == 2) %>%
  filter(sum(class == "low-BRAF") >= 4 & sum(class == "high-BRAF") >= 4) %>%
  nest() %>%
  ungroup() %>%
  mutate(res = map(.x = data, .f = wt_f, x = "class", y = "protein_log2fc")) %>%
  unnest(res) %>%
  arrange(pval) %>%
  mutate(padj = p.adjust(pval, method = "BH"))

diff_phospha_prot %>%
  ggplot(mapping = aes(x = lfc, y = -log10(pval))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw()

diff_phospha_prot %>%
  slice(1) %>%
  unnest(data) %>%
  ggplot(mapping = aes(x = class, y = protein_log2fc)) +
  geom_boxplot() +
  geom_jitter(width = 0.1)


diff_phospha_rna %>%
  select(gene, rna_lfc = lfc, rna_pval = pval) %>%
  inner_join(
    diff_phospha_prot %>%
      select(gene, prot_lfc = lfc, prot_pval = pval), by = "gene") %>%
  ggplot(mapping = aes(x = rna_lfc, y = prot_lfc)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor() +
  theme_bw()

diff_phospha_rna %>%
  select(gene, rna_lfc = lfc, rna_pval = pval) %>%
  inner_join(
    diff_phospha_prot %>%
      select(gene, prot_lfc = lfc, prot_pval = pval), by = "gene") %>%
  ggplot(mapping = aes(x = -log10(rna_pval), y = -log10(prot_pval))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor() +
  theme_bw()


# correlation between the gene expression of phosphatases and kinase activities
phos_kin_corrs <- rna %>%
  semi_join(phosphatases, by = "gene") %>%
  rename(phosphatase = gene) %>%
  inner_join(
    cptac_kin_activity %>%
      filter(gene %in% c("BRAF", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2")) %>%
      rename(kinase = gene),
    by = "sample") %>%
  group_by(phosphatase, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(
    cor = map(
      .x=data,
      .f = ~ broom::tidy(cor.test(.x$rna_log2fpkm, .x$activity)) %>% 
        select(pearson_corr = estimate, pearson_pval = p.value))) %>%
  #select(-data) %>%
  unnest(cols = cor)

# phos_kin_corrs %>%
#   select(-data) %>%
#   left_join(
#     diff_phospha_rna %>%
#       filter(padj < 0.15 & lfc > 0) %>%
#       select(phosphatase = gene) %>%
#       mutate(diff = "1"),
#     by = "phosphatase") %>%
#   mutate(diff = replace_na(diff, 0)) %>%
#   arrange(phosphatase, pearson_corr) %>%
#   write_tsv(file = "~/Downloads/phosphatases.txt")


phos_kin_corrs %>%
  #filter(pearson_pval < 0.05) %>%
  group_by(phosphatase) %>%
  summarise(pearson_corr_mean = mean(pearson_corr), pearson_corr_sd = sd(pearson_corr)) %>%
  ungroup() %>%
  left_join(
    diff_phospha_rna %>%
      filter(padj < 0.15 & lfc > 0) %>%
      select(phosphatase = gene) %>%
      mutate(diff = "1"),
    by = "phosphatase") %>%
  mutate(diff = replace_na(diff, "0")) %>%
  mutate(phosphatase = fct_reorder(phosphatase, pearson_corr_mean, function(x) {x})) %>%
  mutate(diff = fct_relevel(diff, "1")) %>%
  ggplot(mapping = aes(x = phosphatase, y = pearson_corr_mean)) +
  geom_col(mapping = aes(fill = diff)) +
  geom_point(size = 0.5) +
  geom_errorbar(mapping = aes(ymin = pearson_corr_mean - pearson_corr_sd, ymax = pearson_corr_mean + pearson_corr_sd), size = 0.3) +
  scale_fill_manual(values = c("1" = "red", "0" = "grey60"), labels = c("1" = "Differentially expressed high vs low BRAF activity", "0" = "Non differentially expressed")) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  labs(x = "Phosphatase", y = "Pearson's r (BRAF-MAPK signalling pathway)")


phos_kin_corrs %>%
  #filter(pearson_pval < 0.05) %>%
  group_by(phosphatase) %>%
  summarise(pearson_corr_mean = mean(pearson_corr), pearson_corr_sd = sd(pearson_corr)) %>%
  ungroup() %>%
  left_join(
    diff_phospha_rna %>%
      filter(padj < 0.15 & lfc > 0) %>%
      select(phosphatase = gene) %>%
      mutate(diff = "1"),
    by = "phosphatase") %>%
  mutate(diff = replace_na(diff, "0")) %>%
  filter(diff == 1) %>%
  arrange(pearson_corr_mean)

phos_kin_corrs %>%
  filter(phosphatase == "MTMR14") %>%
  unnest(cols = "data") %>%
  ggplot(mapping = aes(x = rna_log2fpkm, y = activity)) +
  geom_point(alpha = 0.5) +
  facet_wrap(facets = vars(kinase)) +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor() +
  theme_bw() + 
  labs(x = "MTMR14 mRNA expression", y = "Kinase activity")


phos_kin_corrs %>%
  filter(pearson_pval < 0.05) %>%
  filter(phosphatase == "PTP4A2") %>%
  unnest(cols = "data") %>%
  ggplot(mapping = aes(x = rna_log2fpkm, y = activity)) +
  geom_point(alpha = 0.5) +
  facet_wrap(facets = vars(kinase)) +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor() +
  theme_bw() + 
  labs(x = "PTP4A2 mRNA expression", y = "Kinase activity")


diff_phospha_rna %>%
  filter(gene == "MTMR14") %>%
  unnest(data) %>%
  ggplot(mapping = aes(x = class, y = rna_log2fpkm)) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_bw()



# only correlations with BRAF kinase
phos_kin_corrs <- rna %>%
  semi_join(phosphatases, by = "gene") %>%
  rename(phosphatase = gene) %>%
  inner_join(
    cptac_kin_activity %>%
      filter(gene %in% c("BRAF")) %>%
      rename(kinase = gene),
    by = "sample") %>%
  group_by(phosphatase, kinase) %>%
  nest() %>%
  ungroup() %>%
  filter(map_lgl(data, ~ nrow(.x) >= 10)) %>%
  mutate(
    cor = map(
      .x=data,
      .f = ~ broom::tidy(cor.test(.x$rna_log2fpkm, .x$activity)) %>% 
        select(pearson_corr = estimate, pearson_pval = p.value))) %>%
  #select(-data) %>%
  unnest(cols = cor)


phos_kin_corrs %>%
  left_join(
    diff_phospha_rna %>%
      filter(padj < 0.15 & lfc > 0) %>%
      select(phosphatase = gene) %>%
      mutate(diff = "1"),
    by = "phosphatase") %>%
  mutate(diff = replace_na(diff, "0")) %>%
  mutate(phosphatase = fct_reorder(phosphatase, pearson_corr, function(x) {x})) %>%
  mutate(diff = fct_relevel(diff, "1")) %>%
  ggplot(mapping = aes(x = phosphatase, y = pearson_corr)) +
  geom_col(mapping = aes(fill = diff)) +
  #geom_point(size = 0.5) +
  scale_fill_manual(values = c("1" = "red", "0" = "grey60"), labels = c("1" = "Differentially expressed high vs low BRAF activity", "0" = "Non differentially expressed")) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  labs(x = "Phosphatase", y = "Pearson's r (BRAF activity)")


#
brafMutKin %>%
  filter(kinase1 == "BRAF" & kinase2 == "BRAF") %>%
  lm(formula = activity ~ V600E, data = .) %>%
  broom::tidy()

brafMutKin %>%
  filter(kinase1 == "BRAF" & kinase2 == "BRAF") %>%
  lm(formula = activity ~ V600E, data = .) %>%
  summary() %>%
  pluck("r.squared")

phosphatases_genes <- phos_kin_corrs %>%
  left_join(
    diff_phospha_rna %>%
      filter(padj < 0.15 & lfc > 0) %>%
      select(phosphatase = gene) %>%
      mutate(diff = "1"),
    by = "phosphatase") %>%
  mutate(diff = replace_na(diff, "0")) %>%
  filter(diff == "1") %>%
  arrange(pearson_corr) %>%
  pull(phosphatase) %>%
  head(2)

brafMutKin %>%
  filter(kinase1 == "BRAF" & kinase2 == "BRAF") %>%
  inner_join(
    rna %>%
      filter(gene %in% phosphatases_genes) %>%
      pivot_wider(names_from = "gene", values_from = "rna_log2fpkm"),
    by = "sample") %>%
  lm(formula = activity ~ V600E + MTMR14 + PTP4A2, data = .) %>%
  broom::tidy()

brafMutKin %>%
  filter(kinase1 == "BRAF" & kinase2 == "BRAF") %>%
  inner_join(
    rna %>%
      filter(gene %in% phosphatases_genes) %>%
      pivot_wider(names_from = "gene", values_from = "rna_log2fpkm"),
    by = "sample") %>%
  lm(formula = activity ~ V600E + MTMR14 + PTP4A2, data = .) %>%
  summary() %>%
  pluck("adj.r.squared")

mod <- brafMutKin %>%
  filter(kinase1 == "BRAF" & kinase2 == "BRAF") %>%
  select(-starts_with("kin")) %>%
  inner_join(
    rna %>%
      filter(gene %in% phosphatases_genes) %>%
      pivot_wider(names_from = "gene", values_from = "rna_log2fpkm"),
    by = "sample") %>%
  as.data.frame() %>% 
  column_to_rownames(var = "sample") %>%
  lm(formula = activity ~ V600E + MTMR14 + PTP4A2, data = .)

tibble(prediction_braf_act = mod$fitted.values, observed_braf_act = mod$model[, "activity"]) %>%
  ggplot(mapping = aes(x = observed_braf_act, y = prediction_braf_act)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  stat_cor()

# brafMutKin %>%
#   filter(kinase1 == "BRAF" & kinase2 == "BRAF") %>%
#   inner_join(
#     rna %>%
#       filter(gene %in% phosphatases_genes) %>%
#       pivot_wider(names_from = "gene", values_from = "rna_log2fpkm"),
#     by = "sample") %>%
#   lm(formula = activity ~ MTMR14 + PTP4A2 + V600E, data = .) %>%
#   summary()


phosphatases_genes <- phos_kin_corrs %>%
  left_join(
    diff_phospha_rna %>%
      filter(padj < 0.15 & lfc > 0) %>%
      select(phosphatase = gene) %>%
      mutate(diff = "1"),
    by = "phosphatase") %>%
  mutate(diff = replace_na(diff, "0")) %>%
  filter(diff == "1" & pearson_corr < 0) %>%
  pull(phosphatase)

brafMutKin %>%
  filter(kinase1 == "BRAF" & kinase2 == "BRAF") %>%
  inner_join(
    rna %>%
      filter(gene %in% phosphatases_genes) %>%
      pivot_wider(names_from = "gene", values_from = "rna_log2fpkm"),
    by = "sample") %>%
  lm(formula = as.formula(str_c("activity ~ V600E", str_c(phosphatases_genes, collapse = " + "), sep = " + ")), data = .) %>%
  broom::tidy()

brafMutKin %>%
  filter(kinase1 == "BRAF" & kinase2 == "BRAF") %>%
  inner_join(
    rna %>%
      filter(gene %in% phosphatases_genes) %>%
      pivot_wider(names_from = "gene", values_from = "rna_log2fpkm"),
    by = "sample") %>%
  lm(formula = as.formula(str_c("activity ~ V600E", str_c(phosphatases_genes, collapse = " + "), sep = " + ")), data = .) %>%
  summary() %>%
  pluck("adj.r.squared")



# 
brafMutTF <- brafMut %>%
  inner_join(cptac_tf_activity, by = "sample") %>%
  rename(kinase = gene.x, tf = gene.y)

brafMutTF_model <- brafMutTF %>%
  group_by(kinase, tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(model = map(.x=data, .f = ~ broom::tidy(lm(formula = activity ~ V600E, data = .x)) %>% filter(term == "V600E"))) %>%
  #select(-data) %>%
  unnest(cols = model) %>%
  filter(!is.na(estimate)) %>%
  mutate(padj = p.adjust(p = p.value, method = "BH")) %>%
  arrange(padj)

volcano <- brafMutTF_model %>%
  mutate(l1 = pmap_chr(.l = ., .f = ~ if(..9 < 0.05){str_c(..1, ..2, sep = "--")}else{NA})) %>%
  mutate(l2 = if_else(padj < 0.05, "signif", "non_signif")) %>%
  mutate(padj = -log10(padj)) %>%
  ggplot(mapping = aes(x = estimate, y = padj, color = l2, label = l1)) +
  geom_point(size = 1, alpha = 0.6) +
  ggrepel::geom_text_repel(size = 3, segment.size = 0.3, color = "black", box.padding = 1, seed = 123) +
  scale_colour_manual(values = c("grey", "red"), guide = F) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 16, colour = "black")) +
  labs(x = expression("Effect size (" ~ beta ~ ")"), y = "Adjusted P-value (-log10)", title = "TF ~ BRAF V600E\nassociations")

ggsave(filename = "cptac_tf_activity_BRAFV600E_mutation_volcano.pdf", plot = volcano, path = "./output/plots/genetic_associations/", width = 3, height = 5)



## TCGA RPPA data


# load RPPA phosphorylation data
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
  inner_join(matched_prot_phosphoprot_list[, c("fil_rppa_phospho_ab", "fil_gene_prot_list")], by = c("feature" = "fil_rppa_phospho_ab")) %>%
  select(psite = feature, protein = fil_gene_prot_list, sample, psite_expr = expr) %>%
  mutate_if(is.factor, as.character)


brafMutRPPA <- brafMut %>%
  inner_join(rppa_pho, by = "sample")

brafMutRPPA_model <- brafMutRPPA %>%
  group_by(gene, protein, psite) %>%
  nest() %>%
  ungroup() %>%
  mutate(model = map(.x=data, .f = ~broom::tidy(lm(data = .x, formula = psite_expr~V600E)) %>% filter(term == "V600E"))) %>%
  unnest(cols = model) %>%
  mutate(padj = p.adjust(p.value, "BH")) %>%
  arrange(padj)

brafMutRPPA_model %>%
  pull(data) %>%
  pluck(1) %>%
  ggplot(aes(x = as.character(V600E), y = (psite_expr))) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  labs( x = "V600E", y = "PRKCA: PKC.alpha_pS657")

brafMutRPPA_model %>% 
  pull(data) %>%
  pluck(2) %>%
  ggplot(aes(x = as.character(V600E), y = (psite_expr))) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  labs( x = "BRAF V600E", y = "EGFR: EGFR_pY1068")




## NCI60/CRC65 data


# load cancer cell lines
nci60_crc65 <- read_tsv(file = "./output/files/nci60_crc65_cell_lines.txt")

shared_cell_lines <- nci60_crc65 %>%
  group_by(batch) %>%
  summarise(cell_line = list(cell_line)) %>%
  pull(cell_line) %>%
  reduce(intersect)


# load NCI60/CRC65 mutations
cells_mutations <- read_tsv(file = "./output/files/nci60_crc65_mutations.txt.gz")


# load kinase activity inference data from NCI60/CRC65 dataset
# select quantifications with more than 3 substrates
cells_kin_activity <- read_tsv(file = "./output/files/nci60_crc65_kinase_activities.txt.gz") %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= 3) %>%
  select(-source_type, -n) %>%
  filter(!(batch == "CRC65" & sample %in% shared_cell_lines)) %>%
  select(sample, gene = kinase, activity = log10P)


# load TF activity inference data from NCI60/CRC65 dataset
tf_activity_crc65 <- read_csv(file = "./data/dorothea/nci60_crc65/TF_activity_log2FC_allsamples_CRC65.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "activity") %>%
  rename(gene=X1) %>%
  filter(!sample %in% shared_cell_lines) %>%
  select(sample, gene, activity)

tf_activity_nci60 <- read_csv(file = "./data/dorothea/nci60_crc65/TF_activity_log2FC_allsamples_NCI60.csv") %>%
  pivot_longer(-X1, names_to = "sample", values_to = "activity") %>%
  rename(gene=X1) %>%
  select(sample, gene, activity)

cells_tf_activity <- bind_rows(tf_activity_crc65, tf_activity_nci60)



brafMut_cells <- cells_mutations %>%
  filter(gene == "BRAF", prot_pos_i == 600, aa_wt == "V", aa_mut == "E") %>%
  select(gene, sample) %>%
  distinct() %>%
  mutate(V600E = 1)

brafNotMut_cells <- tibble(gene = "BRAF", V600E = 0) %>%
  mutate(sample = list(setdiff(unique(cells_mutations$sample), brafMut_cells$sample))) %>%
  unnest(cols = "sample")

brafMut_cells <- brafMut_cells %>%
  bind_rows(brafNotMut_cells)


cells_BRAF_act_BRAFV600E_plot <- brafMut_cells %>%
  inner_join(cells_kin_activity, by = c("gene", "sample")) %>%
  ggplot(mapping = aes(x = V600E, y = activity)) +
  geom_boxplot(mapping = aes(group = as.character(V600E), fill = as.character(V600E)), outlier.shape = NA, show.legend = F) +
  geom_jitter(mapping = aes(alpha = as.character(V600E)), show.legend = F, width = 0.1, size = 1) +
  stat_cor(label.x = 0.1, size = 4.5, label.sep = "\n") +
  theme_classic() +
  scale_x_continuous(breaks = c(0,1)) +
  scale_alpha_manual(values = c(0.1, 0.5)) +
  theme(
    plot.title = element_text(color = "black", size = 13),
    axis.title = element_text(color = "black", size = 14),
    axis.text = element_text(color = "black", size = 13)) +
  labs(x = "BRAF V600E\nmutation", y = "BRAF activity", title = "NCI60/CRC65 panels")

ggsave("nci60_crc65_BRAF_activity_BRAFV600E_mutation_plot.pdf", plot = cells_BRAF_act_BRAFV600E_plot, path = "./output/plots/genetic_associations/", height = 5, width = 2)


brafMutKin_cells <- brafMut_cells %>%
  inner_join(cells_kin_activity, by = "sample") %>%
  rename(kinase1 = gene.x, kinase2 = gene.y)

brafMutKin_model_cells <- brafMutKin_cells %>%
  group_by(kinase1, kinase2) %>%
  nest() %>%
  ungroup() %>%
  mutate(model = map(.x=data, .f = ~ broom::tidy(lm(formula = activity ~ V600E, data = .x)) %>% filter(term == "V600E"))) %>%
  #select(-data) %>%
  unnest(cols = model) %>%
  filter(!is.na(estimate)) %>%
  mutate(padj = p.adjust(p = p.value, method = "BH")) %>%
  arrange(p.value)

brafMutKin_model_cells %>%
  filter(p.value < 0.1) %>%
  mutate(kinase2 = fct_reorder(kinase2, p.value, function(x) x)) %>%
  unnest(cols = data) %>%
  ggplot(aes(x = as.character(V600E), y = activity)) +
  geom_boxplot() +
  #geom_jitter(width = 0.1) +
  facet_wrap(~ kinase2, scales = "free")



brafMutTF_cells <- brafMut_cells %>%
  inner_join(cells_tf_activity, by = "sample") %>%
  rename(kinase = gene.x, tf = gene.y)

brafMutTF_model_cells <- brafMutTF_cells %>%
  group_by(kinase, tf) %>%
  nest() %>%
  ungroup() %>%
  mutate(model = map(.x=data, .f = ~ broom::tidy(lm(formula = activity ~ V600E, data = .x)) %>% filter(term == "V600E"))) %>%
  #select(-data) %>%
  unnest(cols = model) %>%
  filter(!is.na(estimate)) %>%
  mutate(padj = p.adjust(p = p.value, method = "BH")) %>%
  arrange(padj)
