library(tidyverse)
library(KEGGREST)
library(ggrepel)
library(ggpubr)

source("./src/utils/getSamples.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


# load samples annotation
samples_annotation <- read_tsv(file = "./output/files/all_samples_annotation.txt")
samples_annotation <- getSamples(samples_annotation, c("phosphorylation", "mRNA"))


# load kinase-pathways associations
associations <- read_tsv(file = "./output/files/kinaseImputed_pathway_activity_associations.txt.gz")


# list of human genes
genes <- keggList("hsa")
genes <- tibble(kegg_gene = names(genes), gene_name = unname(genes))
genes <- genes %>%
  mutate(gene_name = str_split_fixed(gene_name, ";", 2)[,1]) %>%
  mutate(gene_name = str_split(gene_name, ", ")) %>%
  unnest()

# list of human pathways
pathways <- keggList("pathway", "hsa")
pathways <- tibble(kegg_pathway = names(pathways), pathway_name = unname(pathways))


# human genes linked from each of the KEGG pathway
links <- keggLink("hsa", "pathway")
links <- tibble(pathway = names(links), gene = unname(links))

links <- links %>%
  inner_join(pathways, by = c("pathway" = "kegg_pathway")) %>%
  inner_join(genes, by = c("gene" = "kegg_gene")) %>%
  select(pathway_name, gene = gene_name) %>%
  group_by(pathway_name) %>%
  summarise(gene = list(gene)) %>%
  ungroup()

map <- tibble(
  pathway = c("Androgen", "EGFR", "Estrogen", "Hypoxia", "JAK-STAT", "MAPK", "NFkB", "p53", "PI3K", "TGFb", "TNFa", "Trail", "VEGF", "WNT"),
  pathway_name = c(NA, "EGFR tyrosine kinase inhibitor resistance - Homo sapiens (human)", "Estrogen signaling pathway - Homo sapiens (human)", "HIF-1 signaling pathway - Homo sapiens (human)", "JAK-STAT signaling pathway - Homo sapiens (human)", "MAPK signaling pathway - Homo sapiens (human)", "NF-kappa B signaling pathway - Homo sapiens (human)", "p53 signaling pathway - Homo sapiens (human)", "PI3K-Akt signaling pathway - Homo sapiens (human)", "TGF-beta signaling pathway - Homo sapiens (human)", "TNF signaling pathway - Homo sapiens (human)", "Apoptosis - Homo sapiens (human)", "VEGF signaling pathway - Homo sapiens (human)", "Wnt signaling pathway - Homo sapiens (human)"))

links <- links %>%
  inner_join(map, by = "pathway_name")


# volcano plot of kinase-pathways associations
vplot <- associations %>%
  mutate(l = pmap_chr(.l = ., .f = ~ if(-log10(..5) > 4){str_c(..1, ..2, sep="--")}else{NA})) %>%
  left_join(links[, c("pathway", "gene")], by = "pathway") %>%
  mutate(p = map2_lgl(.x = kinase, .y = gene, .f = ~ if(is.null(.y)){NA}else{.x %in% .y})) %>%
  ggplot(mapping = aes(x = kin_beta, y = -log10(padj), label = l, color = p)) +
  geom_point(size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.5, color = "grey60") + 
  geom_text_repel(size = 2, segment.size = 0.2, color = "black") +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 12),
    legend.key.size = unit(0.4, "cm")) +
  labs(x = expression(beta ~ "(kinase activity)"), y = "Adjusted P-value (-log10)", color = "In pathway")

ggsave(filename = "kinase_pathways_associations_volcano_plot.png", plot = vplot, path = "./output/plots/kinase_TF_pathways_associations/", width = 6, height = 4)
ggsave(filename = "kinase_pathways_associations_volcano_plot.pdf", plot = vplot, path = "./output/plots/kinase_TF_pathways_associations/", width = 6, height = 4)
unlink("kinase_pathways_associations_volcano_plot.png")
unlink("kinase_pathways_associations_volcano_plot.pdf")


# adjusted P-values distribution between kinases inside/outside pathway
plot <- associations %>%
  left_join(links[, c("pathway", "gene")], by = "pathway") %>%
  mutate(p = map2_lgl(.x = kinase, .y = gene, .f = ~ if(is.null(.y)){NA}else{.x %in% .y})) %>%
  filter(!is.na(p)) %>%
  ggplot(mapping = aes(x = p, y = -log10(padj), fill = p)) +
  geom_boxplot(show.legend = F, outlier.size = 0.1, size = 0.3, fatten = 1) +
  coord_flip() +
  theme_classic() +
  theme(axis.title = element_text(size = 8), axis.text = element_text(size = 4), axis.line = element_line(size = 0.1), axis.ticks = element_line(size = 0.1)) +
  labs(x = "In pathway", y = "Adjusted P-value (-log10)")

ggsave(filename = "kinase_pathways_associations_pval_distribution.png", plot = plot, path = "./output/plots/kinase_TF_pathways_associations/", width = 2, height = 1)
ggsave(filename = "kinase_pathways_associations_pval_distribution.pdf", plot = plot, path = "./output/plots/kinase_TF_pathways_associations/", width = 2, height = 1)
unlink("kinase_pathways_associations_pval_distribution.png")
unlink("kinase_pathways_associations_pval_distribution.pdf")


# plot some examples

# load pathway activities
pathways <- data.table::fread("./data/progeny/pathway_activities_log2fpkm.csv") %>%
  as_tibble() %>%
  rename(pathway=V1) %>%
  pivot_longer(-pathway, names_to = "sample", values_to = "pathway_activity") %>%
  mutate(sample = str_replace(sample, "^X{1}", ""))


# load kinase activities
kin_activities <- data.table::fread("./data/Danish/kinaseActMatImputed.tsv") %>%
  as_tibble() %>%
  rename(kinase = V1) %>%
  pivot_longer(-kinase, names_to = "sample", values_to = "kin_activity")

plot <- associations %>%
  filter(-log10(padj) > 4) %>%
  inner_join(kin_activities, by = "kinase") %>%
  inner_join(pathways, by = c("pathway", "sample")) %>%
  inner_join(samples_annotation[, c("sample", "batch")], by = "sample") %>%
  group_by(kinase, pathway) %>%
  nest() %>%
  ungroup() %>%
  mutate(pathway_activity_resd = map(.x = data, .f = ~ residuals(lm(pathway_activity ~ batch, data = .x)))) %>%
  unnest() %>%
  mutate(pair = str_c(kinase, pathway, sep = "--")) %>%
  mutate(pair = fct_reorder(pair, padj, .fun = function(x){unique(x)})) %>%
  ggplot(mapping = aes(x = kin_activity, y = pathway_activity_resd)) +
  geom_point(size = 0.1) +
  geom_smooth(method = "lm", lwd = 0.5) +
  stat_cor(size = 2) +
  facet_wrap(~ pair, scales = "free")

ggsave(filename = "kinase_pathways_associations_examples.png", plot = plot, path = "./output/plots/kinase_TF_pathways_associations/", width = 12, height = 8)
ggsave(filename = "kinase_pathways_associations_examples.pdf", plot = plot, path = "./output/plots/kinase_TF_pathways_associations/", width = 12, height = 8)
unlink("kinase_pathways_associations_examples.png")
unlink("kinase_pathways_associations_examples.pdf")

