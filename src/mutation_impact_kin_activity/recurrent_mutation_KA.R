library(tidyverse)


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load kinase-activity inference data
ka_subst <- read_tsv(file = "./output/files/CPTAC_KA_phosphoQ.txt.gz")
ka_psite <- read_tsv(file = "./output/files/CPTAC_KA_kin_psites_phosphoQ.txt.gz")

k_substN <- 3
k_psiteN <- k_substN

ka_subst <- ka_subst %>%
  filter(source_type == "DB_text-mining") %>%
  filter(n >= k_substN) %>%
  select(-source_type, -n)

ka_psite <- ka_psite %>%
  filter(n >= k_psiteN) %>%
  select(-n) %>%
  rename(kinase=gene)


# load kinase-substrate lists
# select kinases
kin_sub <- read_tsv(file = "./output/files/kinase_substrate_list.txt.gz")
kinases <- kin_sub %>%
  select(kinase) %>%
  distinct() %>%
  pull(kinase)


#' load mutation data
#' calculate recurrent mutations
mutation <- read_tsv(file = "./output/files/mutations_protpos.txt.gz")

recurrent <- mutation %>%
  filter(variant_class == "Missense_Mutation" & variant_type == "SNP") %>%
  select(sample, gene_symbol, prot_pos, aa_wt, aa_mut) %>%
  distinct() %>%
  rename(gene = gene_symbol, pos = prot_pos, ref = aa_wt, alt = aa_mut) %>%
  group_by(gene, pos, ref, alt) %>%
  summarise(sample = list(sample), n = n()) %>%
  ungroup() %>%
  mutate(mut = paste(gene, pos, ref, alt, sep = "_")) %>%
  select(mut, everything()) %>%
  arrange(desc(n))

kin_recurrent_KAsub <- recurrent %>%
  filter(gene %in% kinases & n >= 4) %>%
  unnest() %>%
  inner_join(ka_subst, by = c("gene" = "kinase", "sample")) %>%
  select(-c(gene, pos, ref, alt, n)) %>%
  group_by(mut) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  filter(n>=4) %>%
  select(-n) %>%
  mutate(class = "Substrate-based Kinase Activity")

kin_recurrent_KApsite <- recurrent %>%
  filter(gene %in% kinases & n >= 4) %>%
  unnest() %>%
  inner_join(ka_psite, by = c("gene" = "kinase", "sample")) %>%
  select(-c(gene, pos, ref, alt, n)) %>%
  group_by(mut) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  filter(n>=4) %>%
  select(-n) %>%
  mutate(class = "Phophosite-based Kinase Activity")



kin_recurrent_plot <- kin_recurrent_KAsub %>%
  bind_rows(kin_recurrent_KApsite) %>%
  ggplot(mapping = aes(x = mut, y = log10P)) +
  geom_boxplot(outlier.shape = NA, lwd = 0.2) +
  geom_jitter(mapping = aes(color = mut), width = 0.1, size=0.5) +
  facet_wrap(~class) +
  coord_flip() +
  scale_color_discrete(guide=F) +
  theme(
    axis.title = element_text(colour = "black", size = 12),
    axis.text = element_text(colour = "black", size = 10),
    strip.text = element_text(colour = "black", size = 12)) +
  labs(x="", y="Kinase activation")
  

ggsave(filename = "recurrent_mut_KA_plot.png", plot = kin_recurrent_plot, path = "./output/plots/mutation_impact_kin_activity/", height = 3, width = 8)
unlink("recurrent_mut_KA_plot.png")



#' load phosphorylation data
phospho <- read_tsv(file = "./output/files/phosphoproteomicsQ_Protreg_allSamp_withZtransf.txt.gz") %>%
  pivot_longer(-c(gene, psite, psites), names_to = "sample", values_to = "log2fc") %>%
  filter(!is.na(log2fc))

#' load cptac functional score data
cptac_funscore <- read_tsv(file = "./output/files/cptac_kin_psites_funscoR.txt.gz") %>%
  mutate(psite = paste0(gene, "_", residue, position))


foo <- recurrent %>%
  filter(mut == "BRAF_600_V_E") %>%
  unnest() %>%
  inner_join(phospho, by = c("gene", "sample")) %>%
  inner_join(cptac_funscore, by = c("psite", "gene"))

foo2 <- foo %>%
  select(psite, probabilities, PSP_reg_status) %>%
  distinct()

plot <- foo %>%
  ggplot(mapping = aes(x = psite, y = log2fc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(mapping = aes(color = psite), width = 0.1) +
  geom_point(data = foo2, mapping = aes(x = psite, y = probabilities, size = probabilities), color = "red") +
  coord_flip() +
  scale_color_discrete(guide=F) +
  labs(size = "funscore")

ggsave(filename = "plot.png", plot = plot, path = "~/Downloads/", height = 3, width = 7)
unlink("plot.png")

