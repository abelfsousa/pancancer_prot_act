#' ---
#' title: "Phosphoproteomics data"
#' author: "Abel Sousa"
#' date: "July 29th, 2019"
#' ---


#' load R packages
suppressMessages(library(limma))
suppressMessages(library(viridis))
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))



#' load samples with phosphoproteomics data
phospho_samples <- read_tsv(file = "./output/files/phosphoproteomics_samples.txt") %>%
  mutate(batch = fct_rev(fct_infreq(batch)))



#' ## number of samples in each experimental batch
samples_barplot <- ggplot(data = phospho_samples) +
  geom_bar(mapping = aes(x = batch, fill = cancer), stat = "count") +
  coord_flip() +
  theme_classic() +
  scale_fill_viridis(discrete = T) +
  theme(
    axis.title = element_text(color = "black", size = 16),
    legend.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 14)
  ) +
  scale_y_continuous(limits = c(NA, 200)) +
  labs(x = "Batch", y = "Samples", fill = "Cancer")

#+ fig.width=10, fig.height=5
samples_barplot

ggsave(filename="phospho_samples.png", plot = samples_barplot, path = "./output/plots/pre_process/", width = 10, height = 5)
unlink("phospho_samples.png")



#' load phosphoproteomics data
phospho <- fread(input = "./output/files/phosphoproteomics.txt.gz") %>%
  as_tibble() %>%
  gather(key = "sample", value = "log2FC", -gene, -psite, -psites) %>%
  inner_join(phospho_samples %>% mutate_if(is.factor, as.character), by = "sample")



#' ## number of proteins and phosphosites by experimental batch
proteins_psites_barplot <- phospho %>%
  group_by(batch, gene, psite, psites) %>%
  filter( !(sum(is.na(log2FC)) == n()) ) %>%
  group_by(batch) %>%
  summarise(Proteins = length(unique(gene)), Phosphosites = length(unique(psite))) %>%
  ungroup() %>%
  gather(key = "feature", value = "n", -batch) %>%
  mutate(batch = fct_reorder(.f = batch, .x = n, .fun = sum, .desc = FALSE)) %>%
  ggplot() +
  geom_bar(mapping = aes(x = batch, y = n, fill = feature), stat = "identity") +
  coord_flip() +
  theme_classic() +
  scale_fill_viridis(discrete = T) +
  theme(
    axis.title = element_text(color = "black", size = 16),
    legend.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 14),
    legend.text = element_text(color = "black", size = 14)
  ) +
  scale_y_continuous(limits = c(NA, 80000)) +
  labs(x = "Batch", y = "Number", fill = "Feature")

#+ fig.width=10, fig.height=5
proteins_psites_barplot

ggsave(filename="number_proteins_psites.png", plot = proteins_psites_barplot, path = "./output/plots/pre_process/", width = 10, height = 5)
unlink("number_proteins_psites.png")



#' ## log2FC distribution by sample and experimental batch
phospho_distribution_1 <- ggplot(data = phospho, mapping = aes(x = sample, y = log2FC, fill = cancer)) +
  theme_classic() +
  geom_boxplot(outlier.size = 0.1, lwd = 0.1) +
  facet_wrap(~ batch, nrow = 3, scales = "free") +
  theme(
    axis.title=element_text(colour="black", size=20),
    axis.text.y=element_text(colour="black", size=15),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.background=element_blank(),
    strip.text.x=element_text(colour="black", size=20),
    legend.position = "bottom",
    legend.title = element_text(colour="black", size=20),
    legend.text = element_text(colour="black", size=15)) +
  scale_fill_viridis(discrete = TRUE) +
  labs(x = "Sample", y = "log2 FC") +
  guides(fill = guide_legend(nrow = 2))

#+ fig.width=30, fig.height=20
phospho_distribution_1

ggsave(filename="phospho_distribution_1.png", plot = phospho_distribution_1, path = "./output/plots/pre_process/", width = 30, height = 20)
unlink("phospho_distribution_1.png")



#' ## log2FC distribution by sample and experimental batch (without outliers)
phospho_distribution_2 <- ggplot(data = phospho, mapping = aes(x = sample, y = log2FC, fill = cancer)) +
  theme_classic() +
  geom_boxplot(outlier.shape = NA, lwd = 0.1) +
  facet_wrap(~ batch, nrow = 3, scales = "free") +
  theme(
    axis.title=element_text(colour="black", size=20),
    axis.text.y=element_text(colour="black", size=15),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.background=element_blank(),
    strip.text.x=element_text(colour="black", size=20),
    legend.position = "bottom",
    legend.title = element_text(colour="black", size=20),
    legend.text = element_text(colour="black", size=15)) +
  scale_fill_viridis(discrete = TRUE) +
  scale_y_continuous(limits = c(-5,5)) +
  labs(x = "Sample", y = "log2 FC") +
  guides(fill = guide_legend(nrow = 2))

#+ fig.width=30, fig.height=20
phospho_distribution_2

ggsave(filename="phospho_distribution_2.png", plot = phospho_distribution_2, path = "./output/plots/pre_process/", width = 30, height = 20)
unlink("phospho_distribution_2.png")



#' ## summary statistics by experimental batch
summary_stats <- phospho %>%
  group_by(batch, gene, psite, psites) %>%
  filter( !(sum(is.na(log2FC)) == n()) ) %>%
  ungroup() %>%
  select(batch, sample, log2FC) %>%
  nest(data = c(log2FC)) %>%
  mutate(summary = map(.x = data, .f = ~ broom::tidy(summary(.x$log2FC)))) %>%
  select(-data) %>%
  unnest(cols=summary) %>%
  gather(key = "stat", value = "value", -c(batch, sample)) %>%
  group_by(batch, stat) %>%
  summarise(value = median(value)) %>%
  ungroup() %>%
  mutate(stat = str_replace(stat, "median", "q2")) %>%
  filter(stat != "mean") %>%
  ggplot() +
  geom_bar(mapping = aes(x = batch, y = value, fill = batch), stat = "identity", position = "dodge") +
  facet_wrap(~ stat, scales = "free") +
  coord_flip() +
  theme_classic() +
  scale_fill_viridis(discrete = T, guide=F) +
  theme(
    axis.title = element_text(color = "black", size = 16),
    axis.text.x = element_text(color = "black", size = 14, angle = 45, vjust = 0.5),
    axis.text.y = element_text(color = "black", size = 14),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 16)
  ) +
  labs(x = "Batch", y = "Median", fill = "Batch")

#+ fig.width=15, fig.height=10
summary_stats

ggsave(filename="psites_summary_stats.png", plot = summary_stats, path = "./output/plots/pre_process/", width = 15, height = 10)
unlink("psites_summary_stats.png")



#' ## number of samples whose distribution median exceeds 1/2/3 z-scores
n_outlier_samples1 <- phospho %>%
  group_by(batch, gene, psite, psites) %>%
  filter( !(sum(is.na(log2FC)) == n()) ) %>%
  ungroup() %>%
  select(batch, sample, log2FC) %>%
  nest(data = c(log2FC)) %>%
  mutate(median = map_dbl(.x = data, .f = function(x) median(x$log2FC, na.rm = TRUE))) %>%
  select(-data) %>%
  group_by(batch) %>%
  mutate(z = scale(median)[,1]) %>%
  ungroup() %>%
  mutate(higher_1z = abs(z) > 1, higher_2z = abs(z) > 2, higher_3z = abs(z) > 3) %>%
  group_by(batch) %>%
  mutate(n_samples = n()) %>%
  ungroup() %>%
  select(-sample, -median, -z) %>%
  group_by(batch, n_samples) %>%
  summarise(higher_1z = sum(higher_1z), higher_2z = sum(higher_2z), higher_3z = sum(higher_3z)) %>%
  ungroup() %>%
  gather(key = "times_higher", value = "count", -c("batch", "n_samples")) %>%
  arrange(batch)

n_outlier_samples_plot1 <- n_outlier_samples1 %>%
  ggplot() +
  geom_bar(mapping = aes(x = batch, y = count, fill = times_higher), stat = "identity", position = "dodge") +
  annotate(
    "label",
    label = n_outlier_samples1 %>% select(batch, n_samples) %>% distinct() %>% pull(n_samples) %>% paste("samples", ., sep=": "),
    y = rep(35,10),
    x = 1:10,
    size=3) +
  coord_flip() +
  theme_classic() +
  scale_fill_viridis(discrete = T) +
  theme(
    axis.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 14),
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 14)
  ) +
  labs(x = "Batch", y = "Count", fill = "Times higher")

#+ fig.width=8, fig.height=6
n_outlier_samples_plot1

ggsave(filename="phospho_n_outlier_samples1.png", plot = n_outlier_samples_plot1, path = "./output/plots/pre_process/", width = 8, height = 6)
unlink("phospho_n_outlier_samples1.png")





#' ## number of samples whose distribution median exceeds log2(2)/(3)/(4)
n_outlier_samples2 <- phospho %>%
  group_by(batch, gene, psite, psites) %>%
  filter( !(sum(is.na(log2FC)) == n()) ) %>%
  ungroup() %>%
  select(batch, sample, log2FC) %>%
  nest(data = c(log2FC)) %>%
  mutate(median = map_dbl(.x = data, .f = function(x) median(x$log2FC, na.rm = TRUE))) %>%
  select(-data) %>%
  mutate(median = abs(median)) %>%
  mutate(higher_2_fc = median > log2(2), higher_3_fc = median > log2(3), higher_4_fc = median > log2(4)) %>%
  group_by(batch) %>%
  mutate(n_samples = n()) %>%
  ungroup() %>%
  select(-sample, -median) %>%
  group_by(batch, n_samples) %>%
  summarise(higher_2_fc = sum(higher_2_fc), higher_3_fc = sum(higher_3_fc), higher_4_fc = sum(higher_4_fc)) %>%
  ungroup() %>%
  gather(key = "times_higher", value = "count", -c("batch", "n_samples")) %>%
  arrange(batch)



n_outlier_samples_plot2 <- n_outlier_samples2 %>%
  ggplot() +
  geom_bar(mapping = aes(x = batch, y = count, fill = times_higher), stat = "identity", position = "dodge") +
  annotate(
    "label",
    label = n_outlier_samples2 %>% select(batch, n_samples) %>% distinct() %>% pull(n_samples) %>% paste("samples", ., sep=": "),
    y = rep(20,10),
    x = 1:10,
    size=3) +
  coord_flip() +
  theme_classic() +
  scale_fill_viridis(discrete = T) +
  theme(
    axis.title = element_text(color = "black", size = 16),
    axis.text = element_text(color = "black", size = 14),
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 14)
  ) +
  labs(x = "Batch", y = "Count", fill = "Times higher")

#+ fig.width=8, fig.height=6
n_outlier_samples_plot2

ggsave(filename="phospho_n_outlier_samples2.png", plot = n_outlier_samples_plot2, path = "./output/plots/pre_process/", width = 8, height = 6)
unlink("phospho_n_outlier_samples2.png")



#' ### A: outlier samples with median distribution z-score > 2
samples_z2 <- phospho %>%
  group_by(batch, gene, psite, psites) %>%
  filter( !(sum(is.na(log2FC)) == n()) ) %>%
  ungroup() %>%
  select(batch, sample, log2FC) %>%
  nest(data=log2FC) %>%
  mutate(median = map_dbl(.x = data, .f = function(x) median(x$log2FC, na.rm = TRUE))) %>%
  select(-data) %>%
  group_by(batch) %>%
  mutate(z = scale(median)[,1]) %>%
  ungroup() %>%
  filter(abs(z) > 2) %>%
  pull(sample)

length(samples_z2)


#' ### B: outlier samples with median distribution > 2 fold-change = log2(2) = 1
samples_fc2 <- phospho %>%
  group_by(batch, gene, psite, psites) %>%
  filter( !(sum(is.na(log2FC)) == n()) ) %>%
  ungroup() %>%
  select(batch, sample, log2FC) %>%
  nest(data = log2FC) %>%
  mutate(median = map_dbl(.x = data, .f = function(x) median(x$log2FC, na.rm = TRUE))) %>%
  select(-data) %>%
  mutate(median = abs(median)) %>%
  filter(median > log2(2)) %>%
  pull(sample)

length(samples_fc2)



#' ## remove outlier samples in B (median distribution > 2 fold-change = log2(2) = 1)
phospho <- phospho %>%
  filter( !(sample %in% samples_fc2) )


phospho_samples <- phospho_samples %>%
  filter( !(sample %in% samples_fc2) )

write.table(phospho_samples, "./output/files/phosphoproteomics_samples.txt", sep="\t", quote=F, row.names=F)



#' ## log2FC distribution by sample and experimental batch (without outliers)
phospho_distribution_3 <- ggplot(data = phospho, mapping = aes(x = sample, y = log2FC, fill = cancer)) +
  theme_classic() +
  geom_boxplot(outlier.shape = NA, lwd = 0.1) +
  facet_wrap(~ batch, nrow = 3, scales = "free") +
  theme(
    axis.title=element_text(colour="black", size=20),
    axis.text.y=element_text(colour="black", size=15),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.background=element_blank(),
    strip.text.x=element_text(colour="black", size=20),
    legend.position = "bottom",
    legend.title = element_text(colour="black", size=20),
    legend.text = element_text(colour="black", size=15)) +
  scale_fill_viridis(discrete = TRUE) +
  scale_y_continuous(limits = c(-5,5)) +
  labs(x = "Sample", y = "log2 FC") +
  guides(fill = guide_legend(nrow = 2))

#+ fig.width=30, fig.height=20
phospho_distribution_3

ggsave(filename="phospho_distribution_3.png", plot = phospho_distribution_3, path = "./output/plots/pre_process/", width = 30, height = 20)
unlink("phospho_distribution_3.png")


phospho <- phospho %>%
  select(-batch, -cancer) %>%
  spread(key = "sample", value = "log2FC")


gz1 <- gzfile("./output/files/phosphoproteomics.txt.gz", "w")
write.table(phospho, gz1, sep="\t", quote=F, row.names=F)
close(gz1)



#' ## quantile normalization using normalizeQuantiles limma function
phosphoQ_psites <- phospho[, c(1:3)]
phosphoQ <- normalizeQuantiles( as.matrix(phospho[,-c(1:3)]) )

phosphoQ <- phosphoQ %>%
  as_tibble() %>%
  bind_cols(phosphoQ_psites) %>%
  select(gene, psite, psites, everything())

gz2 <- gzfile("./output/files/phosphoproteomicsQ.txt.gz", "w")
write.table(phosphoQ, gz2, sep="\t", quote=F, row.names=F)
close(gz2)



#' ## log2FC distribution by sample and experimental batch (without outliers)
phospho_distribution_4 <- phosphoQ %>%
  gather(key = "sample", value = "log2FC", -gene, -psite, -psites) %>%
  inner_join(phospho_samples %>% mutate_if(is.factor, as.character), by = "sample") %>%
  ggplot(mapping = aes(x = sample, y = log2FC, fill = cancer)) +
  theme_classic() +
  geom_boxplot(outlier.shape = NA, lwd = 0.1) +
  facet_wrap(~ batch, nrow = 3, scales = "free") +
  theme(
    axis.title=element_text(colour="black", size=20),
    axis.text.y=element_text(colour="black", size=15),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    strip.background=element_blank(),
    strip.text.x=element_text(colour="black", size=20),
    legend.position = "bottom",
    legend.title = element_text(colour="black", size=20),
    legend.text = element_text(colour="black", size=15)) +
  scale_fill_viridis(discrete = TRUE) +
  scale_y_continuous(limits = c(-5,5)) +
  labs(x = "Sample", y = "log2 FC") +
  guides(fill = guide_legend(nrow = 2))

#+ fig.width=30, fig.height=20
phospho_distribution_4

ggsave(filename="phospho_distribution_4.png", plot = phospho_distribution_4, path = "./output/plots/pre_process/", width = 30, height = 20)
unlink("phospho_distribution_4.png")
