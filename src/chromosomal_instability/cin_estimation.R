# Estimation of chromosomal instability


# load R packages
library(tidyverse)


source("./src/utils/getSamples.R")


# set nest and unnest functions to the previous tidyr version
# current versions (tidyr 1.0.0) are very slow
nest <- nest_legacy
unnest <- unnest_legacy


#' load samples annotation
samples_annotation <- read_tsv(file = "./output/files/all_samples_annotation.txt")
metadata <- getSamples(samples_annotation, c("cnv", "protein"))


# load CNV data (GISTIC2 scores)
cnv <- read_tsv("./output/files/cnv.txt.gz") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "score")

cin <- cnv %>%
  group_by(sample) %>%
  summarise(cin = sum(score == 2 | score == -2)) %>%
  ungroup() %>%
  mutate(scaled = (cin - min(cin))/(max(cin)-min(cin)))

write_tsv(cin, "./output/files/cin.txt.gz")
