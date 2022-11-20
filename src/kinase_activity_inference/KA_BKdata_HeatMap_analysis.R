#' ---
#' title: "Heatmaps from kinase-activity inference using benchmark phosphorylation data"
#' author: "Abel Sousa"
#' ---


#+ global options, include=FALSE
knitr::opts_chunk$set(warning = F, message = F, collapse = T, fig.width=6, fig.height=6)



#' load R packages
library(tidyverse)
library(gplots)



#' source R kinase-activity inference pipeline
source("./src/utils/infer_kinase_activity.R")


#' load kinase activation scores
# non-weighted z-test
KAinferences <- read_tsv(file = "./output/files/BKphospho_kinase_activity_ztestNW.txt.gz")



#' # database kinase-substrate pairs
inferKA_db <- KAinferences %>%
  filter(source_type == "database") %>%
  select(-source_type)

db_hm <- kinase_heatmap(inferKA_db)

db_hm()
dev.off()

pdf(file="./output/plots/kinase_activity_inference/bk_phospho_db.pdf", height=6, width=6)
db_hm()
dev.off()



#' # text mining kinase-substrate pairs
inferKA_tm <- KAinferences %>%
  filter(source_type == "text-mining") %>%
  select(-source_type)

tm_hm <- kinase_heatmap(inferKA_tm)

tm_hm()
dev.off()

pdf(file="./output/plots/kinase_activity_inference/bk_phospho_textmining.pdf", height=6, width=6)
tm_hm()
dev.off()


#' # in vivo kinase-substrate pairs
inferKA_in_vv <- KAinferences %>%
  filter(source_type == "in_vivo") %>%
  select(-source_type)

in_vv_hm <- kinase_heatmap(inferKA_in_vv)

in_vv_hm()
dev.off()

pdf(file="./output/plots/kinase_activity_inference/bk_phospho_invv.pdf", height=6, width=6)
in_vv_hm()
dev.off()



#' # in vitro kinase-substrate pairs
inferKA_in_vt <- KAinferences %>%
  filter(source_type == "in_vitro") %>%
  select(-source_type)

in_vt_hm <- kinase_heatmap(inferKA_in_vt)

in_vt_hm()
dev.off()

pdf(file="./output/plots/kinase_activity_inference/bk_phospho_invt.pdf", height=6, width=6)
in_vt_hm()
dev.off()
