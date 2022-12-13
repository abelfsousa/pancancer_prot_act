library(readr)
library(dorothea)
library(omicToolsTest)
library(viper)

dorothea <- as.data.frame(dorothea_hs[dorothea_hs$confidence %in% c("A","B","C"),])

dorothea_viper <- df_to_viper_regulon(dorothea[,c(3,1,4)])

transcriptomics_log2fc <- as.data.frame(read_delim("~/Dropbox/beltrao_CPTAC/data/updated/cptac_data/rna/transcriptomics_log2fpkm.txt", 
                                                   "\t", escape_double = FALSE, trim_ws = TRUE))
names(transcriptomics_log2fc)[1] <- "ID"
row.names(transcriptomics_log2fc) <- transcriptomics_log2fc$ID
transcriptomics_log2fc <- transcriptomics_log2fc[,-1]

TFactivity <- as.data.frame(viper(eset = transcriptomics_log2fc, regulon = dorothea_viper, minsize = 5, adaptive.size = F, eset.filter = F))

write.csv(TFactivity,"~/Dropbox/beltrao_CPTAC/results/PKG_TF_activity_FPKMlog2FC_allsamples.csv")
