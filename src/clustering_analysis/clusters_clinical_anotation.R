library(readr)
library(reshape2)
library(piano)
library(pheatmap)

clinical_data <- as.data.frame(
  read_delim("data/updated/cptac_data/clinical/clinical_data.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE))
clinical_data$bmi <- ifelse(clinical_data$bmi >= 37, "big","smol")
clinical_data <- clinical_data[,-c(12,14,15)]
clinical_data$age <- ifelse(clinical_data$age > 50, "young","venerable")


clin_dat_long <- melt(clinical_data, id.vars = "sample")
clin_dat_long$set <- paste(clin_dat_long$variable, clin_dat_long$value, sep = "_")


test <- table(clin_dat_long$set)
test <- test[test > 5]

clin_dat_long <- clin_dat_long[clin_dat_long$set %in% names(test),]
  
clin_dat_long <- clin_dat_long[,c(1,4)]

sample_clusters <- as.data.frame(read_csv("support/sample_clusters_bcor.csv"))

ORA_res_list <- list()
i <- 1
for(cluster in unique(sample_clusters$condition))
{
  sig_samples <- sample_clusters[sample_clusters$condition == cluster, "sample"]
  
  ORA_res_list[[i]] <- runGSAhyper(sig_samples, universe = sample_clusters$sample, gsc = loadGSC(clin_dat_long),adjMethod = "fdr")
  i <- i+1
}

test <- ORA_res_list[[8]]$resTab
test$pathway <- row.names(ORA_res_list[[8]]$resTab)

View(ORA_res_list[[8]]$resTab)


View(table(sample_clusters$condition))

i <- 1
for(ORA_res in ORA_res_list)
{
  write.csv(ORA_res$resTab, paste0("~/Dropbox/beltrao_CPTAC/results/cluster_clinical_ORA_C",i,".csv"))
  i <- i+1
}


ORA_res_list_todf <- sapply(ORA_res_list, function(x){
  df <- as.data.frame(x$resTab)
  View(df)
  df$pathway <- row.names(df)
  df <- df[,c(7,1)]
  View(df)
  return(df)
}, simplify = F)

ORA_res_df <- merge(ORA_res_list_todf[[1]],ORA_res_list_todf[[2]], by = "pathway")
for(i in 3:length(ORA_res_list_todf))
{
  ORA_res_df <- merge(ORA_res_df, ORA_res_list_todf[[i]], by = "pathway")
}

names(ORA_res_df) <- c("clinical_feature",paste0('cluster_',1:length(ORA_res_list)))
row.names(ORA_res_df) <- ORA_res_df[,1]
ORA_res_df <- ORA_res_df[,-1]
ORA_res_df[ORA_res_df > 0.1] <- NA
ORA_res_df <- ORA_res_df[rowSums(is.na(ORA_res_df)) < 8,]
pheatmap(ORA_res_df, cluster_cols = F, cluster_rows = F, na_col = "grey", colorRampPalette(c("blue", "lightblue"))(50))

top_n <- 10
for(i in 1:8)
{
  if(i == 1)
  {
    ora_res_top <- as.data.frame(row.names(ORA_res_df[order(ORA_res_df[,i], decreasing = F),])[1:top_n])
  } else
  {
    ora_res_top <- as.data.frame(cbind(ora_res_top,row.names(ORA_res_df[order(ORA_res_df[,i], decreasing = F),])[1:top_n]))
  }
}
names(ora_res_top) <- 1:8
