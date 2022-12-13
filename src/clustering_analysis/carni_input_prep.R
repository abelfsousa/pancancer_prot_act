library(readr)
library(readr)
library(viper)
library(pheatmap)
library(CARNIVAL)
library(scales)
library(dorothea)
library(igraph)

#optional clustering
library(factoextra)
library(NbClust)

downstream_neighbours <- function(meta_network, n_steps, input_names)
{
  meta_g <- graph_from_data_frame(meta_network[,c(1,3,2)]) 
  
  dn_nbours <- ego(graph = meta_g, order = n_steps,nodes = input_names,mode = "out")
  
  sub_nodes <- unique(names(unlist(dn_nbours)))
  
  meta_network <- meta_network[meta_network$source %in% sub_nodes & meta_network$target %in% sub_nodes,]
  
  return(meta_network)
}

upstream_neighbours <- function(meta_network, n_steps, input_names)
{
  meta_g <- graph_from_data_frame(meta_network[,c(1,3,2)]) 
  
  dn_nbours <- ego(graph = meta_g, order = n_steps,nodes = input_names,mode = "in")
  
  sub_nodes <- unique(names(unlist(dn_nbours)))
  
  meta_network <- meta_network[meta_network$source %in% sub_nodes & meta_network$target %in% sub_nodes,]
  
  return(meta_network)
}

kinaseActMatImputed <- as.data.frame(
  read_delim("~/Dropbox/beltrao_CPTAC/data/updated/cptac_data/kinase_activities/kinaseActMatImputed_batch_regOut.txt", 
             "\t", escape_double = FALSE, trim_ws = TRUE))




# TFactivity <- as.data.frame(
#   read_delim("~/Dropbox/beltrao_CPTAC/data/updated/cptac_data/rna/TF_activity_log2FC_batch_regOut.txt",
#              "\t", escape_double = FALSE, trim_ws = TRUE)
# )

#from PKG_TF_activity_log2FC_allsamples.csv, regout by Abel
TFactivity <- as.data.frame(
  read_delim("~/Dropbox/beltrao_CPTAC/results/TF_activity_log2FC_batch_regOut.txt", 
             "\t", escape_double = FALSE, trim_ws = TRUE)
)

# url <- paste0(
#   'http://omnipathdb.org/interactions?',
#   'datasets=tfregulons&tfregulons_levels=A,B,C&genesymbols=1&fields=sources,tfregulons_level'
# )
# 
download_omnipath <- function(){
  
  read.table(url, sep = '\t', header = TRUE)
  
}

##Dorothea/viper
# dorothea <- download_omnipath()
# dorothea <- dorothea[,c(4,3,6,7)]
# dorothea$sign <- dorothea$is_stimulation - dorothea$is_inhibition
# dorothea <- dorothea[dorothea$sign != 0,]
# dorothea <- dorothea[,c(1,2,5)]
# dorothea <- as.data.frame(dorothea_hs[dorothea_hs$confidence %in% c("A","B","C"),c(3,1,4)])
dorothea <- as.data.frame(read_csv("~/Dropbox/beltrao_CPTAC/support/dorothea_ABC.csv"))

names(kinaseActMatImputed)[1] <- "ID"
names(TFactivity)[1] <- "ID"

kinaseActMatImputed <- kinaseActMatImputed[,names(kinaseActMatImputed) %in% names(TFactivity),]
TFactivity <- TFactivity[,names(TFactivity) %in% names(kinaseActMatImputed)]

activity_df_combined <- as.data.frame(rbind(kinaseActMatImputed, TFactivity))
row.names(activity_df_combined) <- activity_df_combined$ID
activity_df_combined <- activity_df_combined[,-1]

write.csv(activity_df_combined, "~/Dropbox/beltrao_CPTAC/data/updated/cptac_data/activity_df_combined_batch_corrected.csv")

### Clustering

pheatmap(cor(activity_df_combined, method = "spearman"), filename = "~/Dropbox/beltrao_CPTAC/visualisation/hclust_cormat_bcor.pdf", height = 100, width = 100)
pheatmap(cor(activity_df_combined, method = "spearman"), filename = "~/Dropbox/beltrao_CPTAC/visualisation/hclust_cormat_nolab_bcor.pdf",show_rownames = F, show_colnames = F) 
pheatmap(cor(activity_df_combined, method = "spearman"), filename = "~/Dropbox/beltrao_CPTAC/visualisation/hclust_cormat_nolab_bcor.png",show_rownames = F, show_colnames = F) 

pheatmap(cor(activity_df_combined, method = "spearman"), filename = "~/Dropbox/beltrao_CPTAC/visualisation/hclust_cormat_nolab_bcor_cutree8.png",show_rownames = F, show_colnames = F, cutree_rows = 8) 

cor_mat <- as.data.frame(cor(activity_df_combined, method = "spearman"))
hclust_result <- hclust(dist(cor_mat))

clusters <- cutree(hclust_result, k = 8)

# alternative clustering (these give small number of clusters, but we want to stratify more)
# Elbow method
df <- as.data.frame(t(activity_df_combined))
fviz_nbclust(df, kmeans, method = "wss", k.max = 30) +
  geom_vline(xintercept = 8, linetype = 2)+
  labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(df, kmeans, method = "silhouette", k.max = 30)+
  labs(subtitle = "Silhouette method")
# 
# # Gap statistic
# # nboot = 50 to keep the function speedy. 
# # recommended value: nboot= 500 for your analysis.
# # Use verbose = FALSE to hide computing progression.
# set.seed(123)
# fviz_nbclust(df, kmeans, nstart = 25,  method = "gap_stat", nboot = 50, k.max = 30)+
#   labs(subtitle = "Gap statistic method")

### This is were the cluster - sample association can be found
targets <- as.data.frame(matrix(NA,length(activity_df_combined[1,]),2))
names(targets) <- c("sample","condition")
targets$sample <- names(activity_df_combined)
targets$condition <- paste0("C",clusters)

cluster_count <- table(targets$condition)

View(cluster_count)

write_csv(targets, "~/Dropbox/beltrao_CPTAC/support/sample_clusters_bcor.csv")

cluster_level_activities <- list()
i <- 1
for(cluster in unique(targets$condition))
{
  df <- activity_df_combined[,names(activity_df_combined) %in% targets[targets$condition == cluster,"sample"]]
  meds <- rowMedians(as.matrix(df))
  sds <- apply(df,1,sd)
  med_over_sd <- meds/sds
  cluster_level_activities[[i]] <- med_over_sd
  i <- i+1
}

cluster_level_activities <- as.data.frame(do.call(cbind, cluster_level_activities))

cluster_level_activities_TF <- cluster_level_activities[row.names(cluster_level_activities) %in% dorothea$tf,]

pheatmap(cluster_level_activities)

# plot(density(cluster_level_activities[,10]))

write.csv(cluster_level_activities,"~/Dropbox/beltrao_CPTAC/results/cluster_level_activities_bcor.csv")

cluster_level_activities_rank <- cluster_level_activities

for(i in 1:8)
{
  cluster_level_activities_rank[,i] <- length(cluster_level_activities_rank[,i]) - rank(abs(cluster_level_activities_rank[,i])) + 1
}

cluster_level_activities_top <- cluster_level_activities
cluster_level_activities_top[cluster_level_activities_rank > 5] <- NA
cluster_level_activities_top <- cluster_level_activities_top[rowSums(is.na(cluster_level_activities_top)) < 8, ]
names(cluster_level_activities_top) <- gsub("V","Cluster_", names(cluster_level_activities_top))

pheatmap(cluster_level_activities_top, cluster_cols = F, cluster_rows = F, na_col = "white", colorRampPalette(c("blue", "white", "red"))(50))

### CARNIVAL ####

url <- "http://omnipathdb.org/interactions?types=post_translational,transcriptional&datasets=omnipath,pathwayextra,dorothea&fields=sources,references,curation_effort,dorothea_level,type&genesymbols=yes"
# url <- "http://omnipathdb.org/interactions?types=post_translational,transcriptional&datasets=omnipath,dorothea&fields=sources,references,curation_effort,dorothea_level,type&genesymbols=yes"


omni_interactions <- download_omnipath()
omni_interactions$sign <- omni_interactions$consensus_stimulation - omni_interactions$consensus_inhibition
omni_interactions <- omni_interactions[omni_interactions$sign != 0,]
omni_interactions <- omni_interactions[,c(3,16,4)]

cluster_level_activities <- cluster_level_activities[(row.names(cluster_level_activities) %in% omni_interactions$source_genesymbol | row.names(cluster_level_activities) %in% omni_interactions$target_genesymbol),]

cluster_level_activities_TF <- cluster_level_activities[row.names(cluster_level_activities) %in% dorothea$tf,]
cluster_level_activities_kinases <- cluster_level_activities[!(row.names(cluster_level_activities) %in% dorothea$tf),]

carni_input_meas_list <- list()
for(i in 1:length(cluster_level_activities_kinases[1,]))
{
  print(i)
  carni_input_1 <- cluster_level_activities_kinases[order(abs(cluster_level_activities_kinases[,i]), decreasing = T),]
  carni_input_1 <- carni_input_1[,i,drop =F]
  print('plop')
  carni_input_1[,1] <- (carni_input_1[,1] - mean(carni_input_1[,1])) / sd(carni_input_1[,1])
  plot(density(carni_input_1[,1]))
  print('plop')
  carni_input_1 <- carni_input_1[abs(carni_input_1[,1]) > 2,,drop = F]
  # carni_input_1 <- carni_input_1[1:5,,drop = F]
  print(dim(carni_input_1))
  carni_input_1 <- as.data.frame(t(carni_input_1))
  carni_input_1 <- sign(carni_input_1)
  names(carni_input_1) <- gsub("[-,;.()]",'__',names(carni_input_1))
  
  carni_meas_1 <- cluster_level_activities_TF[order(abs(cluster_level_activities_TF[,i]), decreasing = T),]
  carni_meas_1 <- carni_meas_1[,i,drop =F]
  carni_meas_1[,1] <- (carni_meas_1[,1] - mean(carni_meas_1[,1])) / sd(carni_meas_1[,1])
  plot(density(carni_meas_1[,1]))
  carni_meas_1 <- carni_meas_1[abs(carni_meas_1[,1]) > 2,,drop = F]
  # carni_meas_1 <- carni_meas_1[1:12,,drop = F]
  print(dim(carni_meas_1))
  carni_meas_1 <- as.data.frame(t(carni_meas_1))
  names(carni_meas_1) <- gsub("[-,;.()]",'__',names(carni_meas_1))
  
  carni_input_meas_list[[i]] <- list()
  carni_input_meas_list[[i]][[1]] <- carni_input_1
  carni_input_meas_list[[i]][[2]] <- carni_meas_1
}

saveRDS(carni_input_meas_list,'~/Dropbox/beltrao_CPTAC/data/cptac_data/final_carni_input/carni_input_meas_list.rds')

carni_input_1 <- cluster_level_activities_kinases[order(abs(cluster_level_activities_kinases$V1), decreasing = T),]
carni_input_1 <- carni_input_1[,1,drop =F]
carni_input_1$V1 <- (carni_input_1$V1 - mean(carni_input_1$V1)) / sd(carni_input_1$V1)
plot(density(carni_input_1$V1))
carni_input_1 <- carni_input_1[abs(carni_input_1$V1) > 2,,drop = F]
carni_input_1 <- as.data.frame(t(carni_input_1))
carni_input_1 <- sign(carni_input_1)

carni_meas_1 <- cluster_level_activities_TF[order(abs(cluster_level_activities_TF$V1), decreasing = T),]
carni_meas_1 <- carni_meas_1[,1,drop =F]
carni_meas_1$V1 <- (carni_meas_1$V1 - mean(carni_meas_1$V1)) / sd(carni_meas_1$V1)
plot(density(carni_meas_1$V1))
carni_meas_1 <- carni_meas_1[abs(carni_meas_1$V1) > 2,,drop = F]
carni_meas_1 <- as.data.frame(t(carni_meas_1))

names(omni_interactions) <- c("source","Effect","target")
omni_interactions$source <- gsub("[-,;.()]",'__',omni_interactions$source)
omni_interactions$target <- gsub("[-,;.()]",'__',omni_interactions$target)

write_csv(omni_interactions, '~/Dropbox/beltrao_CPTAC/support/omnipath_for_carnival.csv')
