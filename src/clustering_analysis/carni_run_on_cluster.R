library(CARNIVAL)
library(igraph)
library(readr)

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

args = commandArgs(trailingOnly=TRUE)

i <- as.numeric(args[1])

print(i)

setwd(paste0('/home/ad234505/beltrao_CPTAC/',i))

inputs_meas <- readRDS('../carni_input_meas_list.rds')

omnipath_for_carnival <- as.data.frame(read_csv("../omnipath_for_carnival.csv"))

carni_input <- inputs_meas[[i]][[1]]

carni_meas <- inputs_meas[[i]][[2]]

carni_input <- carni_input[,names(carni_input) %in% omnipath_for_carnival$source | 
                             names(carni_input) %in% omnipath_for_carnival$target]

omni_int <- downstream_neighbours(omnipath_for_carnival, 10, names(carni_input))

carni_meas <- carni_meas[,names(carni_meas) %in% omni_int$source | names(carni_meas) %in% omni_int$target]



omni_int <- upstream_neighbours(omni_int, 10, names(carni_meas))

carni_input <- carni_input[,names(carni_input) %in% omnipath_for_carnival$source | 
                             names(carni_input) %in% omnipath_for_carnival$target]

print(carni_input)
print(carni_meas)                             

carni_clust_1 <- runCARNIVAL(inputObj = carni_input, 
                             measObj = carni_meas, 
                             netObj = omni_int, 
                             solverPath = "../cplex", 
                             solver = "cplex", 
                             timelimit = 10800, 
                             mipGAP = 0.2) ##gap 0

carni_clust_1_net <- as.data.frame(carni_clust_1$weightedSIF)
carni_clust_1_att <- as.data.frame(carni_clust_1$nodesAttributes)
carni_clust_1_att <- carni_clust_1_att[carni_clust_1_att$AvgAct != 0,]
carni_clust_1_att$NodeType <- ifelse(carni_clust_1_att$NodeType == "", "intermediate",carni_clust_1_att$NodeType)

write_csv(carni_clust_1_net, "sif.csv"))
write_csv(carni_clust_1_att, "att.csv"))