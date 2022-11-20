# implementation of the z-test method for kinase activity inference
# based on claudia's method
# pubmed: 28200105, 23532336

# abelsousa


library(tidyverse)
library(parallel)


z_test <- function(all_psites, kin_psites){
  
  muBg <- mean(all_psites, na.rm = T)
  sdBg <- sd(all_psites, na.rm = T)
  
  X <- mean(kin_psites, na.rm = T)
  N <- length(kin_psites)
  
  z_score <- ((X - muBg) * sqrt(N))/sdBg
  
  p_value <- 2*pnorm(-abs(z_score))
  if(p_value == 0) p_value <- 0.0001
  
  if (z_score < 0) activity <- log10(p_value)
  if (z_score > 0) activity <- -log10(p_value)
  #if (X < 0) activity <- log10(p_value)
  #if (X > 0) activity <- -log10(p_value)
  
  results <- tibble(n=N, log10P=activity)
  
  return(results)
}


infer_by_kinase <- function(kinase, psites, phospho, test){
  
  # get the phosphosites of a given kinase that were quantified on a given sample
  psites <- psites %>%
    mutate(gene = kinase) %>%
    select(gene, everything()) %>%
    inner_join(phospho, by = c("gene", "position", "residue"))
  
  # remove the kinase phosphosites from the all phosphosites quantified
  phospho <- phospho %>%
    anti_join(psites, by = c("gene", "position", "residue"))
  
  # get the statistical test
  stat_test <- get(test)
  
  # if no psites were quantified for this kinase return NA
  if(nrow(psites) == 0){
    p <- tibble(n=0, log10P=NA)
  }else{
    p <- stat_test(phospho$log2fc, psites$log2fc)
  }
  return(p)
}


infer_by_sample <- function(phospho, kin_list, test){
  
  phosphoS <- phospho
  
  # infer the activation of each kinase on a given sample
  ka <- kin_list %>%
    mutate(data2 = map2(.x = gene, .y = data, .f=infer_by_kinase, phospho=phosphoS, test = test)) %>%
    select(-data) %>%
    unnest(cols = data2)
  
  return(ka)
}


infer_kinase_activity <- function(phospho, kin_list, test = "z_test", multi_core=F){
  
  # tidy phosphorylation data
  # remove NAs
  # for each sample build a tibble with all proteins, positions, residues and log2FC
  # store them in a list-column
  phospho <- phospho %>%
    pivot_longer(names_to = "sample", values_to = "log2fc", -c(position, gene, residue)) %>%
    select(sample, gene, everything()) %>%
    filter(!is.na(log2fc)) %>%
    group_by(sample) %>%
    nest() %>%
    ungroup()
  
  # for each kinase build a tibble with all the phosphorylated positions and residues
  # store them in a list-column
  kin_list <- kin_list %>%
    nest(data = c(position, residue))
  
  # infer kinase activation across samples
  if(!multi_core){
    inferKA <- phospho %>%
      mutate(data2 = map(.x = data, .f = infer_by_sample, kin_list = kin_list, test = test))
  }else{
    cores = detectCores()-1
    cl <- makeCluster(cores)
    clusterExport(cl, c("infer_by_kinase", "z_test"))
    clusterEvalQ(cl, library(tidyverse))
    data2 <- parLapply(
      cl = cl,
      X = phospho$data,
      fun = infer_by_sample,
      kin_list = kin_list,
      test = test)
    stopCluster(cl)
    inferKA <- phospho %>%
      mutate(data2 = data2)
  }
  inferKA <- inferKA %>%
    select(-data) %>%
    unnest(cols = data2)
  return(inferKA)
}

