library(tidyverse)
library(parallel)
library(Hmisc)


adjacencyMat_rcorr <- function(samples, datExpr, size = 10, type = "pearson"){
  
  # select respective samples
  datExpr <- datExpr[, c("gene", samples)]
  
  # select genes with expression in at least "size" samples
  datExpr <- datExpr %>%
    pivot_longer(cols = -gene, names_to = "sample", values_to = "value") %>%
    group_by(gene) %>%
    filter(sum(!is.na(value)) >= size) %>%
    ungroup() %>%
    pivot_wider(names_from = "sample", values_from = "value")
  
  # transpose the expression data
  datExpr <- datExpr %>%
    as.data.frame() %>%
    column_to_rownames(var = "gene") %>%
    t()
  
  # calculate the adjacency matrix
  adjMat <- rcorr(x = datExpr, type = type)
  
  return(adjMat)
}


adjacency_mapper <- function(expression, metadata, group = FALSE, group_var = NULL, multi_core = FALSE, size = 10, type = "pearson"){
  
  expr_data <- expression
  
  if(group){
    variable <- group_var
    
    dataF <- metadata %>%
      group_by(!! sym(variable)) %>%
      summarise(samples = list(sample)) %>%
      ungroup()
  }else{
    dataF <- metadata %>%
      summarise(samples = list(sample))
  }
  
  if(!multi_core){
    matrices <- dataF %>%
      mutate(matrices = map(.x = samples, .f = adjacencyMat_rcorr, datExpr = expr_data, size = size, type = type))
  }else{
    cores = detectCores()-1
    cl <- makeCluster(cores, type="FORK")
    #cl <- makeCluster(cores, type="PSOCK")
    #clusterExport(cl, c("expr_data", "type"))
    #clusterEvalQ(cl, library(tidyverse))
    mat_list <- parLapply(
      cl = cl,
      X = dataF$samples,
      fun =  adjacencyMat_rcorr, datExpr = expr_data, size = size, type = type)
    stopCluster(cl)
    
    matrices <- dataF %>%
      mutate(matrices = mat_list)
  }
  return(matrices)
}


# set up a function to convert a squared matrix into a data frame
matrix_to_tb <- function(x, y, nm){
  
  mat <- x
  mat_data <- y
  df_names <- nm
  
  mat[upper.tri(mat)] <- Inf
  df <- mat %>%
    as.data.frame() %>%
    rownames_to_column(var=df_names[1]) %>%
    as_tibble() %>%
    pivot_longer(-c(df_names[1]), names_to=df_names[2], values_to=mat_data) %>%
    #filter(!is.infinite(!! rlang::syms(mat_data)[[1]]))
    filter_at(.vars = mat_data, .vars_predicate = all_vars(!is.infinite(.)))
  
  df
}


# set up a function to join the rcorr output into a single tibble
joinMatrices <- function(x, nm = c("protein_A", "protein_B")){
  
  rcorr_matrices <- x
  rcorr_names <- names(rcorr_matrices)
  
  df_names <- nm
  
  dfs <- map2(.x = rcorr_matrices, .y = rcorr_names, .f = matrix_to_tb, nm = df_names)
  
  tb <- reduce(.x = dfs, .f = ~ inner_join(.x, .y, by = df_names))
  
  tb <- tb %>%
    filter(!is.na(P))
  
  tb
}
