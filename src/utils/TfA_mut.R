# Estimate the association between TF activity and the mutational state of genes through linear modelling


library(tidyverse)
library(parallel)


linear_model <- function(df, covs = NULL, mutN = 5, NonMutN = 0, tfV = "activity", mutV = "mutated"){
  dfM <- df
  
  if(ncol(dfM) < 3) stop("Few columns. You need at least 3 columns in the data frame.")
  
  # additional covariates to regress-out
  if(!is.null(covs)){
    dfM <- dfM %>%
      inner_join(covs, by = "sample")
  }
  
  # remove sample column
  dfM <- dfM %>%
    select(-sample)
  
  # remove factors with less then 2 levels
  dfM <- dfM %>%
    mutate_at(mutV, as.numeric) %>%
    mutate_if(.predicate = is.character, .funs = as.factor) %>%
    mutate_if(.predicate = is.factor, .funs = fct_drop) %>%
    select_if(.predicate = ~ if(!is.factor(.x)){TRUE}else{if(nlevels(.x) == 1){FALSE}else{TRUE}})
  
  # calculate number of samples with and without mutations
  mut_n = table(dfM[[mutV]])
  if(length(mut_n) == 1 & "0" %in% names(mut_n)){
    mut_n = c(mut_n, setNames(0,"1"))
  }
  if(length(mut_n) == 1 & "1" %in% names(mut_n)){
    #stop("all samples are mutated!")
    mut_n = c(mut_n, setNames(0,"0"))
  }
  
  # if the number of samples with mutation is less than mutN or the number of samples without mutation is less than NonMutN return NAs
  if(!(mut_n["1"] >= mutN & mut_n["0"] >= NonMutN)){
    res <- tibble(estimate = NA, p.value = NA)
  } else {
    # define response (TF activity) and explanatory variables (mutation and other covariates)
    resp <- tfV
    expl <- rev(colnames(dfM)[colnames(dfM) != resp])
    
    # set up the formula
    f <- as.formula(paste0(resp, "~", paste(expl, collapse = "+")))
    
    # fit the model
    reg <- lm(formula = f, data = dfM)
    res <- broom::tidy(reg) %>%
      filter(term == mutV) %>%
      select(estimate, p.value)
  }
  return(res)
}


TfA_mut <- function(mutMat, tfA, covars = NULL, multi_core = T, in_tf = F, samp_mut = 20, mut_tfN = 5, NonMut_tfN = 0){
  
  # filter mutation matrix
  mutMat <- mutMat %>%
    filter(n >= samp_mut) %>%
    select(-n) %>%
    pivot_longer(cols = -c(gene), names_to = "sample", values_to = "mutated")
  
  # prepare data to model
  if(in_tf){
    dataMod <- tfA %>%
      inner_join(mutMat, by = c("tf" = "gene", "sample")) %>%
      select(tf, sample, activity, mutated) %>%
      group_by(tf) %>%
      nest() %>%
      ungroup()
  }else{
    dataMod <- tfA %>%
      #filter(tf == "AHR") %>%
      inner_join(mutMat, by = c("sample")) %>%
      select(tf, gene, sample, activity, mutated) %>%
      group_by(tf, gene) %>%
      nest() %>%
      ungroup()
  }
  
  # perform the linear modelling
  if(!multi_core){
    dataMod <- dataMod %>%
      mutate(model = map(.x = data, .f = linear_model, covs = covars, mutN = mut_tfN, NonMutN = NonMut_tfN))
  }else{
    cores = detectCores()-1
    cl <- makeCluster(cores)
    clusterEvalQ(cl, library(tidyverse))
    
    models <- parLapply(
      cl = cl,
      X = dataMod$data,
      fun = linear_model,
      covs = covars,
      mutN = mut_tfN,
      NonMutN = NonMut_tfN)
    
    stopCluster(cl)
    
    dataMod <- dataMod %>%
      mutate(model = models)
  }
  
  dataMod <- dataMod %>%
    select(-data) %>%
    unnest() %>%
    filter(!is.na(estimate)) %>%
    mutate(p.adjust = p.adjust(p.value, "BH"))
  
  return(dataMod)
}


plotTF <- function(mutMat, tfA, TF, mut_gene){
  
  mutMat <- mutMat %>%
    select(-n) %>%
    filter(gene == mut_gene) %>%
    pivot_longer(cols = -c(gene), names_to = "sample", values_to = "mutated")
  
  tfA <- tfA %>%
    filter(tf == TF)
  
  to_plot <- inner_join(tfA, mutMat, by = "sample") %>%
    select(tf, gene, everything()) %>%
    mutate(mutated = as.character(mutated))
  
  plot <- to_plot %>%
    ggplot(mapping = aes(x = mutated, y = activity)) +
    geom_boxplot(notch = TRUE, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.2) +
    ggpubr::stat_compare_means() +
    theme_classic() +
    theme(
      axis.title = element_text(colour = "black", size = 12),
      axis.text = element_text(colour = "black", size = 10)) +
    labs(x = paste(mut_gene, "mutation status", sep = " "), y = paste(TF, "activity", sep = " "))
  
  return(plot)
}
