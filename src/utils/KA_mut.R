# Estimate the association between kinase-activity and the mutational state of genes through linear modelling


library(tidyverse)
library(parallel)
library(ggpubr)
library(lmtest)


#' function to compute likelihood ratio tests between two nested linear models
#' @param df must a data.frame or tibble
#' @param dep the dependent variable, ind the independent variable(s), ind_drop the independent variable(s) to drop out
#' @param mutP minimal percentage of samples with mutations to perform the association (10% by default)
#' @param mutN minimal number of samples with mutations to perform the association (5 by default)
lr_test <- function(df, covs, dep, ind, ind_drop, mutP = 0.01, mutN = 5){
  
  dfM <- df
  # additional covariates to include in the model
  if(!is.null(covs)){
    dfM <- dfM %>%
      inner_join(covs, by = "sample")
  }
  
  # remove NAs from the dataset
  dfM <- na.exclude(dfM)
  
  # remove sample column
  dfM <- dfM %>%
    select(-sample)
  
  # calculate percentage of samples with and without mutations
  mut_n = table(dfM[[2]])
  if(length(mut_n) == 1 & "0" %in% names(mut_n)) mut_n = c(mut_n, setNames(0,"1"))
  if(length(mut_n) == 1 & "1" %in% names(mut_n)) stop("all samples are mutated!")
  
  mut_prop = mut_n/sum(mut_n)
  
  if(!(mut_prop["1"] >= mutP & mut_n["1"] >= mutN)){
    res <- tibble(estimate = NA, p.value = NA)
  }else{
    ind <- ind[!ind %in% ind_drop]
    
    fl <- as.formula(str_c(dep, str_c(ind, collapse = "+"), sep = "~"))
    lmod1 <- lm(fl, dfM)
    fl <- update.formula(fl, str_c("~ . +", str_c(ind_drop, collapse = "+")))
    lmod2 <- lm(fl, dfM)
    
    lrt <- lrtest(lmod2, lmod1)
    
    res <- broom::tidy(lmod2) %>%
      filter(term == ind_drop) %>%
      select(estimate) %>%
      mutate(p.value = lrt[2,5])
  }
  return(res)
}




#' function to fit a linear model between kinase activity and mutational status of cancer samples
#' @param df data.frame with at least 3 columns:
#' the first with the sample identifiers, the second the kinase activity and the third the mutational status (0 or 1)
#' @param covs other covariates to be passed. Covs is a data.frame with one column "sample" containing the same sample identifiers as df
#' @param mutN minimal number of samples with mutations to perform the association (5 by default)
linear_model <- function(df, covs = NULL, mutN = 5, NonMutN = 0, kaV = "log10P", mutV = "mutated"){
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
    # define response (KA) and explanatory variables (mutation and other covariates)
    resp <- kaV
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






#' function to define the kinase-gene pairs that will be used to fit linear models
#' @param mutMat mutation matrix
#' @param ka kinase-activity inference data
#' @param multi_core if TRUE the processes will be parallelized (socket approach)
#' @param unsignPval if TRUE the absolute value of the kinase activity score will be taken
#' @param in_kinase if TRUE the associations will be performed with the mutational status of the same kinase. All kinases against all genes otherwise
#' @param samp_mut minimal number of samples with mutations
#' @param kinase_list list to select the kinases from
#' @param k_subN minimal number of kinase substrates
#' @param k_sampN minimal number of samples with inference of kinase activity
#' @param mut_kaN minimal number of samples with mutations required to perform associations
KA_mut <- function(mutMat, ka, covars = NULL, multi_core = T, unsignPval = F, in_kinase = F, samp_mut = 20, kinase_list = "DB_text-mining", k_subN = 3, k_sampN = 10, mut_kaN = 5, NonMut_kaN = 0){
  
  # filter mutation matrix
  mutMat <- mutMat %>%
    filter(n >= samp_mut) %>%
    select(-n) %>%
    pivot_longer(cols = -c(gene), names_to = "sample", values_to = "mutated")

  # filter kinase-activity inference data
  ka <- ka %>%
    filter(source_type %in% kinase_list) %>%
    filter(n >= k_subN) %>%
    select(-source_type, -n) %>%
    group_by(kinase) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    filter(n > k_sampN) %>%
    select(-n)
  
  # if TRUE unsign the p-value
  if(unsignPval) ka <- ka %>% mutate(log10P = abs(log10P))
  
  # prepare data to model
  if(in_kinase){
    dataMod <- ka %>%
      inner_join(mutMat, by = c("kinase" = "gene", "sample")) %>%
      select(kinase, sample, log10P, mutated) %>%
      group_by(kinase) %>%
      nest() %>%
      ungroup()
  }else{
    dataMod <- ka %>%
      #filter(kinase == "CAMK2A") %>%
      inner_join(mutMat, by = c("sample")) %>%
      select(kinase, gene, sample, log10P, mutated) %>%
      group_by(kinase, gene) %>%
      nest() %>%
      ungroup()
  }
  
  # perform the linear modelling
  if(!multi_core){
    dataMod <- dataMod %>%
      mutate(model = map(.x = data, .f = linear_model, covs = covars, mutN = mut_kaN, NonMutN = NonMut_kaN))
  }else{
    cores = detectCores()-1
    cl <- makeCluster(cores)
    clusterEvalQ(cl, library(tidyverse))
    
    models <- parLapply(
      cl = cl,
      X = dataMod$data,
      fun = linear_model,
      covs = covars,
      mutN = mut_kaN,
      NonMutN = NonMut_kaN)
    
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



plotKinase <- function(mutMat, ka, k, g, kinase_list = "DB_text-mining", unsignPval = F, k_subN = 3){
  
  mutMat <- mutMat %>%
    select(-n) %>%
    filter(gene == g) %>%
    pivot_longer(cols = -c(gene), names_to = "sample", values_to = "mutated")
  
  if(unsignPval) ka <- ka %>% mutate(log10P = abs(log10P))
  
  ka <- ka %>%
    filter(source_type %in% kinase_list) %>%
    filter(kinase == k) %>%
    filter(n >= k_subN) %>%
    select(-source_type, -n)
  
  to_plot <- inner_join(ka, mutMat, by = "sample") %>%
    select(kinase, gene, everything()) %>%
    mutate(mutated = as.character(mutated))
  
  plot <- to_plot %>%
    ggplot(mapping = aes(x = mutated, y = log10P)) +
    geom_boxplot(notch = FALSE, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.2) +
    stat_compare_means() +
    theme_classic() +
    theme(
      axis.title = element_text(colour = "black", size = 12),
      axis.text = element_text(colour = "black", size = 10)) +
    labs(x = paste(g, "mutation status", sep = " "), y = paste(k, "log10P (activation score)", sep = " "))
  
  return(plot)
}
