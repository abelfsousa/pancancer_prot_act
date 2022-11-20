library(tidyverse)


# function to regress-out a set of covariates from a variable fitting a linear model
# the first argument must be a data.frame with at least 3 columns:
# the first with the sample identifiers, the second the response variable and the third the covariate to regress-out
# df can have more covariates to regress-out
# other covariates can be passed using the covs argument
# covs is a data.frame with one column "sample" containing the same sample identifiers as df
regress_outCovs <- function(df, covs = NULL, ztransf = FALSE){
  
  if(ncol(df) < 3) stop("Few columns. You need at least 3 columns in the data frame.")
  
  # additional covariates to regress-out
  if(!is.null(covs)){
    df <- df %>%
      inner_join(covs, by = "sample")
  }
  
  # remove NAs from the dataset
  dfM <- na.exclude(df)
  
  # if the number of observations is equal to or less than 1 return residuals as NA values
  if(nrow(dfM) <= 1){
    res = tibble(residual = rep(NA, nrow(df)))
  }else{
    # remove sample column
    samples <- dfM$sample
    dfM <- dfM %>%
      select(-sample)
    
    # define response variable
    resp1 <- colnames(dfM)[1]
    
    # remove factors with less then 2 levels
    # remove numeric variables with null variance
    dfM <- dfM %>%
      mutate_if(.predicate = is.character, .funs = as.factor) %>%
      mutate_if(.predicate = is.factor, .funs = fct_drop) %>%
      select_if(.predicate = ~ if(!is.factor(.x)){TRUE}else{if(nlevels(.x) == 1){FALSE}else{TRUE}}) %>%
      select_if(.predicate = ~ if(!is.numeric(.x)){TRUE}else{if(var(.x, na.rm=T) > 0){TRUE}else{FALSE}})
    
    # if the number of variables is less than or equal to 1 or the response variable is not in the data frame:
    # return residuals as NA values
    if(ncol(dfM) <= 1 | !(resp1 %in% colnames(dfM))){
      res = tibble(residual = rep(NA, nrow(df)))
    }else{
      # define response and explanatory variables
      resp2 <- colnames(dfM)[1]
      expl <- colnames(dfM)[-1]
      if(resp1 != resp2) stop("Response variables are different!")
      
      # set up the formula and define model matrix
      f <- as.formula(paste0(resp2, "~", paste(expl, collapse = "+")))
      model_mat <- model.matrix(object=f, data=dfM)
      
      # if the number of variables in the model is equal to or higher than the number of observations:
      # return residuals as NA values in these cases
      if(ncol(model_mat) >= nrow(model_mat)){
        res = tibble(residual = rep(NA, nrow(df)))
      }else{
        # z-score transform numeric variables if specified
        # re-calculate model matrix if so
        if(ztransf){
          dfM <- dfM %>%
            mutate_if(.predicate = is.numeric, .funs = ~ scale(.x)[,1])
          model_mat <- model.matrix(object=f, data=dfM)
        }
        
        # fit the model and get the residuals
        reg <- lm.fit(x = model_mat, y = dfM[[resp2]])
        reg_residual <- reg$residuals
        
        res = tibble(sample = samples, residual = reg_residual)
        res <- res %>%
          right_join(df, by = "sample") %>%
          #mutate(residual2 = map2_dbl(.x=residual, .y=y, .f=~if(is.na(.x) & !is.na(.y)){.y}else{.x})) %>%
          #select(-residual, residual=residual2) %>%
          select(residual)
      }
    }
  }
  return(res)
}


# simulated data
# sim1 <- tibble(
#   sample = paste0("s",1:14),
#   y = c(rnorm(12),rep(NA,2)),
#   x = c(rnorm(5),NA,rnorm(5),NA,rep(NA,2)))
# 
# sim2 <- tibble(
#   sample = paste0("s",1:14),
#   y = c(1:12,rep(NA,2)),
#   x = c(1:5,NA,7:11,NA,rep(NA,2)))
# 
# sim3 <- tibble(
#   sample = paste0("s",1:14),
#   y = c(rep(NA, 13), rnorm(1)),
#   x = c(rnorm(14)))
# 
# sim4 <- tibble(
#   sample = paste0("s",1:14),
#   y = rnorm(14),
#   x = rnorm(14))
# 
# sim5 <- tibble(
#   sample = paste0("s",1:14),
#   age = sample(x=20:80, 14, replace = T),
#   gender = sample(x=c("f", "m"), 14, replace = T),
#   batch1 = rep("a",14))
# 
# sim_list <- list(sim1, sim2, sim3, sim4, sim5)
# 
# expr <- sim_list[[1]]
# meta <- sim_list[[5]]
# 
# res <- regress_outCovs(expr, meta, ztransf = T)
# 
# 
# expr <- sim_list[[4]]
# expr1 <- expr %>% select(-x) %>% mutate(y = scale(y)[,1]) %>% pivot_wider(names_from = "sample", values_from = "y")
# expr1 <- as.matrix(expr1)
# 
# expr2 <- expr %>% select(x) %>% mutate(x = scale(x)[,1])
# expr2 <- as.matrix(expr2)
# 
# res1 <- limma::removeBatchEffect(x=expr1, covariates=expr2)
# res2 <- regress_outCovs(df = expr, ztransf = T)
# cor(res1[1,], res2$residual)



