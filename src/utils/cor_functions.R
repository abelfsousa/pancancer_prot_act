library(tidyverse)


ecor_test <- function(x, y, method = "pearson"){
  if(length(x) != length(y)) stop("x and y must have the same length")
  
  cor <- tryCatch(
    
    broom::tidy(cor.test(x, y, method=method)),
    
    error=function(cond){
      x <- broom::tidy(cor.test(rnorm(10), rnorm(10), method=method))
      x <- x %>%
        mutate_if(is.numeric, function(x) x=0) %>%
        mutate_if(is.character, function(x) x="NA")
      return(x)})
  
  cor <- cor %>%
    rename_all(.funs = ~ str_c(method, .x, sep = "_"))
  
  return(cor)
}

cor_test <- function(x, y, method = "pearson"){
  if(length(x) != length(y)) stop("x and y must have the same length")
  
  cor <- broom::tidy(cor.test(x, y, method=method))
  
  cor <- cor %>%
    rename_all(.funs = ~ str_c(method, .x, sep = "_"))
  
  return(cor)
}

corr <- function(df, x, y, method = "pearson"){
  a = df[[x]]
  b = df[[y]]
  cor <- broom::tidy(cor.test(a, b, method = method))
  
  res <- cor %>%
    select(estimate, p.value)
  
  return(res)
}

