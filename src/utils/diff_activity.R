library(tidyverse)
library(parallel)



kinase_diff <- function(KA, covar, Ns, Ng){
  
  act <- KA %>%
    group_by_at(.vars = covar) %>%
    summarise(log10P = list(log10P), nSamp = n()) %>%
    ungroup() %>%
    filter(nSamp >= Ns)
  
  if(!(nrow(act) >= Ng)){
    nms <- c(covar, "N", "nSamp", "totalMedianAct", "medianAct", "coeff", "coeffPvalue", "anovaPvalue")
    res <- as_tibble(matrix(ncol = 8, nrow = 1)) %>%
      rename_all(.funs = ~ nms)
  } else {
    
    act <- unnest(act)
    
    fm <- as.formula(str_c("log10P", "~", covar, "+0", sep = ""))
    
    lreg <- lm(fm, data = act)
    av <- anova(lreg)
    
    var_levels <- levels(as.factor(act[[covar]]))
    dummies <- colnames(contrasts(as.factor(act[[covar]])))
    out_dummy <- var_levels[!(var_levels %in% dummies)]
    
    lreg <- lreg %>%
      broom::tidy() %>%
      select(term, coeff = estimate, coeffPvalue = p.value) %>%
      mutate(term = str_replace(term, covar, "")) %>%
      mutate(term = str_replace(term, "\\(Intercept\\)", out_dummy)) %>%
      rename_at(.vars = "term", .funs = ~ covar)
    
    av <- av %>%
      broom::tidy()
    
    res <- act %>%
      mutate(totalMedianAct = median(log10P)) %>%
      group_by_at(.vars = covar) %>%
      summarise(nSamp = unique(nSamp), log10P = list(log10P), totalMedianAct = unique(totalMedianAct)) %>%
      ungroup() %>%
      mutate(N = n()) %>%
      select(!! sym(covar), N, everything()) %>%
      mutate(medianAct = map_dbl(.x = log10P, .f = ~ median(.x, na.rm = T))) %>%
      #mutate(signal = sign(medianAct), actLog2fc = log2(abs(medianAct)/abs(totalMedianAct))*signal) %>%
      #select(-signal, -log10P) %>%
      select(-log10P) %>%
      inner_join(lreg, by = covar) %>%
      mutate(anovaPvalue = av[av$term == covar, "p.value", drop = T])
  }
  
  return(res)
}




differential_activity <- function(kinAct, metadata, cov, group_vars = NULL, kin_sub = 3, samples_group = 10, groups = 2, multi_core = FALSE){
  
  kin <- kinAct %>%
    filter(n >= kin_sub) %>%
    inner_join(metadata, by = "sample") %>%
    group_by_at(.vars = c("kinase", group_vars)) %>%
    nest(.key = "activities") %>%
    ungroup()
  
  if(!multi_core){
    kin <- kin %>%
      mutate(diff = map(.x = activities, .f = kinase_diff, covar = cov, Ns = samples_group, Ng = groups))
  }else{
    cores = detectCores()-1
    diffList <- mclapply(
      X = kin$activities,
      FUN =  kinase_diff,
      mc.cores = cores,
      covar = cov, Ns = samples_group, Ng = groups)
    
    kin <- kin %>%
      mutate(diff = diffList)
  }
  
  kin <- kin %>%
    select(-activities) %>%
    unnest() %>%
    filter_all(all_vars(!is.na(.))) %>%
    group_by_at(.vars = c("kinase", group_vars)) %>%
    nest() %>%
    ungroup() %>%
    mutate(anovaPvalue = map_dbl(.x = data, .f = ~ unique(.x$anovaPvalue))) %>%
    mutate(anovaPadj = p.adjust(anovaPvalue, "BH")) %>%
    select(-anovaPvalue) %>%
    unnest() %>%
    select(-anovaPadj, anovaPadj)
  
  return(kin)
}
