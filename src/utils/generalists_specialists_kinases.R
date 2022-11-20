library(tidyverse)


kinReg <- function(kinAct, regCutoff, nSamples, token){
  
  kin <- kinAct %>%
    mutate(N = nSamples) %>%
    group_by(kinase, N) %>%
    summarise(n = sum(abs(log10P) > regCutoff)) %>%
    ungroup() %>%
    mutate(p = n/N) %>%
    rename_if(.predicate = is.numeric, .funs = ~ str_c(.x, token, sep = ""))
  
  kin
}


kinRegRand <- function(kinAct, NrandSamples, regCutoff, token){
  
  kinAct <- kinAct %>%
    group_by(sample) %>%
    nest() %>%
    ungroup() %>%
    sample_n(NrandSamples) %>%
    unnest()
  
  nSamples = kinAct$sample %>% unique() %>% length()
    
  kinReg <- kinReg(kinAct, regCutoff, nSamples, token)
  
  kinReg
}


gs <- function(kinActCancer, kinActPert, cutoffReg, topReg, z, permt, Npermt){
  
  # total number of cancer samples/perturbations
  tumours <- length(unique(kinActCancer$sample))
  pertbs <- length(unique(kinActPert$sample))
  
  if(tumours > pertbs & permt) permt = FALSE
  
  # count number of cancer samples in which a kinase is regulated
  cancer <- kinReg(kinActCancer, cutoffReg, tumours, "T")
  
  # count number of perturbations in which a kinase is regulated
  # permute over the perturbations randomly if specified
  if(permt){
    pertb <- tibble(perm = 1:Npermt, data = list(kinActPert)) %>%
      mutate(reg = map(.x = data, .f = kinRegRand, NrandSamples = tumours, regCutoff = cutoffReg, token = "P")) %>%
      select(-data) %>%
      unnest() %>%
      group_by(kinase) %>%
      summarise(NP = mean(NP), nP = mean(nP), pP = mean(pP)) %>%
      ungroup()
  }else{
    pertb <- kinReg(kinActPert, cutoffReg, pertbs, "P")
  }
  
  
  # join cancer/perturbation data
  KAs <- cancer %>%
    inner_join(pertb, by = "kinase") %>%
    filter(!(nP == 0 & nT == 0))
    #filter(nT > 0 & nP > 0)
  
  
  # top regulated kinases
  cancer_top <- top_n(KAs, topReg, nT) %>%
    mutate(topKin = "cancer")
  
  pertb_top <- top_n(KAs, topReg, nP) %>%
    mutate(topKin = "perturbation")
  
  topKin <- cancer_top %>%
    bind_rows(pertb_top) %>%
    group_by(kinase) %>%
    nest() %>%
    ungroup() %>%
    mutate(data = map(.x = data, .f = ~ if(nrow(.x) > 1){.x %>% mutate(topKin = "both") %>% distinct()}else{.x})) %>%
    unnest()
  
  
  # add regulated token to all kinases
  KAs <- KAs %>%
    left_join(topKin, by = c("kinase", "nP", "nT", "pP", "pT", "NP", "NT")) %>%
    mutate(topKin = replace_na(topKin, "others")) %>%
    mutate(topKin = fct_relevel(topKin, "others", "cancer", "perturbation", "both"))
  
  
  # estimate the cancer/perturbation specific kinases
  # use the linear regression residuals
  #lreg <- lm(nT ~ nP, KAs)
  lreg <- lm(pT ~ pP, KAs)
  
  KAs <- KAs %>%
    mutate(resd = residuals(lreg)) %>%
    mutate(zscore = scale(resd)[,1]) %>%
    mutate(class = if_else(zscore > z, "tumour", if_else(zscore < -z, "perturbation", "both")))
  
  
  return(KAs)
}




gs_mapper <- function(kinActCancer, kinActPert, cutoffNSubs = 3, cutoffKinReg = 1.75, topReg = 20, zeta = 2, group = FALSE, metadataCancer = NULL, group_var = NULL, perm = F, perm_n = NULL){
  
  # filter kinase measurements by number of substrates
  kinActCancer <- kinActCancer %>%
    filter(n >= cutoffNSubs) %>%
    select(-n)
  
  kinActPert <- kinActPert %>%
    filter(n >= cutoffNSubs) %>%
    select(-n)
  
  
  # perform the calculation by groups if specified
  if(group){
    variable <- enquo(group_var)
    
    metadataCancer <- metadataCancer %>%
      select(sample, !! variable)
    
    kinActCancer <- kinActCancer %>%
      inner_join(metadataCancer, by = "sample") %>%
      select(!! variable, everything()) %>%
      group_by(!! variable) %>%
      nest() %>%
      ungroup()
  }else{
    kinActCancer <- kinActCancer %>%
      nest()
  }
  
  reg <- kinActCancer %>%
    mutate(reg = map(.x = data, .f = gs, kinActPert = kinActPert, cutoffReg = cutoffKinReg, topReg = topReg, z = zeta, permt = perm, Npermt = perm_n)) %>%
    select(-data) %>%
    unnest()
  
  return(reg)
}

