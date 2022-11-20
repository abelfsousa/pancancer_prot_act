library(tidyverse)


studies <- set_names(
  x = c("tcga-brca", "tcga-coread", "tcga-ov", "cbttc", "discovery-ccrcc", "discovery-luad", "discovery-ucec", "colon-opportunities", "ccle", "eogc-proteogenomics", "hcc-proteogenomics"),
  nm = c(rep("f1", 8), "f2", "f3", "f3"))


aa_code <- read_tsv(file = "./data/protein/protein_sequences/aminoacids_code.txt", skip = c(1)) %>%
  rbind(tibble(single_letter_code = "Ter", three_letter_code = "Ter",  name = "Termination codon"))


# tcga, cbttc, discovery and colon opportunities datasets
f1 <- function(df){
  
  if_indels <- c("In_Frame_Del", "In_Frame_Ins")
  fs_indels <- c("Frame_Shift_Del", "Frame_Shift_Ins")
  
  if (df$variant_class %in% if_indels){
    aa_pos = str_extract(df$HGVSp, "^p\\.[A-Za-z]+[0-9]+")
    pos = str_extract(aa_pos, "[0-9]+")
    aa = str_extract(aa_pos, "^p\\.[A-Za-z]+")
    aa1 = str_replace(aa, "p.", "")
    aa1 = as.character(aa_code[aa_code$three_letter_code == aa1, "single_letter_code"])
    res = tibble(prot_pos = pos, aa_wt = aa1, aa_mut = NA)
  }
  else if (df$variant_class %in% fs_indels){
    aa_pos = str_extract(df$HGVSp, "^p\\.[A-Za-z]+[0-9]+")
    pos = str_extract(aa_pos, "[0-9]+")
    aa = str_extract(aa_pos, "^p\\.[A-Za-z]+")
    aa1 = str_replace(aa, "p.", "")
    aa1 = as.character(aa_code[aa_code$three_letter_code == aa1, "single_letter_code"])
    aa2 = str_extract(df$HGVSp, "[0-9]+[A-Za-z]{3}")
    if(!is.na(aa2)){
      aa2 = str_extract(aa2, "[A-Za-z]+")
      aa2 = as.character(aa_code[aa_code$three_letter_code == aa2, "single_letter_code"])
    }
    res = tibble(prot_pos = pos, aa_wt = aa1, aa_mut = aa2)
  }
  else if(df$variant_class == "Missense_Mutation"){
    aa_pos = str_extract(df$HGVSp, "^p\\.[A-Za-z]+[0-9]+")
    pos = str_extract(aa_pos, "[0-9]+")
    aa = str_extract(aa_pos, "^p\\.[A-Za-z]+")
    aa1 = str_replace(aa, "p.", "")
    aa2 = str_extract(df$HGVSp, "[A-Za-z]+$")
    aa1 = as.character(aa_code[aa_code$three_letter_code == aa1, "single_letter_code"])
    aa2 = as.character(aa_code[aa_code$three_letter_code == aa2, "single_letter_code"])
    res = tibble(prot_pos = pos, aa_wt = aa1, aa_mut = aa2)
  }
  else if(df$variant_class == "Nonstop_Mutation"){
    stop_pos = str_extract(df$HGVSp, "^p\\.([A-Za-z]+|\\*{1})[0-9]+")
    stop_pos = str_extract(stop_pos, "[0-9]+")
    aa2 = str_extract(df$HGVSp, "[0-9]+[A-Za-z]{3}")
    aa2 = str_extract(aa2, "[A-Za-z]+")
    aa2 = as.character(aa_code[aa_code$three_letter_code == aa2, "single_letter_code"])
    res = tibble(prot_pos = stop_pos, aa_wt = "Ter", aa_mut = aa2)
  }
  else if(df$variant_class == "Nonsense_Mutation"){
    aa_pos = str_extract(df$HGVSp, "^p\\.[A-Za-z]+[0-9]+")
    pos = str_extract(aa_pos, "[0-9]+")
    aa = str_extract(aa_pos, "^p\\.[A-Za-z]+")
    aa1 = str_replace(aa, "p.", "")
    aa1 = as.character(aa_code[aa_code$three_letter_code == aa1, "single_letter_code"])
    res = tibble(prot_pos = pos, aa_wt = aa1, aa_mut = "Ter")
  }
  else if(df$variant_class == "Splice_Site"){
    if(is.na(df$HGVSp)){
      res = tibble(prot_pos = NA, aa_wt = NA, aa_mut = NA)
    }else{
      aa_pos = str_extract(df$HGVSp, "^p\\.[A-Za-z]+[0-9]+")
      pos = as.numeric(str_extract(aa_pos, "[0-9]+"))
      aa1 = str_extract(aa_pos, "^p.[A-Za-z]+")
      aa1 = str_replace(aa1, "p.", "")
      aa1 = as.character(aa_code[aa_code$three_letter_code == aa1, "single_letter_code"])
      equal = as.numeric(str_sub(df$seq, pos, pos) == aa1)
      pos = as.character(pos)
      if(str_detect(df$HGVSp, "del|dup")){
        aa2 = NA
      }else{
        aa2 = str_extract(df$HGVSp, "[A-Za-z]{3}[0-9]+[A-Za-z]{3}")
        aa2 = str_extract(aa2, "[A-Za-z]{3}$")
        aa2 = as.character(aa_code[aa_code$three_letter_code == aa2, "single_letter_code"])
      }
      res = tibble(prot_pos = pos, aa_wt = aa1, aa_mut = aa2, same = equal)
    }
  }
  return(res)
}


# ccle dataset
f2 <- function(df){

  fs_indels <- c("Frame_Shift_Del", "Frame_Shift_Ins")
  
  if(df$variant_class == "Missense_Mutation"){
    if(str_detect(df$HGVSp, "_", negate = T)){
      pos = str_extract(df$HGVSp, "[0-9]+")
      aa1 = str_extract(df$HGVSp, "^p\\.[A-Z]{1}")
      aa1 = str_replace(aa1, "p.", "")
      aa2 = str_extract(df$HGVSp, "[A-Z]{1}$")
      res = tibble(prot_pos = pos, aa_wt = aa1, aa_mut = aa2)
    }else{
      pos = str_extract(df$HGVSp, "[0-9]+_[0-9]+")
      pos = as.character(str_split_fixed(pos, "_", 2))[1]
      aas = str_extract(df$HGVSp, "[A-Z]+>[A-Z]+")
      aas = as.character(str_split(aas, ">", simplify = T))
      res = tibble(prot_pos = pos, aa_wt = aas[1], aa_mut = aas[2])
    }
  }
  else if(df$variant_class == "Nonsense_Mutation"){
    if(str_detect(df$HGVSp, "_", negate = T)){
      pos = str_extract(df$HGVSp, "[0-9]+")
      aa1 = str_extract(df$HGVSp, "^p\\.[A-Z]{1}")
      aa1 = str_replace(aa1, "p.", "")
      res = tibble(prot_pos = pos, aa_wt = aa1, aa_mut = "Ter")
    }else{
      pos = str_extract(df$HGVSp, "[0-9]+_[0-9]+")
      pos = as.character(str_split_fixed(pos, "_", 2))[1]
      aa1 = str_extract(df$HGVSp, "[A-Z]+>")
      aa1 = str_replace(aa1, ">", "")
      res = tibble(prot_pos = pos, aa_wt = aa1, aa_mut = "Ter")
    }
  }
  else if(df$variant_class %in% fs_indels){
    pos = str_extract(df$HGVSp, "[0-9]+")
    aa1 = str_extract(df$HGVSp, "^p\\.[A-Z]+")
    aa1 = str_replace(aa1, "p.", "")
    res = tibble(prot_pos = pos, aa_wt = aa1, aa_mut = NA)
  }
  else if(df$variant_class == "In_Frame_Ins"){
    pos = str_extract(df$HGVSp, "[0-9]+_[0-9]+")
    pos = as.character(str_split_fixed(pos, "_", 2))[1]
    aas = str_extract(df$HGVSp, "[A-Z]+>[A-Z]+")
    aas = as.character(str_split(aas, ">", simplify = T))
    res = tibble(prot_pos = pos, aa_wt = aas[1], aa_mut = aas[2])
  }
  else if(df$variant_class == "In_Frame_Del"){
    if(str_detect(df$HGVSp, "del")){
      pos = str_extract(df$HGVSp, "[0-9]+")
      aa1 = str_extract(df$HGVSp, "^p\\.[A-Z]+")
      aa1 = str_replace(aa1, "p.", "")
      res = tibble(prot_pos = pos, aa_wt = aa1, aa_mut = NA)
    }else{
      pos = str_extract(df$HGVSp, "[0-9]+_[0-9]+")
      pos = as.character(str_split_fixed(pos, "_", 2))[1]
      aas = str_extract(df$HGVSp, "[A-Z]+>[A-Z]+")
      aas = as.character(str_split(aas, ">", simplify = T))
      res = tibble(prot_pos = pos, aa_wt = aas[1], aa_mut = aas[2])
    }
  }
  else if(df$variant_class == "Nonstop_Mutation"){
    pos = str_extract(df$HGVSp, "[0-9]+")
    aa2 = str_extract(df$HGVSp, "[A-Z]{1}$")
    res = tibble(prot_pos = pos, aa_wt = "Ter", aa_mut = aa2)
  }
  else if(df$variant_class == "Splice_Site"){
    if(is.na(df$HGVSp)){
      res = tibble(prot_pos = NA, aa_wt = NA, aa_mut = NA)
    }
    else{
      pos = str_extract(df$HGVSp, "[0-9]+")
      aa1 = str_extract(df$HGVSp, "^p\\.[A-Z]{1}")
      aa1 = str_replace(aa1, "p.", "")
      aa2 = str_extract(df$HGVSp, "[A-Za-z*]+$")
      if(aa2 == "*") aa2 = "Ter"
      if(aa2 %in% c("fs", "del")) aa2 = NA
      res = tibble(prot_pos = pos, aa_wt = aa1, aa_mut = aa2)
    }
  }
  return(res)
}


# eogc and hcc datasets
f3 <- function(df){
  
  fs_indels <- c("Frame_Shift_Del", "Frame_Shift_Ins")
  
  if(df$variant_class == "Missense_Mutation"){
    pos = str_extract(df$HGVSp, "[0-9]+")
    aa1 = str_extract(df$HGVSp, "^p\\.[A-Z]{1}")
    aa1 = str_replace(aa1, "p.", "")
    aa2 = str_extract(df$HGVSp, "[A-Z]{1}$")
    res = tibble(prot_pos = pos, aa_wt = aa1, aa_mut = aa2)
  }
  else if(df$variant_class == "Nonsense_Mutation"){
    pos = str_extract(df$HGVSp, "[0-9]+")
    aa1 = str_extract(df$HGVSp, "^p\\.[A-Z]{1}")
    aa1 = str_replace(aa1, "p.", "")
    res = tibble(prot_pos = pos, aa_wt = aa1, aa_mut = "Ter")
  }
  else if(df$variant_class %in% fs_indels){
    pos = str_extract(df$HGVSp, "[0-9]+")
    aa1 = str_extract(df$HGVSp, "^p\\.[A-Z]{1}")
    aa1 = str_replace(aa1, "p.", "")
    res = tibble(prot_pos = pos, aa_wt = aa1, aa_mut = NA)
  }
  else if(df$variant_class == "In_Frame_Ins"){
    pos = str_extract(df$HGVSp, "[0-9]+")
    aa1 = str_extract(df$HGVSp, "^p\\.[A-Z]{1}")
    aa1 = str_replace(aa1, "p.", "")
    aa2 = str_extract(df$HGVSp, "[A-Z]+$")
    res = tibble(prot_pos = pos, aa_wt = aa1, aa_mut = aa2)
  }
  else if(df$variant_class == "Nonstop_Mutation"){
    pos = str_extract(df$HGVSp, "[0-9]+")
    aa2 = str_extract(df$HGVSp, "[A-Z]{1}$")
    res = tibble(prot_pos = pos, aa_wt = "Ter", aa_mut = aa2)
  }
  else if(df$variant_class == "Splice_Site"){
    res = tibble(prot_pos = NA, aa_wt = NA, aa_mut = NA)
  }
  return(res)
}


read_hgvsp <- function(df){
  
  selFun <- get(names(studies[studies == df$batch]))
  
  res <- selFun(df)
  
  return(res)
}



