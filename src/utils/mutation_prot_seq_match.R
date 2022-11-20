library(tidyverse)


# functions to check if the amino-acid changed is the same on the protein sequenced



# function used with cbioportal (TCGA and CBTTC), discovery and colon datasets
cbio_datasets <- function(df, code){
  
  mutations <- c(
    "Missense_Mutation",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Nonsense_Mutation")
    
  if (df$variant_class %in% mutations){
    if (df$HGVSp == ""){
      same = 0
    } else{
    aa_pos = str_extract(df$HGVSp, "^p\\.[A-Za-z]+[0-9]+")
    aa_pos = str_extract(aa_pos, "[0-9]+")
    aa_original = str_extract(df$HGVSp, "^p.[A-Za-z]+")
    aa_original = str_replace(aa_original, "p.", "")
    aa_original = as.character(code[code$three_letter_code == aa_original, "single_letter_code"])
    same = as.numeric(str_sub(df$seq, aa_pos, aa_pos) == aa_original)
    }
  }
  else if(df$variant_class == "Nonstop_Mutation"){
    if (df$HGVSp == ""){
      same = 0
    } else{
    stop_pos = str_extract(df$HGVSp, "^p\\.([A-Za-z]+|\\*{1})[0-9]+")
    stop_pos = str_extract(stop_pos, "[0-9]+")
    same = as.numeric(nchar(df$seq)+1 == stop_pos)
    }
  }
  else if(df$variant_class == "Splice_Site"){
    same = 1
  }
  
  return(same)
}



# function for CCLE depmap dataset
ccle_dataset <- function(df){
  #print(as.data.frame(df))
  mutations <- c(
    "Missense_Mutation",
    "Nonsense_Mutation")
  
  mutations_fs <- c(
    "Frame_Shift_Del",
    "Frame_Shift_Ins")
  
  if (df$variant_class %in% mutations){
    if (df$HGVSp == ""){
      same = 0
    }
    else if(str_detect(df$HGVSp, "_")){
      aa_pos = str_extract(df$HGVSp, "[0-9]+_[0-9]+")
      aa_pos = as.character(str_split_fixed(aa_pos, "_", 2))
      aa_pos = as.numeric(aa_pos)
      aa_originals = str_extract(df$HGVSp, "([A-Z]*>)")
      aa_originals = str_replace(aa_originals, ">", "")
      same = as.numeric(str_sub(df$seq, aa_pos[1], aa_pos[2]) == aa_originals)
    }else{
      aa_pos = str_extract(df$HGVSp, "^p\\.[A-Z]{1}[0-9]+")
      aa_pos = str_extract(aa_pos, "[0-9]+")
      aa_original = str_extract(df$HGVSp, "^p\\.[A-Z]{1}")
      aa_original = str_replace(aa_original, "p.", "")
      same = as.numeric(str_sub(df$seq, aa_pos, aa_pos) == aa_original)
    }
  }
  else if (df$variant_class %in% mutations_fs){
    if (df$HGVSp == ""){
      same = 0
    } else{
      aa_pos = str_extract(df$HGVSp, "^p\\.[A-Z]+[0-9]+")
      aa_pos = str_extract(aa_pos, "[0-9]+")
      aa_pos = as.numeric(aa_pos)
      aa_originals = str_extract(df$HGVSp, "^p\\.[A-Z]+")
      aa_originals = str_replace(aa_originals, "p.", "")
      aa_originals_l = nchar(aa_originals)
      same = as.numeric(str_sub(df$seq, aa_pos, aa_pos+(aa_originals_l-1)) == aa_originals)
    }
  }
  else if(df$variant_class == "In_Frame_Ins"){
    if (df$HGVSp == "" | str_detect(df$HGVSp, ">", negate = T)){
      same = 0
    } else {
      aa_pos = str_extract(df$HGVSp, "^p\\.[0-9]+")
      aa_pos = str_extract(aa_pos, "[0-9]+")
      aa_original = str_extract(df$HGVSp, "[A-Z]*>")
      aa_original = str_replace(aa_original, ">", "")
      same = as.numeric(str_sub(df$seq, aa_pos, aa_pos) == str_sub(aa_original,1,1))
    }
  }
  else if(df$variant_class == "In_Frame_Del"){
    if (df$HGVSp == ""){
      same = 0
    } else if (str_detect(df$HGVSp, "del")){
      aa_pos = str_extract(df$HGVSp, "^p\\.[A-Z]+[0-9]+")
      aa_pos = str_extract(aa_pos, "[0-9]+")
      aa_pos = as.numeric(aa_pos)
      aa_originals = str_extract(df$HGVSp, "^p\\.[A-Z]+")
      aa_originals = str_replace(aa_originals, "p.", "")
      aa_originals_l = nchar(aa_originals)
      same = as.numeric(str_sub(df$seq, aa_pos, aa_pos+(aa_originals_l-1)) == aa_originals)
    }
    else{
      aa_pos = str_extract(df$HGVSp, "[0-9]+_[0-9]+")
      aa_pos = as.character(str_split_fixed(aa_pos, "_", 2))
      aa_pos = as.numeric(aa_pos)
      aa_originals = str_extract(df$HGVSp, "([A-Z]*>)")
      aa_originals = str_replace(aa_originals, ">", "")
      same = as.numeric(str_sub(df$seq, aa_pos[1], aa_pos[2]) == aa_originals)
    }
  }
  else if(df$variant_class == "Nonstop_Mutation"){
    if (df$HGVSp == ""){
      same = 0
    } else{
      stop_pos = str_extract(df$HGVSp, "^p\\.([A-Za-z]+|\\*{1})[0-9]+")
      stop_pos = str_extract(stop_pos, "[0-9]+")
      same = as.numeric(nchar(df$seq)+1 == stop_pos)
    }
  }
  else if(df$variant_class == "Splice_Site"){
    if (df$HGVSp == ""){
      same = 1
    } else{
      aa_pos = str_extract(df$HGVSp, "^p\\.[A-Z]+[0-9]+")
      aa_pos = str_extract(aa_pos, "[0-9]+")
      aa_pos = as.numeric(aa_pos)
      aa_originals = str_extract(df$HGVSp, "^p\\.[A-Z]+")
      aa_originals = str_replace(aa_originals, "p.", "")
      aa_originals_l = nchar(aa_originals)
      same = as.numeric(str_sub(df$seq, aa_pos, aa_pos+(aa_originals_l-1)) == aa_originals)
    }
  }
  
  return(same)
}




# function for GC and HCC datasets
gc_hcc_datasets <- function(df){
  
  mutations <- c(
    "Missense_Mutation",
    "Nonsense_Mutation")
  
  mutations_fs <- c(
    "Frame_Shift_Del",
    "Frame_Shift_Ins")
  
  if (df$variant_class %in% mutations){
    if (df$HGVSp == ""){
      same = 0
    }
    else{
      aa_pos = str_extract(df$HGVSp, "^p\\.[A-Z]{1}[0-9]+")
      aa_pos = str_extract(aa_pos, "[0-9]+")
      aa_original = str_extract(df$HGVSp, "^p\\.[A-Z]{1}")
      aa_original = str_replace(aa_original, "p.", "")
      same = as.numeric(str_sub(df$seq, aa_pos, aa_pos) == aa_original)
    }
  }
  else if (df$variant_class %in% mutations_fs){
    if (df$HGVSp == "" | str_detect(df$HGVSp, "del")){
      same = 0
    }else{
      aa_pos = str_extract(df$HGVSp, "^p\\.[A-Z]{1}[0-9]+")
      aa_pos = str_extract(aa_pos, "[0-9]+")
      aa_original = str_extract(df$HGVSp, "^p\\.[A-Z]{1}")
      aa_original = str_replace(aa_original, "p.", "")
      same = as.numeric(str_sub(df$seq, aa_pos, aa_pos) == aa_original)
    }
  }
  else if(df$variant_class == "In_Frame_Del"){
    same = 0
  }
  else if(df$variant_class == "In_Frame_Ins"){
    aa_pos = str_extract(df$HGVSp, "^p\\.[A-Z]{1}[0-9]+")
    aa_pos = str_extract(aa_pos, "[0-9]+")
    aa_original = str_extract(df$HGVSp, "^p\\.[A-Z]{1}")
    aa_original = str_replace(aa_original, "p.", "")
    same = as.numeric(str_sub(df$seq, aa_pos, aa_pos) == aa_original)
  }
  else if(df$variant_class == "Nonstop_Mutation"){
    if (df$HGVSp == ""){
      same = 0
    } else{
      stop_pos = str_extract(df$HGVSp, "^p\\.([A-Za-z]+|\\*{1})[0-9]+")
      stop_pos = str_extract(stop_pos, "[0-9]+")
      same = as.numeric(nchar(df$seq)+1 == stop_pos)
    }
  }
  else if(df$variant_class == "Splice_Site"){
    same = 1
  }
  
  return(same)
}
