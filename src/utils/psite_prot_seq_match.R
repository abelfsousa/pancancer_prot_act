library(tidyverse)


# function to check if in the protein sequences the aa on a given position is the phosphosite from the phosphorylation data

map_psite <- function(df){
  
  psites = df$psites
  psite = df$psite
  seq = df$seq
  
  if(psites == 1){
    aa = str_extract(psite, "[A-Z]{1}")
    aa_pos = str_extract(psite, "[0-9]+")
    
    same = as.numeric(str_sub(seq, aa_pos, aa_pos) == aa)
  }
  else{
    aas = str_extract_all(psite, "[A-Z]{1}")[[1]]
    ass_pos = str_extract_all(psite, "[0-9]+")[[1]]
    
    same = c()
    for (i in 1:length(aas)){
      same = c(same, as.numeric(str_sub(seq, ass_pos[i], ass_pos[i]) == aas[i]))
    }
    
    same = if(sum(same) == length(aas)) 1 else 0
    
  }
  
  return(same)
}


# function to check if in the protein sequences the aa on a given position is the phosphosite from kinase substrates lists

map_psite2 <- function(df){
  
  residue = df$residue
  position = df$position
  seq = df$seq
  
  same = as.numeric(str_sub(seq, position, position) == residue)

  return(same)
}
