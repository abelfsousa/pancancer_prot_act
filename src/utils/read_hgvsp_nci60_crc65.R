read_hgvsp <- function(mutation_type, hgvsp){
  
  fun <- get(mutation_type)
  
  res <- fun(hgvsp)
  
  return(res)
}


# -- functions for the NCI60 dataset

# missense mutations
Missense <- function(hgvsp){
  pos <- as.numeric(str_extract(hgvsp, "[0-9]+"))
  aa_wt <- str_extract(hgvsp, "^[A-Z]{1}")
  aa_mut <- str_extract(hgvsp, "[A-Z]{1}$")
  res = tibble(prot_pos_i = pos, prot_pos_f = prot_pos_i, aa_wt = aa_wt, aa_mut = aa_mut)
  res
}

# nonsense mutations
Nonsense <- function(hgvsp){
  if(str_detect(hgvsp, "_", negate = T)){
    pos <- as.numeric(str_extract(hgvsp, "[0-9]+"))
    aa_wt <- str_extract(hgvsp, "^[A-Z]{1}")
    aa_mut <- str_extract(hgvsp, "[A-Z]+$")
    res = tibble(prot_pos_i = pos, prot_pos_f = prot_pos_i, aa_wt = aa_wt, aa_mut = aa_mut)
  } else {
    pos <- str_extract(hgvsp, "^[A-Z]+[0-9]+_[A-Z]+[0-9]+")
    pos <- as.numeric(str_extract_all(pos, "[0-9]+")[[1]])
    aa_wt <- str_extract(hgvsp, "^[A-Z]+[0-9]+_[A-Z]+[0-9]+")
    aa_wt <- str_extract_all(aa_wt, "[A-Z]+")[[1]]
    aa_wt <- str_c(aa_wt, collapse = "")
    aa_mut <- str_extract(hgvsp, "[A-Z]+$")
    res = tibble(prot_pos_i = pos[1], prot_pos_f = pos[2], aa_wt = aa_wt, aa_mut = aa_mut)
  }
  res
}

# initiation loss mutations
Initiation_loss <- function(hgvsp){
  pos <- as.numeric(str_extract(hgvsp, "[0-9]{1}"))
  aa_wt <- str_extract(hgvsp, "^[A-Z]{1}")
  aa_mut <- str_extract(hgvsp, "[A-Z]+$")
  res = tibble(prot_pos_i = pos, prot_pos_f = prot_pos_i, aa_wt = aa_wt, aa_mut = aa_mut)
  res
}

# read through mutations
Read_through <- function(hgvsp){
  pos <- as.numeric(str_extract(hgvsp, "[0-9]+"))
  aa_wt <- str_extract(hgvsp, "^[A-Z]+")
  aa_mut <- str_extract(hgvsp, "[A-Z]+$")
  res = tibble(prot_pos_i = pos, prot_pos_f = prot_pos_i, aa_wt = aa_wt, aa_mut = aa_mut)
  res
}

# frameshift mutations
Frameshift <- function(hgvsp){
  if(str_detect(hgvsp, "_", negate = T)){
    pos <- as.numeric(str_extract(hgvsp, "[0-9]+"))
    aa_wt <- str_extract(hgvsp, "^[A-Z]+")
    aa_mut <- str_extract(hgvsp, "[A-Z]+$")
    res = tibble(prot_pos_i = pos, prot_pos_f = prot_pos_i, aa_wt = aa_wt, aa_mut = aa_mut)
  } else {
    pos <- as.numeric(str_extract_all(hgvsp, "[0-9]+")[[1]])
    aa_wt <- str_extract(hgvsp, "^[A-Z]+")
    aa_mut <- str_extract(hgvsp, "[A-Z]+$")
    res = tibble(prot_pos_i = pos[1], prot_pos_f = pos[2], aa_wt = aa_wt, aa_mut = aa_mut)
  }
  res
}

# nonframeshift mutations
Nonframeshift <- function(hgvsp){
  if(str_detect(hgvsp, "_", negate = T)){
    pos <- as.numeric(str_extract(hgvsp, "[0-9]+"))
    aa_wt <- str_extract(hgvsp, "^[A-Z]+")
    aa_mut <- str_extract(hgvsp, "[A-Z]+$")
    res = tibble(prot_pos_i = pos, prot_pos_f = prot_pos_i, aa_wt = aa_wt, aa_mut = aa_mut)
  } else {
    pos <- as.numeric(str_extract_all(hgvsp, "[0-9]+")[[1]])
    aa_wt <- str_extract(hgvsp, "^[A-Z]+")
    aa_mut <- str_extract(hgvsp, "[A-Z]+$")
    res = tibble(prot_pos_i = pos[1], prot_pos_f = pos[2], aa_wt = aa_wt, aa_mut = aa_mut)
  }
  res
}


# -- functions for the CRC65 dataset

# missense mutations
Missense_Mutation <- function(hgvsp){
  if(str_detect(hgvsp, "_", negate = T)){
    pos <- as.numeric(str_extract(hgvsp, "[0-9]+"))
    aa_wt <- str_extract(hgvsp, "^[A-Z]+")
    aa_mut <- str_extract(hgvsp, "[A-Z]+$")
    res = tibble(prot_pos_i = pos, prot_pos_f = prot_pos_i, aa_wt = aa_wt, aa_mut = aa_mut)
  } else {
    pos <- str_extract(hgvsp, "^[0-9]+_[0-9]+")
    pos <- as.numeric(str_extract_all(pos, "[0-9]+")[[1]])
    aa <- str_extract_all(hgvsp, "[A-Z]+")[[1]]
    res = tibble(prot_pos_i = pos[1], prot_pos_f = pos[2], aa_wt = aa[1], aa_mut = aa[2])
  }
  res
}

# nonsense mutations
Nonsense_Mutation <- function(hgvsp){
  if(str_detect(hgvsp, "_", negate = T)){
    pos <- as.numeric(str_extract(hgvsp, "[0-9]+"))
    aa_wt <- str_extract(hgvsp, "^[A-Z]+")
    aa_mut <- str_extract(hgvsp, "[A-Z*]+$")
    aa_mut <- str_replace(aa_mut, "\\*", "X")
    res = tibble(prot_pos_i = pos, prot_pos_f = prot_pos_i, aa_wt = aa_wt, aa_mut = aa_mut)
  } else {
    pos <- str_extract(hgvsp, "^[0-9]+_[0-9]+")
    pos <- as.numeric(str_extract_all(pos, "[0-9]+")[[1]])
    aa <- str_extract_all(hgvsp, "[A-Z*]+")[[1]]
    if(length(aa) == 1 & aa[1] == "*") aa <- c(NA, "*")
    aa <- str_replace(aa, "\\*", "X")
    res = tibble(prot_pos_i = pos[1], prot_pos_f = pos[2], aa_wt = aa[1], aa_mut = aa[2])
  }
  res
}

# initiation loss mutations
Start_Codon_Del <- function(hgvsp){
  res = tibble(prot_pos_i = NA, prot_pos_f = NA, aa_wt = NA, aa_mut = NA)
  res
}

Start_Codon_Ins <- function(hgvsp){
  res = tibble(prot_pos_i = NA, prot_pos_f = NA, aa_wt = NA, aa_mut = NA)
  res
}

Start_Codon_SNP <- function(hgvsp){
  pos <- as.numeric(str_extract(hgvsp, "[0-9]{1}"))
  aa_wt <- str_extract(hgvsp, "^[A-Z]{1}")
  aa_mut <- str_extract(hgvsp, "[A-Z]{1}$")
  res = tibble(prot_pos_i = pos, prot_pos_f = prot_pos_i, aa_wt = aa_wt, aa_mut = aa_mut)
  res
}

# read through mutations
Nonstop_Mutation <- function(hgvsp){
  pos <- as.numeric(str_extract(hgvsp, "[0-9]+"))
  aa_wt <- str_extract(hgvsp, "^[A-Z*]+")
  aa_wt <- str_replace(aa_wt, "\\*", "X")
  aa_mut <- str_extract(hgvsp, "[A-Z]+$")
  res = tibble(prot_pos_i = pos, prot_pos_f = prot_pos_i, aa_wt = aa_wt, aa_mut = aa_mut)
  res
}

Stop_Codon_Del <- function(hgvsp){
  res = tibble(prot_pos_i = NA, prot_pos_f = NA, aa_wt = NA, aa_mut = NA)
  res
}

Stop_Codon_Ins <- function(hgvsp){
  res = tibble(prot_pos_i = NA, prot_pos_f = NA, aa_wt = NA, aa_mut = NA)
  res
}


# frameshift mutations
Frame_Shift_Del <- function(hgvsp){
  pos <- as.numeric(str_extract(hgvsp, "[0-9]+"))
  aa_wt <- str_extract(hgvsp, "^[A-Z]+")
  aa_mut <- str_extract(hgvsp, "[A-Z]+$")
  res = tibble(prot_pos_i = pos, prot_pos_f = prot_pos_i, aa_wt = aa_wt, aa_mut = aa_mut)
  res
}

Frame_Shift_Ins <- function(hgvsp){
  pos <- as.numeric(str_extract(hgvsp, "[0-9]+"))
  aa_wt <- str_extract(hgvsp, "^[A-Z]+")
  aa_mut <- str_extract(hgvsp, "[A-Z]+$")
  res = tibble(prot_pos_i = pos, prot_pos_f = prot_pos_i, aa_wt = aa_wt, aa_mut = aa_mut)
  res
}

# nonframeshift mutations
In_Frame_Del <- function(hgvsp){
  if(str_detect(hgvsp, "_", negate = T)){
    pos <- as.numeric(str_extract(hgvsp, "[0-9]+"))
    aa_wt <- str_extract(hgvsp, "^[A-Z]+")
    aa_mut <- str_extract(hgvsp, "[A-Z]+$")
    res = tibble(prot_pos_i = pos, prot_pos_f = prot_pos_i, aa_wt = aa_wt, aa_mut = aa_mut)
  } else {
    pos <- str_extract(hgvsp, "^[0-9]+_[0-9]+")
    pos <- as.numeric(str_extract_all(pos, "[0-9]+")[[1]])
    aa <- str_extract_all(hgvsp, "[A-Z]+")[[1]]
    res = tibble(prot_pos_i = pos[1], prot_pos_f = pos[2], aa_wt = aa[1], aa_mut = aa[2])
  }
  res
}

In_Frame_Ins <- function(hgvsp){
  if(str_detect(hgvsp, "_", negate = T)){
    pos <- as.numeric(str_extract(hgvsp, "[0-9]+"))
    aa_wt <- str_extract(hgvsp, "^[A-Z]+")
    aa_mut <- str_extract(hgvsp, "[A-Z]+$")
    res = tibble(prot_pos_i = pos, prot_pos_f = prot_pos_i, aa_wt = aa_wt, aa_mut = aa_mut)
  } else {
    pos <- str_extract(hgvsp, "^[0-9]+_[0-9]+")
    pos <- as.numeric(str_extract_all(pos, "[0-9]+")[[1]])
    aa <- str_extract_all(hgvsp, "[A-Z*]+")[[1]]
    if(length(aa) == 1) aa <- c(NA, aa)
    aa <- str_replace(aa, "\\*", "X")
    res = tibble(prot_pos_i = pos[1], prot_pos_f = pos[2], aa_wt = aa[1], aa_mut = aa[2])
  }
  res
}


# -- functions for the NCI60 and CRC65 datasets
Silent <- function(hgvsp){
  pos <- as.numeric(str_extract(hgvsp, "[0-9]+"))
  aa_wt <- str_extract(hgvsp, "^[A-Z]+")
  aa_mut <- str_extract(hgvsp, "[A-Z]+$")
  res = tibble(prot_pos_i = pos, prot_pos_f = prot_pos_i, aa_wt = aa_wt, aa_mut = aa_mut)
  res
}
