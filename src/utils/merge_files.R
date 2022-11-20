
merge_files <- function(files_info, files_dir, gene_ids){
  
  # this function binds together several gene expression files per sample into a single table
  # use with GDC gene expression files
  
  expression_file <- gene_ids
  
  for(i in 1:nrow(files_info)){
    
    file_i <- paste0(files_dir, as.character(files_info[i,"File.Name"]))
    file_i <- data.table::fread( file = file_i, col.names = c("gene", str_replace_all(as.character(files_info[i,"Case.ID"]), "-", ".")) ) %>%
      as_tibble()
    
    expression_file <- inner_join(expression_file, file_i, by = "gene")
    
  }
  
  return(expression_file)
  
}


merge_files_gc <- function(files_info, files_dir, gene_ids){
  
  # this function binds together several gene expression files per sample into a single table
  # use with the GC GEO gene expression files
  
  expression_file <- gene_ids
  
  for(i in 1:nrow(files_info)){
    
    file_i <- paste0(files_dir, as.character(files_info[i,"file"]))
    file_i <- data.table::fread( file = file_i, skip = 1, col.names = c("gene_id", as.character(files_info[i,"geo_ID"])) ) %>%
      as_tibble()
    
    expression_file <- inner_join(expression_file, file_i, by = "gene_id")
    
  }
  
  return(expression_file)
  
}
