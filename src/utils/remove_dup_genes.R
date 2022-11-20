# this function select the duplicated gene with the highest average gene expression value across samples

remove_dup_genes <- function(table_obj, gene_col1, gene_col2, sample_col, expr_col){
  
  gene_col1_i <- enquo(gene_col1)
  gene_col2_i <- enquo(gene_col2)
  sample_col_i <- enquo(sample_col)
  expr_col_i <- enquo(expr_col)
  expr_var <- rlang::as_name(expr_col_i)
  
  table_obj2 <- table_obj %>%
    nest(data = c(!! sample_col_i, !! expr_col_i)) %>%
    mutate(avg_expr = map_dbl(.x = data, .f = function(x) mean(as.data.frame(x)[, expr_var], na.rm = T))) %>%
    group_by(!! gene_col2_i) %>%
    top_n(n = 1, wt = avg_expr) %>%
    ungroup() %>%
    select(-avg_expr) %>%
    unnest(cols = data) %>%
    spread(key = rlang::as_name(sample_col_i), value = rlang::as_name(expr_col_i))
  
  return(table_obj2)
  
}
