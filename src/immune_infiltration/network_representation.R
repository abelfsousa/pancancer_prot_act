library(tidyverse)
library(igraph)
library(networkD3)


immune <- read_tsv("./output/files/kinase_tf_immune_cells_associations.txt") %>%
  filter(kin_cell_padj < 0.01)


cells <- unique(immune$cell)
kinases <- unique(immune$kinase)
tfs <- unique(immune$tf)

cell_kin <- immune %>%
  select(v1=cell, v2=kinase, weight=kin_cell_beta)

cell_tf <- immune %>%
  select(v1=cell, v2=tf, weight=tf_cell_beta)

kin_tf <- immune %>%
  select(v1=kinase, v2=tf, weight=kin_tf_beta)


g <- bind_rows(cell_kin, cell_tf, kin_tf) %>%
  graph_from_data_frame(directed = F)

v <- attr(V(g), "names")

V(g)$type <- if_else(v %in% cells, "cell", if_else(v %in% kinases, "kinase", "tf"))
V(g)$color <- if_else(v %in% cells, "red", if_else(v %in% kinases, "blue", "green"))
#V(g)$size <- degree(g)

edge_color <- function(edge){
  edge <- str_split_fixed(edge, "\\|", 2)
  v1 <- edge[, 1]
  v2 <- edge[, 2]
  if((v1 %in% cells & v2 %in% tfs) | (v2 %in% cells & v1 %in% tfs)){
    "grey"
  }
  else if((v1 %in% cells & v2 %in% kinases) | (v2 %in% cells & v1 %in% kinases)){
    "black"
  } else {
    "brown"
  }
}

E(g)$color <- map_chr(attr(E(g), "vnames"), edge_color)

plot(g)

