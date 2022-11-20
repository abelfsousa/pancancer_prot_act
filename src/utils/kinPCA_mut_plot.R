kinPCA_mut_plot <- function(pcs, loadings, mutations, metadata, kinMut, x, y, gx, gy, lim){
  
  # first scatterplot
  df1 <- pcs %>%
    filter(pc %in% c(x, y)) %>%
    pivot_wider(names_from = "pc", values_from = "value") %>%
    inner_join(mutations[mutations$gene == gx, ], by = "sample") %>%
    mutate(mutated = as.character(mutated)) %>%
    inner_join(metadata[, c("sample", "batch", "cancer")], by = "sample")
  
  scatterplot1 <- df1 %>%
    ggplot() +
    theme_classic() +
    geom_point(mapping = aes(x = get(x), y = get(y), color = mutated), alpha = 0.7) +
    geom_hline(yintercept = 0, lwd = .5, color = "grey") +
    geom_vline(xintercept = 0, lwd = .5, color = "grey") +
    geom_segment(aes(x = 10, y = -.2, xend = 10, yend = .2), lwd = .5, color = "grey") +
    geom_segment(aes(x = 5, y = -.2, xend = 5, yend = .2), lwd = .5, color = "grey") +
    geom_segment(aes(x = -10, y = -.2, xend = -10, yend = .2), lwd = .5, color = "grey") +
    geom_segment(aes(x = -5, y = -.2, xend = -5, yend = .2), lwd = .5, color = "grey") +
    geom_segment(aes(x = -.2, y = 5, xend = .2, yend = 5), lwd = .5, color = "grey") +
    geom_segment(aes(x = -.2, y = 10, xend = .2, yend = 10), lwd = .5, color = "grey") +
    geom_segment(aes(x = -.2, y = -5, xend = .2, yend = -5), lwd = .5, color = "grey") +
    geom_segment(aes(x = -.2, y = -10, xend = .2, yend = -10), lwd = .5, color = "grey") +
    scale_y_continuous(limits = c(-lim, lim)) +
    scale_x_continuous(limits = c(-lim, lim)) +
    scale_color_manual(values = c("#67a9cf", "red"), labels = c("No", "Yes")) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      plot.title = element_text(hjust = 0.5)) +
    labs(x = x, y = y, title = str_c(gx, " (", x, " associated)"),  color = str_c(gx,  " LoF mutation"))
  
  mutGx <- kinMut[kinMut$gene == gx, "kinase", drop=T]
  
  if(length(mutGx) > 0){
    df1_loadings <- loadings %>%
      filter(pc %in% c(x, y)) %>%
      pivot_wider(names_from = "pc", values_from = "value") %>%
      rename(kinase = gene) %>%
      filter(kinase %in% mutGx) %>%
      mutate(exp = map2(.x = get(x), .y = get(y), .f = expand_coord, limit = lim)) %>%
      unnest()
    
    scatterplot1 <- scatterplot1 +
      geom_segment(data = df1_loadings, aes(x = 0, y = 0, xend = x, yend = y), arrow = arrow(length = unit(1/2, "picas"))) +
      geom_text_repel(data = df1_loadings, mapping = aes(x = x, y = y, label = kinase), box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines"), size = 4, colour="black")
  }
  
  barplot1 <- df1 %>%
    mutate(cancer = fct_infreq(cancer)) %>%
    ggplot() +
    theme_classic() +
    coord_flip() +
    geom_bar(mapping = aes(x = cancer, fill = mutated), position = "fill") +
    scale_fill_manual(values = c("#67a9cf", "red"), labels = c("No", "Yes"), guide = F) +
    theme(axis.text = element_text(color = "black"), axis.title = element_text(color = "black")) +
    labs(x = "Cancer", y = "%", fill = "Mutated")
  
  plot1 <- ggdraw(scatterplot1) +
    draw_plot(barplot1, .66, .65, .33, .3)
  
  
  
  # second scatterplot
  df2 <- pcs %>%
    filter(pc %in% c(x, y)) %>%
    pivot_wider(names_from = "pc", values_from = "value") %>%
    inner_join(mutations[mutations$gene == gy, ], by = "sample") %>%
    mutate(mutated = as.character(mutated)) %>%
    inner_join(metadata[, c("sample", "batch", "cancer")], by = "sample")
  
  scatterplot2 <- df2 %>%
    ggplot() +
    theme_classic() +
    geom_point(mapping = aes(x = get(x), y = get(y), color = mutated), alpha = 0.7) +
    geom_hline(yintercept = 0, lwd = .5, color = "grey") +
    geom_vline(xintercept = 0, lwd = .5, color = "grey") +
    geom_segment(aes(x = 10, y = -.2, xend = 10, yend = .2), lwd = .5, color = "grey") +
    geom_segment(aes(x = 5, y = -.2, xend = 5, yend = .2), lwd = .5, color = "grey") +
    geom_segment(aes(x = -10, y = -.2, xend = -10, yend = .2), lwd = .5, color = "grey") +
    geom_segment(aes(x = -5, y = -.2, xend = -5, yend = .2), lwd = .5, color = "grey") +
    geom_segment(aes(x = -.2, y = 5, xend = .2, yend = 5), lwd = .5, color = "grey") +
    geom_segment(aes(x = -.2, y = 10, xend = .2, yend = 10), lwd = .5, color = "grey") +
    geom_segment(aes(x = -.2, y = -5, xend = .2, yend = -5), lwd = .5, color = "grey") +
    geom_segment(aes(x = -.2, y = -10, xend = .2, yend = -10), lwd = .5, color = "grey") +
    scale_y_continuous(limits = c(-lim, lim)) +
    scale_x_continuous(limits = c(-lim, lim)) +
    scale_color_manual(values = c("#67a9cf", "red"), labels = c("No", "Yes")) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_blank(),
      plot.title = element_text(hjust = 0.5)) +
    labs(x = x, y = y, title = str_c(gy, " (", y, " associated)"),  color = str_c(gy,  " LoF mutation"))
  
  mutGy <- kinMut[kinMut$gene == gy, "kinase", drop=T]
  
  if(length(mutGy) > 0){
    df2_loadings <- loadings %>%
      filter(pc %in% c(x, y)) %>%
      pivot_wider(names_from = "pc", values_from = "value") %>%
      rename(kinase = gene) %>%
      filter(kinase %in% mutGy) %>%
      mutate(exp = map2(.x = get(x), .y = get(y), .f = expand_coord, limit = lim)) %>%
      unnest()
    
    scatterplot2 <- scatterplot2 +
      geom_segment(data = df2_loadings, aes(x = 0, y = 0, xend = x, yend = y), arrow = arrow(length = unit(1/2, "picas"))) +
      geom_text_repel(data = df2_loadings, mapping = aes(x = x, y = y, label = kinase), box.padding = unit(0.3, "lines"), point.padding = unit(0.3, "lines"), size = 4, colour="black")
  }
  
  barplot2 <- df2 %>%
    mutate(cancer = fct_infreq(cancer)) %>%
    ggplot() +
    theme_classic() +
    coord_flip() +
    geom_bar(mapping = aes(x = cancer, fill = mutated), position = "fill") +
    scale_fill_manual(values = c("#67a9cf", "red"), labels = c("No", "Yes"), guide = F) +
    theme(axis.text = element_text(color = "black"), axis.title = element_text(color = "black")) +
    labs(x = "Cancer", y = "%", fill = "Mutated")
  
  plot2 <- ggdraw(scatterplot2) +
    draw_plot(barplot2, .66, .65, .33, .3)
  
  scatterplots <- plot_grid(plot1, plot2)
  
  return(scatterplots)
}


expand_coord <- function(x, y, limit, factors = seq(10, 1000, 10)){
  limit <- limit*0.5
  xx <- x
  yy <- y
  i <- 1
  while(abs(xx) < limit & abs(yy) < limit){
    if(abs(x * factors[i]) > limit | abs(y  * factors[i]) > limit){
      break
    }else{
      xx = x * factors[i]
      yy = y  * factors[i]
      
      i = i + 1
    }
  }
  res <- tibble(x = xx, y = yy)
  return(res)
}


