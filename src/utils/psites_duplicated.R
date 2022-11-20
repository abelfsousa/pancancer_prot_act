library(tidyverse)


# function to check if the measurements across duplicated phosphosites are the same

psites_duplicated <- function(df, indx){
  dat <- as.data.frame(df)
  dat <- dat[, c(indx:ncol(dat))]
  
  for(i in 1:nrow(dat)) dat[i,] <- str_replace_na(dat[i,])
  
  same = c()
  for(i in 1:(nrow(dat)-1)) same <- c(same, sum(dat[1,] == dat[i+1,]) == ncol(dat))
  same <- sum(same) == (nrow(dat)-1)
  
  return(sum(same))
}
