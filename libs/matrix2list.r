
matrix2list <- function(x, dim = 2) lapply(apply(x, dim, list), unlist) 
