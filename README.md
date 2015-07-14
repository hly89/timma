# timma

usage:

1. loading R functions into environment
   file.sources = list.files(pattern="*.R")
   sapply(file.sources,source,.GlobalEnv)
2. loading Rcpp functions into environment
   sourceCpp("max_min.cpp")