# timma

usage:

1. loading R functions into environment__
   file.sources = list.files(pattern="*.R")__ 
   sapply(file.sources,source,.GlobalEnv)
2. loading Rcpp functions into environment__
   sourceCpp("max_min.cpp")