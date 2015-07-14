# timma

usage:
121
1. loading R functions into environment <br />
   file.sources = list.files(pattern="*.R") <br />
   sapply(file.sources,source,.GlobalEnv)
2. loading Rcpp functions into environment <br />
   library(Rcpp) <br />
   sourceCpp("max_min.cpp")