QueryDrugPairMat <- function(data, y) {
	# start the clock
	ptm <- proc.time()
    drug.names <- rownames(data)
    drug.pair.mat <- matrix(0, nrow = nrow(data), ncol = nrow(data))
    # find the drugs without any targets
    drug.no.targets <- which(rowSums(as.matrix(tmp$x)) == 0)
    drug.idx <- which(rowSums(as.matrix(tmp$x)) != 0)
    drug.pair.mat[drug.no.targets, ] <- NA
    drug.pair.mat[, drug.no.targets] <- NA
    for (i in drug.idx) {
      cat("i is ", i, "\n")
      for (j in drug.idx) {
        drug.pair.mat[i, j] <- QueryDrugPair(data, y, c(i, j))$sens.pair.pred
      }
    
    
    }
    rownames(drug.pair.mat) <- drug.names
    colnames(drug.pair.mat) <- drug.names
    write.table(drug.pair.mat, file = "drug.pair.mat.csv", sep = ",")
	print(proc.time()-ptm)
}

QueryDrugPairMat1 <- function(data, y) {
	# start the clock
	ptm <- proc.time()
    drug.names <- rownames(data)
    row.num <- nrow(data)
    drug.pair.mat <- matrix(0, nrow = row.num, ncol = row.num)
    for (i in seq_len(row.num)) {
  	  drug.pair.mat[c(i:row.num), i] <- sapply(c(i:row.num), function(x) QueryDrugPair(data, y, c(i, x)))
    }
    rownames(drug.pair.mat) <- drug.names
    colnames(drug.pair.mat) <- drug.names
    write.table(drug.pair.mat, file = "drug.pair.mat.csv", sep = ",")
	print(proc.time()-ptm)
}
