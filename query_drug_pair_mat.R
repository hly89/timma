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
	print(proc.time() - ptm)
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
	print(proc.time() - ptm)
}

QueryDrugPairList <- function (data, y) {
	# start the clock
	ptm <- proc.time()
	drug.names <- rownames(data)
	row.num <- nrow(data)
	list.options <- combn(c(1:row.num), 2)
	drug.pair.list <- matrix(0, nrow = ncol(list.options), ncol = 6)
	colnames(drug.pair.list) <- c("Drug1", "Drug2", "Prediction Default", "Prediction HSA", "Drug1 Sensitivity", "Drug2 Sensitivity")
	for (i in seq_len(ncol(list.options))) {
		drug.pair.list[i, 1:2] <- drug.names[list.options[, i]]
		drug.pair.list[i, 5:6] <- y[list.options[, i]]
		sens.pred <- QueryDrugPair(data, y, drug.pair = list.options[, i])
		drug.pair.list[i, 3:4] <- c(sens.pred$sens.pair.pred, sens.pred$sens.pred.hsa)
	}
	
	# same drug pair like c(1,1)
	drug.pair.same <- matrix(0, nrow = row.num, ncol = 6)
	colnames(drug.pair.same) <- colnames(drug.pair.list)
	for (i in seq_len(row.num)) {
		drug.pair.same[i, 1:2] <- drug.names[i]
		drug.pair.same[i, 5:6] <- y[i]
		sens.pred <- QueryDrugPair(data, y, drug.pair = c(i, i))
		drug.pair.same[i, 3:4] <- c(sens.pred$sens.pair.pred, sens.pred$sens.pred.hsa)
	}
	drug.pair.list <- rbind(drug.pair.list, drug.pair.same)
	write.table(drug.pair.list, file = "drug.pair.list.csv", sep = ",", row.names = FALSE)
	write.table(drug.pair.list, file = "drug.pair.list.txt", sep = "\t", row.names = FALSE)
	print(proc.time() - ptm)
}