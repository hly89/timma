#' Query the predicted sensitivities for all possible drug pairs in the data
#'
#' A function to query the predicted sensitivities for all possible drug pairs in the data
#' @param data a drug-target interaction matrix. Row names are drug names and column names are target names.
#' @param y a normalized drug sensitivity vector.
#'
#' @return Two files are generated:
#' \item drug.pair.list.csv a csv file contains the prediction for all possible drug pairs 
#' \item drug.pair.list.txt a txt file contains the prediction for all possible drug pairs
#'
#' @author Liye He \email{liye.he@@helsinki.fi}
#' @examples
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' results <- QueryAllDrugPairs(tyner_interaction_binary, tyner_sensitivity[,1])
QueryAllDrugPairs <- function (data, y) {
	# start the clock
	ptm <- proc.time()
	drug.names <- rownames(data)
	row.num <- nrow(data)
	list.options <- combn(c(1:row.num), 2)
	drug.pair.list <- matrix(0, nrow = ncol(list.options), ncol = 7)
	colnames(drug.pair.list) <- c("Group", "Drug1", "Drug2", "Prediction Default", "Prediction HSA", "Drug1 Sensitivity", "Drug2 Sensitivity")
	for (i in seq_len(ncol(list.options))) {
		drug1.name <- drug.names[list.options[1, i]]
		drug2.name <- drug.names[list.options[2, i]]
		if (substr(drug1.name, 1, 1) < substr(drug2.name, 1, 1)) {
			name.group <- paste(drug1.name, drug2.name)
		} else {
			name.group <- paste(drug2.name, drug1.name)
		}
		drug.pair.list[i, 1] <- name.group
		drug.pair.list[i, 2:3] <- drug.names[list.options[, i]]
		drug.pair.list[i, 6:7] <- y[list.options[, i]]
		sens.pred <- QuerySingleDrugPair(data, y, drug.pair = list.options[, i])
		drug.pair.list[i, 4:5] <- c(sens.pred$sens.pair.pred, sens.pred$sens.pred.hsa)
	}
	
	# same drug pair like c(1,1)
	drug.pair.same <- matrix(0, nrow = row.num, ncol = 7)
	colnames(drug.pair.same) <- colnames(drug.pair.list)
	for (i in seq_len(row.num)) {
		drug.pair.same[i, 1] <- paste(drug.names[i], drug.names[i])
		drug.pair.same[i, 2:3] <- drug.names[i]
		drug.pair.same[i, 6:7] <- y[i]
		sens.pred <- QuerySingleDrugPair(data, y, drug.pair = c(i, i))
		drug.pair.same[i, 4:5] <- c(sens.pred$sens.pair.pred, sens.pred$sens.pred.hsa)
	}
	drug.pair.list <- rbind(drug.pair.list, drug.pair.same)
	# sort the list by group
	order.idx <- order(drug.pair.list[, 1])
	drug.pair.list <- drug.pair.list[order.idx, ]
	write.table(drug.pair.list, file = "drug.pair.list.csv", sep = ",", row.names = FALSE)
	write.table(drug.pair.list, file = "drug.pair.list.txt", sep = "\t", row.names = FALSE)
	print(proc.time() - ptm)
}

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
        drug.pair.mat[i, j] <- QuerySingleDrugPair(data, y, c(i, j))$sens.pair.pred
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
  	  drug.pair.mat[c(i:row.num), i] <- sapply(c(i:row.num), function(x) QuerySingleDrugPair(data, y, c(i, x)))
    }
    rownames(drug.pair.mat) <- drug.names
    colnames(drug.pair.mat) <- drug.names
    write.table(drug.pair.mat, file = "drug.pair.mat.csv", sep = ",")
	print(proc.time() - ptm)
}

