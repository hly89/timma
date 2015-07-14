QueryDrugPair <- function(data, sens, drug.pair = c()) {
	if (length(drug.pair) != 2 ) {
		stop("Error: not a pair of drugs.")
	}
	
	if (class(drug.pair) == "character") {
		drug.pair.idx <- match(drug.pair, rownames(data))
		if (length(drug.pair.idx) != 2) {
			stop("The queried drug pair is not found.")
		}
	} else if (class(drug.pair) == "numeric") {
		drug.pair.idx <- drug.pair
	} else {
		stop("Error: unknown data type for drug pair")
	}
	# the targets names
	targets.names <- dimnames(data)[[2]]
	# get the targets for the first drug
	drug1.targets.idx <- which(data[drug.pair.idx[1],] == 1)
	if (length(drug1.targets.idx) == 0) {
		stop("No targets for the first drug.")
	}
	# get the targets for the second drug
	drug2.targets.idx <- which(data[drug.pair.idx[2],] == 1)
	if (length(drug2.targets.idx) == 0) {
		stop("No targets for the second drug.")
	}
	
	sens.pred <- c()
	for (i in drug1.targets.idx) {
		for (j in drug2.targets.idx){
			if (i == j) {
				timma.efficacy.mat <- GetEfficacyMat(data[, i], sens)
				sens.pred <- c(sens.pred, timma.efficacy.mat[1,2])
			} else {
				timma.efficacy.mat <- GetEfficacyMat(data[, c(i,j)], sens)
				sens.pred <- c(sens.pred, timma.efficacy.mat[2,2])
			}
		
			
		}
	}

	#score.test <- c()
	# get sens for single drugs
	sens.drug1 <- sens[drug.pair.idx[1]]
	sens.drug2 <- sens[drug.pair.idx[2]]
	targets.list <- unique(c(drug1.targets.idx, drug2.targets.idx))
	graycode.mat <- GetGrayCodeMat(length(targets.list))
	row.idx <- which.max(graycode.mat$dec.row)
	col.idx <- which.min(graycode.mat$dec.col)
	score.test <- GetEfficacyMat(data[, targets.list], sens)[row.idx,col.idx]
	return(list(sens.pair.pred = mean(sens.pred), sens.drug1 = sens.drug1, sens.drug2 = sens.drug2, score.test = score.test))
	
	
		
	
}