#' Query the predicted sensitivity for a singe drug pair
#'
#' A function to query the predicted sensitity for a single drug pair
#' @param data a drug-target interaction matrix. Row names are drug names and column names are target names.
#' @param sens a normalized drug sensitivity vector.
#' @param drug.pair a vector with either drug names or drug index in the drug-target interaction matrix.
#' 
#' @return A list containing the following components:
#' \item {sens.pair.pred} defualt predicted sensitivity for single drug pair
#' \item {sens.drug1} the sensitivity for the first drug
#' \item {sens.drug2} the sensitivty for the second drug
#' \item {sens.pred.hsa} the predicted sensitivity with hsa model for single drug pair
#'
#' @author Liye He \email{liye.he@@helsinki.fi}
#' @examples
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' results <- QuerySingleDrugPair(tyner_interaction_binary, tyner_sensitivity[,1], c(1,2))
QuerySingleDrugPair <- function(data, sens, drug.pair = c()) {
	if (length(drug.pair) != 2 ) {
		stop("Error: not a pair of drugs.")
	}
	
	if (class(drug.pair) == "character") {
		drug.pair.idx <- match(drug.pair, rownames(data))
		if (length(drug.pair.idx) != 2) {
			stop("The queried drug pair is not found.")
		}
	} else if (class(drug.pair) == "numeric" || class(drug.pair) == "integer") {
		drug.pair.idx <- drug.pair
	} else {
		stop("Error: unknown data type for drug pair")
	}
	# the targets names
	targets.names <- dimnames(data)[[2]]
	# get the targets for the first drug
	drug1.targets.idx <- which(data[drug.pair.idx[1],] == 1)
	if (length(drug1.targets.idx) == 0) {
		cat("No targets for the first drug. \n")
	  sens.drug1 <- sens[drug.pair.idx[1]]
	  sens.drug2 <- sens[drug.pair.idx[2]]
	  return(list(sens.pair.pred = NA, sens.drug1 = sens.drug1, sens.drug2 = sens.drug2))
	}
	# get the targets for the second drug
	drug2.targets.idx <- which(data[drug.pair.idx[2],] == 1)
	if (length(drug2.targets.idx) == 0) {
		cat("No targets for the second drug. \n")
	  sens.drug1 <- sens[drug.pair.idx[1]]
	  sens.drug2 <- sens[drug.pair.idx[2]]
	  return(list(sens.pair.pred = NA, sens.drug1 = sens.drug1, sens.drug2 = sens.drug2))
	}
	# get sens for single drugs
	sens.drug1 <- sens[drug.pair.idx[1]]
	sens.drug2 <- sens[drug.pair.idx[2]]
	
	sens.pred <- c()
	sens.pred.hsa <- c()
	for (i in drug1.targets.idx) {
		for (j in drug2.targets.idx){
			if (i == j) {
				timma.efficacy.mat <- GetEfficacyMat(data[, i], sens)
				sens.pred <- c(sens.pred, timma.efficacy.mat[1,2])
				sens.pred.hsa <- c(sens.pred.hsa, max(sens.drug1, sens.drug2))
			} else {
			  
				timma.efficacy.mat <- GetEfficacyMat(data[, c(i,j)], sens)
				sens.pred <- c(sens.pred, timma.efficacy.mat[2,2])
				# using the HSA model to score
				hsa.score <- timma.efficacy.mat[2,2] - max(timma.efficacy.mat[2,1], timma.efficacy.mat[1,2])                                                                      
				sens.pred.hsa <- c(sens.pred.hsa, hsa.score + max(sens.drug1, sens.drug2))
			}
		
			
		}
	}

	#score.test <- c()
	
	#targets.list <- unique(c(drug1.targets.idx, drug2.targets.idx))
	#graycode.mat <- GetGrayCodeMat(length(targets.list))
	#row.idx <- which.max(graycode.mat$dec.row)
	#col.idx <- which.min(graycode.mat$dec.col)
	#score.test <- GetEfficacyMat(data[, targets.list], sens)[row.idx,col.idx]
	return(list(sens.pair.pred = mean(sens.pred), sens.drug1 = sens.drug1, sens.drug2 = sens.drug2, sens.pred.hsa = mean(sens.pred.hsa)))
	
	
		
	
}