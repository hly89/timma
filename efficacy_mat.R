#' Computing the predicted efficacy matrix 
#' 
#' A function to compute the predicted efficacy matrix. 
#' 
#' @param data the drug-target interaction data. See \code{\link{timma}}.
#' @param sens a drug sensitivity vector.
#' @param averaging a parameter to specify which one of the averaging algorithms will be applied in the model construction. 
#' By default, averaging = "one.sided", which is the original model construction algorithm. When averaging = "two.sided", 
#' a modified averaging algorithm will be used. These two variants only differ for the case where the minimization and 
#' maximization rules are not simultaneously satisfied. For example, for a queried target set if the supersets but not the 
#' subsets can be found in the training data, the one.sided algorithm will take the prediction from the averages on the 
#' supersets sensitivities using the minimization rule. The two.sided algorithm, however, will lower the predicted 
#' sensitivity by averaging it with 0, which is the theoretical lower boundary of the sensitivities that could be 
#' obtained in the subsets.
#' @return The predicted efficacy matrix without headers are returned and the predicted efficacy matrix with headers are saved as "efficacy.mat.csv".
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @references Tang J, Karhinen L, Xu T, Szwajda A, Yadav B, Wennerberg K, Aittokallio T. 
#' Target inhibition networks: predicting selective combinations of druggable targets to block cancer 
#' survival pathways. PLOS Computational Biology 2013; 9: e1003226.
#' @examples
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' results <- GetEfficacyMat(tyner_interaction_binary[, 1:6], tyner_sensitivity[, 1])
GetEfficacyMat <- function(data, sens, averaging = "one.sided") {
	# if the averaging method parameter is correct
	if (!(averaging %in% c("one.sided", "two.sided"))) {
		stop("Unkown averaging method!")
	}
    data <- as.matrix(data)
    # get target names
    target.names <- colnames(data)
    # get drug number
    drug.num <- nrow(data)
    # get target number
    target.num <- ncol(data)
    # get all possible gray code decimal
    graycode.dec <- GetGrayCodeDec(target.num)
    rows <- graycode.dec[[1]]
    cols <- graycode.dec[[2]]
    
    IM.d <- array(NA, dim = c(rows, cols, drug.num))
    IM.subset <- array(Inf, dim = c(rows, cols, drug.num))
    IM.superset <- array(-Inf, dim = c(rows, cols, drug.num))
    # index for the drug
    drug.idx <- rep(0, drug.num)
    data <- matrix(data, nrow = drug.num, ncol = target.num)
    for (i in 1:drug.num) {
        # get the decimal
        dec <- strtoi(paste(data[i, ], collapse = ""), base = 2)
        drug.idx[i] <- which(graycode.dec[[3]] == dec)
        mat.temp <- IM.d[, , i]
        mat.temp[drug.idx[i]] <- 1 * sens[i]
        IM.d[, , i] <- mat.temp
        
        # get the binary set: superset and subset
        super.sub.sets <- GetBinSet(data[i, ])
        # ismember function R version: match
        subset.idx <- match(super.sub.sets$subset, graycode.dec[[3]])
        subset.mat <- IM.subset[, , i]
        subset.mat[subset.idx] <- 1 * sens[i]
        IM.subset[, , i] <- subset.mat
        
        superset.idx <- match(super.sub.sets$superset, graycode.dec[[3]])
        superset.mat <- IM.superset[, , i]
        superset.mat[superset.idx] <- 1 * sens[i]
        IM.superset[, , i] <- superset.mat
    }
    M.d <- sumcpp(IM.d, rows, cols, drug.num)
    maxval <- maxcpp(IM.superset, rows, cols, drug.num)
    minval <- mincpp(IM.subset, rows, cols, drug.num)
    subset.min <- minval$min
    subset.min.idx <- minval$min.idx
    superset.max <- maxval$max
    superset.max.idx <- maxval$max.idx
    # find cell which needs maximization averaging
    cell.max.avg <- is.nan(M.d) & is.finite(superset.max)
    cell.max.avg <- which(cell.max.avg == TRUE)
    if (length(cell.max.avg) != 0) {
        for (i in cell.max.avg) {
            row <- ((i - 1)%%rows) + 1
            col <- floor((i - 1)/rows) + 1
            # the drug sets that are the subset of the cell
            drug.sub.cell <- !is.infinite(IM.superset[row, col, ])
            # the drug index which achieves max sensitivity
            index <- superset.max.idx[i]
            # the correspongding gray code for the drug with max sensitivity
            graycode.idx <- which(IM.d[, , index] >= 0, arr.ind = TRUE)
            # find the supersets of S(index,:) in S that has smaller sensitivity
            supersets.small <- IM.subset[graycode.idx[1], graycode.idx[2], ] < superset.max[i]
            # find common item with drug.sub.cell and supersets.small
            cell.common <- which(drug.sub.cell & supersets.small)
            if (length(cell.common) != 0) {
                # max averaging
                k <- 1
                for (j in cell.common) {
                  superset.max[i] <- (superset.max[i] * k + sens[j])/(k + 1)
                  k <- k + 1
                }
            }
        }
    }
    
    cell.min.avg <- is.nan(M.d) & is.finite(subset.min)
    cell.min.avg <- which(cell.min.avg == TRUE)
    if (length(cell.min.avg) != 0) {
        for (i in cell.min.avg) {
            row <- ((i - 1)%%rows) + 1
            col <- floor((i - 1)/rows) + 1
            # the drug sets that are the superset of the cell
            drug.sub.cell <- !is.infinite(IM.subset[row, col, ])
            # the drug index which achieves min sensitivity
            index <- subset.min.idx[i]
            # the correspongding gray code for the drug with max sensitivity
            graycode.idx <- which(IM.d[, , index] >= 0, arr.ind = TRUE)
            # find the subsets of S(index,:) in S that has higher sensitivity
            subsets.high <- IM.superset[graycode.idx[1], graycode.idx[2], ] > subset.min[i]
            # find common item with drug.sub.cell and subsets.high
            if (length(subsets.high) == 0) {
                cell.common <- vector("numeric")
            } else {
                cell.common <- which(drug.sub.cell & subsets.high)
            }
            if (length(cell.common) != 0) {
                # min averaging
                k <- 1
                for (j in cell.common) {
                  subset.min[i] <- (subset.min[i] * k + sens[j])/(k + 1)
                  k <- k + 1
                }
            }
        }
    }
    
    M <- M.d
	if (averaging == "one.sided") {
	    M[cell.max.avg] <- superset.max[cell.max.avg]
	    M[cell.min.avg] <- subset.min[cell.min.avg]
	} else {
	    M[cell.max.avg] <- (superset.max[cell.max.avg] + 1)/2
	    M[cell.min.avg] <- (subset.min[cell.min.avg] + 0)/2
	}
    
    
    # cels that not only have lower boundery and also have upper boundary
    max.min.avg.idx <- intersect(cell.max.avg, cell.min.avg)
    M[max.min.avg.idx] <- (superset.max[max.min.avg.idx] + subset.min[max.min.avg.idx])/2
	  # output the efficacy to csv file
    graycode <- GetGrayCodeMat(target.num)
    graycode.names <- GetGrayCodeNames(target.num, target.names, graycode$graycode.row, graycode$graycode.col)
    nr <- graycode.names$graycode.names.row
    nc <- t(graycode.names$graycode.names.col)
    mat.row <- nrow(nr) + nrow(nc)
    mat.col <- ncol(nr) + ncol(nc)
    efficacy.mat <- array("", dim = c(mat.row, mat.col))
    efficacy.mat[(nrow(nc) + 1):mat.row, 1:ncol(nr)] <- nr
    efficacy.mat[1:nrow(nc), (ncol(nr) + 1):mat.col] <- nc
    efficacy.mat[(nrow(nc) + 1):mat.row, (ncol(nr) + 1):mat.col] <- M
    write.table(efficacy.mat, file = "efficacy.mat.csv", sep = ",", col.names = FALSE, row.names = FALSE)
    dir <- getwd()
    cat("The efficacy matrix is saved as efficacy.mat.csv in ", dir, "\n")
    return(M)
}