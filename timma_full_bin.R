#' Predicting drug sensitivity with binary drug-target interaction data
#' 
#' A function to predict the drug sensitivity with binary drug-target interaction data using the 
#' one.sided TIMMA model
#' 
#' @param data the drug-target interaction data. See \code{\link{timma}}.
#' @param sens a drug sensitivity vector.
#' @param loo a logical value indicating whether to use the leave-one-out cross-validation in the model
#' selection process. By default, loo = TRUE. 
#' @return A list containing the following components:
#' \item{efficacy.mat}{the predicted efficacy matrix}
#' \item{error}{the prediction errors}
#' \item{prediction}{predicted drug sensitivity}
#' The difference between \code{\link{timmaModel}} and \code{\link{timmaBinary}} is \code{\link{timmaModel}} 
#' returns the predicted efficacy matrix of all possible target combinations while \code{\link{timmaBinary}}
#' not.
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @references Tang J, Karhinen L, Xu T, Szwajda A, Yadav B, Wennerberg K, Aittokallio T. 
#' Target inhibition networks: predicting selective combinations of druggable targets to block cancer 
#' survival pathways. PLOS Computational Biology 2013; 9: e1003226.
#' @examples
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' results <- RunTimmaFullBin(tyner_interaction_binary[, 1:6], tyner_sensitivity[, 1])

RunTimmaFullBin <- function(data, sens, loo = TRUE, averaging = "one.sided") {
    # parameter 1: data, drug with selected target profile    
    # parameter 2: sens, the actual efficacy for the drugs    
    # parameter 3: loo, flag for applying Leave-one-out or not
	  # paramtert 4: averaging, the averaging method
	
	# if the averaging method parameter is correct
	if (!(averaging %in% c("one.sided", "two.sided"))) {
		stop("Unkown averaging method!")
	}
    
	data <- as.matrix(data)
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
    
    # predicted error
    err.pred <- rep(NA, drug.num)
    # predicted efficacy
    pred <- rep(NA, drug.num)
    if (loo == FALSE) {
        pred <- M[drug.idx]
        err.pred <- abs(pred - sens)
        
    } else {
        for (i in 1:drug.num) {
            # remove drug i
            
            # get the dim info
            dim.IMd <- dim(IM.d)
            dim.IMd[3] <- dim.IMd[3] - 1
            
            IM.d.loo <- array(IM.d[, , -i], dim = dim.IMd)
            IM.subset.loo <- array(IM.subset[, , -i], dim = dim.IMd)
            IM.superset.loo <- array(IM.superset[, , -i], dim = dim.IMd)
            sens.loo <- sens[-i]
            
            
            M.d.loo <- sumcpp(IM.d.loo, dim.IMd[1], dim.IMd[2], dim.IMd[3])
            M.loo <- M.d.loo
            
            maxval.loo <- maxcpp(IM.superset.loo, dim.IMd[1], dim.IMd[2], dim.IMd[3])
            superset.max.loo <- maxval.loo$max
            superset.max.idx.loo <- maxval.loo$max.idx
            minval.loo <- mincpp(IM.subset.loo, dim.IMd[1], dim.IMd[2], dim.IMd[3])
            subset.min.loo <- minval.loo$min
            subset.min.idx.loo <- minval.loo$min.idx
            
            cell.max.avg <- is.nan(M.d.loo) & is.finite(superset.max.loo)
            cell.max.avg <- which(cell.max.avg == TRUE)
            
            cell.min.avg <- is.nan(M.d.loo) & is.finite(subset.min.loo)
            cell.min.avg <- which(cell.min.avg == TRUE)
            
            # does the removed drug need max averaging
            max.eq.identical <- which(cell.max.avg == drug.idx[i])
            # does the removed drug need min averaging
            min.eq.identical <- which(cell.min.avg == drug.idx[i])
            
            if (length(max.eq.identical) != 0 && length(min.eq.identical) == 0) {
                # index for the cell
                cell.max.idx <- cell.max.avg[max.eq.identical]
                row <- ((cell.max.idx - 1)%%rows) + 1
                col <- floor((cell.max.idx - 1)/rows) + 1
                drug.sub.cell <- !is.infinite(IM.superset.loo[row, col, ])
                # the drug index which achieves max sensitivity
                index <- superset.max.idx.loo[cell.max.idx]
                # the correspongding gray code for the drug with max sensitivity
                graycode.idx <- which(IM.d.loo[, , index] >= 0, arr.ind = TRUE)
                # find the supersets of S(index,:) in S that has smaller sensitivity
                supersets.small <- IM.subset.loo[graycode.idx[1], graycode.idx[2], ] < superset.max.loo[cell.max.idx]
                # find common item with drug.sub.cell and supersets.small
                cell.common <- which(drug.sub.cell & supersets.small)
                # max averaging
                if (length(cell.common) != 0) {
                  k <- 1
                  for (j in cell.common) {
                    superset.max.loo[cell.max.idx] <- (superset.max.loo[cell.max.idx] * k + sens.loo[j])/(k + 1)
                    k <- k + 1
                  }
                }
				if (averaging == "one.sided") {
					pred[i] <- superset.max.loo[drug.idx[i]]
				} else {
					pred[i] <- (superset.max.loo[drug.idx[i]] + 1)/2
				}
                
                
                err.pred[i] <- abs(pred[i] - sens[i])
                
            } else if (length(max.eq.identical) == 0 && length(min.eq.identical) != 0) {
                cell.min.idx <- cell.min.avg[min.eq.identical]
                row <- ((cell.min.idx - 1)%%rows) + 1
                col <- floor((cell.min.idx - 1)/rows) + 1
                drug.sub.cell <- !is.infinite(IM.subset.loo[row, col, ])
                index <- subset.min.idx.loo[cell.min.idx]
                # the correspongding gray code for the drug with max sensitivity
                graycode.idx <- which(IM.d.loo[, , index] >= 0, arr.ind = TRUE)
                # find the subsets of S(index,:) in S that has higher sensitivity
                subsets.high <- IM.superset.loo[graycode.idx[1], graycode.idx[2], ] > subset.min.loo[cell.min.idx]
                # find common item with drug.sub.cell and subsets.high
                cell.common <- which(drug.sub.cell & subsets.high)
                if (length(cell.common) != 0) {
                  # min averaging
                  k <- 1
                  for (j in cell.common) {
                    subset.min.loo[cell.min.idx] <- (subset.min.loo[cell.min.idx] * k + sens.loo[j])/(k + 1)
                    k <- k + 1
                  }
                }
				if (averaging == "one.sided") {
					pred[i] <- subset.min.loo[drug.idx[i]]
				} else {
					pred[i] <- (subset.min.loo[drug.idx[i]] + 0)/2
				}
                
                
                
                err.pred[i] <- abs(pred[i] - sens[i])
            } else if (length(max.eq.identical) != 0 && length(min.eq.identical) != 0) {
                # index for the cell
                cell.max.idx <- cell.max.avg[max.eq.identical]
                row <- ((cell.max.idx - 1)%%rows) + 1
                col <- floor((cell.max.idx - 1)/rows) + 1
                drug.sub.cell <- !is.infinite(IM.superset.loo[row, col, ])
                # the drug index which achieves max sensitivity
                index <- superset.max.idx.loo[cell.max.idx]
                # the correspongding gray code for the drug with max sensitivity
                graycode.idx <- which(IM.d.loo[, , index] >= 0, arr.ind = TRUE)
                # find the supersets of S(index,:) in S that has smaller sensitivity
                supersets.small <- IM.subset.loo[graycode.idx[1], graycode.idx[2], ] < superset.max.loo[cell.max.idx]
                # find common item with drug.sub.cell and supersets.small
                cell.common <- which(drug.sub.cell & supersets.small)
                if (length(cell.common) != 0) {
                  # max averaging
                  k <- 1
                  for (j in cell.common) {
                    superset.max.loo[cell.max.idx] <- (superset.max.loo[cell.max.idx] * k + sens.loo[j])/(k + 1)
                    k <- k + 1
                  }
                }
                cell.min.idx <- cell.min.avg[min.eq.identical]
                row <- ((cell.min.idx - 1)%%rows) + 1
                col <- floor((cell.min.idx - 1)/rows) + 1
                drug.sub.cell <- !is.infinite(IM.subset.loo[row, col, ])
                index <- subset.min.idx.loo[cell.min.idx]
                # the correspongding gray code for the drug with max sensitivity
                graycode.idx <- which(IM.d.loo[, , index] >= 0, arr.ind = TRUE)
                # find the subsets of S(index,:) in S that has higher sensitivity
                subsets.high <- IM.superset.loo[graycode.idx[1], graycode.idx[2], ] > subset.min.loo[cell.min.idx]
                # find common item with drug.sub.cell and subsets.high
                cell.common <- which(drug.sub.cell & subsets.high)
                if (length(cell.common) != 0) {
                  # min averaging
                  k <- 1
                  for (j in cell.common) {
                    subset.min.loo[cell.min.idx] <- (subset.min.loo[cell.min.idx] * k + sens.loo[j])/(k + 1)
                    
                    k <- k + 1
                  }
                }
                pred[i] <- (superset.max.loo[drug.idx[i]] + subset.min.loo[drug.idx[i]])/2
                err.pred[i] <- abs(pred[i] - sens[i])
                
            } else {
                # length(max.eq.identical)==0 && length(min.eq.identical)==0
                pred[i] <- M.loo[drug.idx[i]]
                err.pred[i] <- abs(pred[i] - sens[i])
            }
        }
        
    }
    return(list(efficacy.mat = M, error = err.pred, prediction = pred))
}
 
