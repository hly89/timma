#'Prediction in the search space with one.sided TIMMA model 
#'
#'A function to return the prediction error in the search space for sffs
#'
#'@param profile_k current selected drug-target interaction data
#'@param space the search space returned by \code{\link{searchSpace}} function
#'@param sens drug sensitivity data
#'@param loo a logical value indicating whether to use the leave-one-out cross-validation in the model
#' selection process. By default, loo = TRUE. 
#'@return the prediction error
#'@author Liye He \email{liye.he@@helsinki.fi}
#' @examples 
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' num<-length(tyner_sensitivity[,1])
#' k_set<-rep(0, dim(tyner_interaction_binary)[2])
#' k_set[c(1,2,3)]<-1
#' space<-searchSpace(num, k_set, tyner_interaction_binary, tyner_sensitivity[,1])
#' profile_k<-tyner_interaction_binary[, which(k_set==1)]
#' error<-GetPredErr(profile_k, space, tyner_sensitivity[,1])

GetPredErr <- function(data, space, sens, loo = TRUE, averaging = "one.sided") {
	# if the averaging method parameter is correct
	if (!(averaging %in% c("one.sided", "two.sided"))) {
		stop("Unkown averaging method!")
	}
	
    # space is 3d array
    dim.info <- dim(space$space.identical)
    rows <- dim.info[1]
    cols <- dim.info[2]
    
    IM.d <- array(NA, dim = dim.info[1:2])
    IM.superset <- array(-Inf, dim = dim.info[1:2])
    IM.subset <- array(Inf, dim = dim.info[1:2])
    
    idx.identical <- rep(0, rows)
    for (i in 1:rows) {
        index <- data[i] + 1
        IM.d[i, ] <- space$space.identical[i, , index]
        IM.superset[i, ] <- space$space.superset[i, , index]
        IM.subset[i, ] <- space$space.subset[i, , index]
        idx.identical[i] <- which((!is.na(IM.d[i, ])) == TRUE)
    }
    
    M.d <- sumcpp1(IM.d, rows, cols)
    maxval <- maxcpp1(IM.superset, rows, cols)
    minval <- mincpp1(IM.subset, rows, cols)
    subset.min <- minval$min
    min.index <- minval$min.idx
    superset.max <- maxval$max
    max.index <- maxval$max.idx
    
    
    
    # find cell which needs maximization averaging
    cell.max.avg <- is.nan(M.d) & is.finite(superset.max)  # is.nan or is.na????????
    cell.max.avg <- which(cell.max.avg == TRUE)
    if (length(cell.max.avg) != 0) {
        for (i in cell.max.avg) {
            # row<-((i-1) %% rows) + 1 col<-floor((i-1) / rows)+1 print(ii) the drug sets that are the subset of the
            # cell
            drug.sub.cell <- !is.infinite(IM.superset[, i])
            # the drug index which achieves max sensitivity
            index <- max.index[i]
            # the dec of the drug with max sensitivity
            maxsens.dec.idx <- idx.identical[index]
            
            # find the supersets of S(index,:) in S that has smaller sensitivity
            supersets.small <- IM.subset[, maxsens.dec.idx] < superset.max[i]
            
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
    
	# find cell which needs minimization averaging
    cell.min.avg <- is.nan(M.d) & is.finite(subset.min)
    cell.min.avg <- which(cell.min.avg == TRUE)
    if (length(cell.min.avg) != 0) {
        for (i in cell.min.avg) {
            # row<-((i-1) %% rows) + 1 col<-floor((i-1) / rows)+1 the drug sets that are the superset of the cell
            drug.sub.cell <- !is.infinite(IM.subset[, i])
            # the drug index which achieves min sensitivity
            index <- min.index[i]
            
            # the dec of the drug with min sensitivity
            minsens.dec.idx <- idx.identical[index]
            
            # find the subsets of S(index,:) in S that has higher sensitivity
            subsets.high <- IM.superset[, minsens.dec.idx] > subset.min[i]
            # find common item with drug.sub.cell and supersets.small
            if (length(subsets.high) == 0) {
                cell.common <- vector("numeric") # why?
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
    
    # cells that not only have lower boundery and also have upper boundary
    max.min.avg.idx <- intersect(cell.max.avg, cell.min.avg)
    M[max.min.avg.idx] <- (superset.max[max.min.avg.idx] + subset.min[max.min.avg.idx])/2
    
    # predicted error
    err.pred <- rep(NA, rows)
    # predicted efficacy
    pred <- rep(NA, rows)
    if (loo == FALSE) {
        pred <- M[idx.identical]
        err.pred <- abs(pred - sens)
        
    } else {
        for (i in 1:rows) {
            # remove drug i, namely remove the i-th row
            
            # get the dim info dim.IMd<-dim(IM_d) dim.IMd[3]<-dim.IMd[3]-1
            dim.IMd <- c(rows - 1, cols)
            
            IM.d.loo <- array(IM.d[-i, ], dim = dim.IMd)
            IM.subset.loo <- array(IM.subset[-i, ], dim = dim.IMd)
            IM.superset.loo <- array(IM.superset[-i, ], dim = dim.IMd)
            
            sens.loo <- sens[-i]
            drug.idx.loo <- idx.identical[-i]
            
            
            M.d.loo <- sumcpp1(IM.d.loo, rows - 1, cols)
            M.loo <- M.d.loo

            
            maxval <- maxcpp1(IM.superset.loo, rows - 1, cols)
            minval <- mincpp1(IM.subset.loo, rows - 1, cols)
            subset.min.loo <- minval$min
            min.idx.loo <- minval$min.idx
            superset.max.loo <- maxval$max
            max.idx.loo <- maxval$max.idx
            
            
            cell.max.avg <- is.nan(M.d.loo) & is.finite(superset.max.loo)
            cell.max.avg <- which(cell.max.avg == TRUE)
            
            cell.min.avg <- is.nan(M.d.loo) & is.finite(subset.min.loo)
            cell.min.avg <- which(cell.min.avg == TRUE)
            
            # does the removed drug need max averaging max.eq.identical<-which(cell.max.avg==drug_index[i])
            max.eq.identical <- which(cell.max.avg == idx.identical[i])
            # does the removed drug need min averaging min.eq.identical<-which(cell.min.avg==drug_index[i])
            min.eq.identical <- which(cell.min.avg == idx.identical[i])
            
            if (length(max.eq.identical) != 0 && length(min.eq.identical) == 0) {
                # index for the cell
                cell.max.idx <- cell.max.avg[max.eq.identical]
                
                drug.sub.cell <- !is.infinite(IM.superset.loo[, cell.max.idx])
                # the drug index which achieves max sensitivity
                drug.maxsens.idx <- max.idx.loo[cell.max.idx]
                
                # the index of the dec of the drug with max sensitivity
                maxsens.dec.idx <- drug.idx.loo[drug.maxsens.idx]
                
                # find the supersets of S(index,:) in S that has smaller sensitivity
                supersets.small <- IM.subset.loo[, maxsens.dec.idx] < superset.max.loo[cell.max.idx]
                
                # find common item with drug.sub.cell and supersets.small
                cell.common <- which(drug.sub.cell & supersets.small)
                
                if (length(cell.common) != 0) {
                  k <- 1
                  for (j in cell.common) {
                    superset.max.loo[cell.max.idx] <- (superset.max.loo[cell.max.idx] * k + sens.loo[j])/(k + 
                      1)
                    k <- k + 1
                  }
                }
                if (averaging == "one.sided") {
                	pred[i] <- superset.max.loo[idx.identical[i]]
                } else {
                	pred[i] <- (superset.max.loo[idx.identical[i]] + 1)/2
                }
                
                err.pred[i] <- abs(pred[i] - sens[i])
                
            } else if (length(max.eq.identical) == 0 && length(min.eq.identical) != 0) {
                cell.min.idx <- cell.min.avg[min.eq.identical]
                
                drug.sub.cell <- !is.infinite(IM.subset.loo[, cell.min.idx])
                
                index <- min.idx.loo[cell.min.idx]
                minsens.dec.idx <- drug.idx.loo[index]
                
                # find the supersets of S(index,:) in S that has smaller sensitivity
                supersets.small <- IM.superset.loo[, minsens.dec.idx] > subset.min.loo[cell.min.idx]
                # find common item with drug.sub.cell and supersets.small
                cell.common <- which(drug.sub.cell & supersets.small)
                if (length(cell.common) != 0) {
                  # max averaging
                  k <- 1
                  for (j in cell.common) {
                    subset.min.loo[cell.min.idx] <- (subset.min.loo[cell.min.idx] * k + sens.loo[j])/(k + 1)
                    k <- k + 1
                  }
                }
                if (averaging == "one.sided") {
                	pred[i] <- subset.min.loo[idx.identical[i]]
                } else {
                	pred[i] <- (subset.min.loo[idx.identical[i]] + 0)/2
                }
                
                err.pred[i] <- abs(pred[i] - sens[i])
            } else if (length(max.eq.identical) != 0 && length(min.eq.identical) != 0) {
                
                
                cell.max.idx <- cell.max.avg[max.eq.identical]
                
                drug.sub.cell <- !is.infinite(IM.superset.loo[, cell.max.idx])
                # the drug index which achieves max sensitivity
                drug.maxsens.idx <- max.idx.loo[cell.max.idx]
                
                # the dec of the drug with max sensitivity
                maxsens.dec.idx <- drug.idx.loo[drug.maxsens.idx]
                
                # find the supersets of S(index,:) in S that has smaller sensitivity
                supersets.small <- IM.subset.loo[, maxsens.dec.idx] < superset.max.loo[cell.max.idx]
                
                # find common item with drug.sub.cell and supersets.small
                cell.common <- which(drug.sub.cell & supersets.small)
                
                if (length(cell.common) != 0) {
                  k <- 1
                  for (j in cell.common) {
                    superset.max.loo[cell.max.idx] <- (superset.max.loo[cell.max.idx] * k + sens.loo[j])/(k + 
                      1)
                    k <- k + 1
                  }
                }
                
                cell.min.idx <- cell.min.avg[min.eq.identical]
                
                drug.sub.cell <- !is.infinite(IM.subset.loo[, cell.min.idx])
                
                drug.minsens.idx <- min.idx.loo[cell.min.idx]
                minsens.dec.idx <- drug.idx.loo[drug.minsens.idx]
                
                # find the supersets of S(index,:) in S that has smaller sensitivity
                supersets.small <- IM.superset.loo[, minsens.dec.idx] > subset.min.loo[cell.min.idx]
                # find common item with drug.sub.cell and supersets.small
                cell.common <- which(drug.sub.cell & supersets.small)
                if (length(cell.common) != 0) {
                  # max averaging
                  k <- 1
                  for (j in cell.common) {
                    subset.min.loo[cell.min.idx] <- (subset.min.loo[cell.min.idx] * k + sens.loo[j])/(k + 1)
                    k <- k + 1
                  }
                }
                
                
                pred[i] <- (superset.max.loo[idx.identical[i]] + subset.min.loo[idx.identical[i]])/2
                err.pred[i] <- abs(pred[i] - sens[i])
                
            } else {
                # length(max.eq.identical)==0 && length(min.eq.identical)==0
                pred[i] <- M.loo[idx.identical[i]]
                err.pred[i] <- abs(pred[i] - sens[i])
            }
        }
        
    }
    return(err.pred)
} 
