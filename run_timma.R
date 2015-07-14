 #' Predicting drug sensitivity with binary drug-target interaction data
#' 
#' A function to predict the drug sensitivity with binary drug-target interaction data using the 
#' original maximization and minimization rules
#' 
#' @param data the drug-target interaction data. See \code{\link{timma}}.
#' @param sens a drug sensitivity vector.
#' @param loo a logical value indicating whether to use the leave-one-out cross-validation in the model
#' selection process. By default, loo = TRUE. 
#' @return A list containing the following components:
#' \item{efficacy.mat}{the predicted efficacy for target combinations that can be found from the training data}
#' \item{error}{the prediction errors}
#' \item{prediction}{predicted drug sensitivity}
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @references Tang J, Karhinen L, Xu T, Szwajda A, Yadav B, Wennerberg K, Aittokallio T. 
#' Target inhibition networks: predicting selective combinations of druggable targets to block cancer 
#' survival pathways. PLOS Computational Biology 2013; 9: e1003226.
#' @examples 
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' results<-RunTimma(tyner_interaction_binary[, 1:6], tyner_sensitivity[,1])
RunTimma <- function(data, sens, loo = TRUE, averaging = "one.sided", class = 2, weighted = FALSE) {
    # parameter 1: data, drug with selected target profile    
    # parameter 2: sens, the actual efficacy for the drugs    
    # parameter 3: loo, flag for applying Leave-one-out or not
  
  
	# if the averaging method parameter is correct
	if (!(averaging %in% c("one.sided", "two.sided"))) {
		stop("Unkown averaging method!")
	}  
    # get target numbers
    target.num <- ncol(as.matrix(data))
    # get drug numbers
    drug.num <- nrow(as.matrix(data))
    
    data <- matrix(data, nrow = drug.num, ncol = target.num)
    data.unique <- unique(data)
    dec.data.unique <- apply(data.unique, 1, function(x) strtoi(paste(x, collapse = ""), base = class))
    dec <- apply(data, 1, function(x) strtoi(paste(x, collapse = ""), base = class))
    # for identical
    col.num <- length(dec.data.unique)
    # index for the drug
    idx.identical <- sapply(dec, function(x) which(dec.data.unique == x))
    IM.d <- array(NA, dim = c(drug.num, col.num))
    IM.subset <- array(Inf, dim = c(drug.num, col.num))
    IM.superset <- array(-Inf, dim = c(drug.num, col.num))
	if (weighted) {
		# weight for the subsets
		IM.sub.weight <- array(0, dim = c(drug.num, col.num))
		# weight for the supersets
		IM.sup.weight <- array(0, dim = c(drug.num, col.num))
	}
    for (i in 1:drug.num) {
        # get the decimal
        IM.d[i, idx.identical[i]] <- 1 * sens[i]
        
        # get the binary set: superset and subset
		if (class == 2) {
			super.sub.sets <- GetBinSet(data[i, ])
		} else if (class >2 && weighted == FALSE) {
			super.sub.sets <- GetBinSetMulti(data[i, ], data)
		} else {
			super.sub.sets <- GetBinSetMulti(data[i, ], data, weighted = TRUE)
		}
        
        # ismember function R version: match subset.idx<-match(super.sub.sets@subset, dec_graycode[[3]])
		if (class == 2) {
	        subset.idx <- dec.data.unique %in% super.sub.sets$subset
	        IM.subset[i, subset.idx] <- sens[i]
		} else {
			if(length(super.sub.sets$subset)!=0){
				subset.idx <- dec.data.unique %in% dec[super.sub.sets$subset]
				IM.subset[i, subset.idx] <- sens[i]
				if (weighted) {
					for(each in which(subset.idx == TRUE)){
						IM.sub.weight[i, each] <- super.sub.sets$subw[which(dec[super.sub.sets$subset] == dec.data.unique[each])[1]]
					}
				}
			}
		}
        
        
		if (class == 2) {
	        superset.idx <- dec.data.unique %in% super.sub.sets$superset
	        IM.superset[i, superset.idx] <- sens[i]
		} else {
	        if(length(super.sub.sets$superset)!=0){
		        superset.idx <- dec.data.unique %in% dec[super.sub.sets$superset]
		        IM.superset[i, superset.idx] <- sens[i]
				if (weighted) {
					for(each in which(superset.idx == TRUE)){
						IM.sup.weight[i, each] <- super.sub.sets$supw[which(dec[super.sub.sets$superset] == dec.data.unique[each])[1]]
					}
				}
	        }
		}
        
        
    }
	if (class == 2) {
		M.d <- sumcpp1(IM.d, drug.num, col.num)
	} else {
		M.d <- colMeans(IM.d, na.rm=TRUE)
	}
    
    maxval <- maxcpp1(IM.superset, drug.num, col.num)
    minval <- mincpp1(IM.subset, drug.num, col.num)
    subset.min <- minval$min
    subset.min.idx <- minval$min.idx
    superset.max <- maxval$max
    superset.max.idx <- maxval$max.idx
    
    
    
    # find cell which needs maximization averaging
    cell.max.avg <- is.nan(M.d) & is.finite(superset.max)
    cell.max.avg <- which(cell.max.avg == TRUE)
    if (length(cell.max.avg) != 0) {
        for (i in cell.max.avg) {
            
            # the drug sets that are the subset of the cell
            drug.sub.cell <- !is.infinite(IM.superset[, i])
            # the drug index which achieves max sensitivity
            index <- superset.max.idx[i]
            # the dec of the drug with max sensitivity
            maxsens.dec.idx <- idx.identical[index]
            
            # find the supersets of S(index,:) in S that has smaller sensitivity
            supersets.small <- IM.subset[, maxsens.dec.idx] < superset.max[i]
            
            # find common item with drug.sub.cell and supersets_small
            cell.common <- which(drug.sub.cell & supersets.small)
            
            if (length(cell.common) != 0) {
                # max averaging
				if (weighted) {
					total <- sum(1/IM.sup.weight[cell.common, i]) + 1/IM.sup.weight[index, i]
					superset.max[i] <- 1/IM.sup.weight[index, i]/total*sens[index] + sum(1/IM.sup.weight[cell.common, i]/total*sens[cell.common])
				} else {
	                k <- 1
	                for (j in cell.common) {
	                  superset.max[i] <- (superset.max[i] * k + sens[j])/(k + 1)
	                  k <- k + 1
	                }
				}
                
            }
        }
    }
    
    cell.min.avg <- is.nan(M.d) & is.finite(subset.min)
    cell.min.avg <- which(cell.min.avg == TRUE)
    if (length(cell.min.avg) != 0) {
        for (i in cell.min.avg) {
            # the drug sets that are the superset of the cell
            drug.sub.cell <- !is.infinite(IM.subset[, i])
            # the drug index which achieves min sensitivity
            index <- subset.min.idx[i]
            
            # the dec of the drug with min sensitivity
            minsens.dec.idx <- idx.identical[index]
            
            # find the subsets of S(index,:) in S that has higher sensitivity
            subsets.high <- IM.superset[, minsens.dec.idx] > subset.min[i]
            # find common item with drug.sub.cell and supersets_small
            if (length(subsets.high) == 0) {
                cell.common <- vector("numeric")
            } else {
                cell.common <- which(drug.sub.cell & subsets.high)
            }
            if (length(cell.common) != 0) {
                # min averaging
				if (weighted) {
					total <- sum(1/IM.sub.weight[cell.common, i]) + 1/IM.sub.weight[index, i]
					subset.min[i] <- 1/IM.sub.weight[index, i]/total*sens[index] + sum(1/IM.sub.weight[cell.common, i]/total*sens[cell.common])
				} else {
	                k <- 1
	                for (j in cell.common) {
	                  subset.min[i] <- (subset.min[i] * k + sens[j])/(k + 1)
	                  k <- k + 1
	                }
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
        pred <- M[idx.identical]
        err.pred <- abs(pred - sens)
        
    } else {
        for (i in 1:drug.num) {
            # remove drug i, namely remove the i-th row
            
            # get the dim info
            dim.IMd <- c(drug.num - 1, col.num)
            
            IM.d.loo <- array(IM.d[-i, ], dim = dim.IMd)
            IM.subset.loo <- array(IM.subset[-i, ], dim = dim.IMd)
            IM.superset.loo <- array(IM.superset[-i, ], dim = dim.IMd)
			if (weighted) {
				IM.sub.weight.loo <- array(IM.sub.weight[-i, ], dim = dim.IMd)
				IM.sup.weight.loo <- array(IM.sup.weight[-i, ], dim = dim.IMd)
			}
            sens.loo <- sens[-i]
            drug.idx.loo <- idx.identical[-i]
			if (class == 2) {
				M.d.loo <- sumcpp1(IM.d.loo, drug.num - 1, col.num)
			} else {
	            M.d.loo <- M.d
				M.d.loo[idx.identical[i]] <- mean(IM.d.loo[, idx.identical[i]], na.rm=TRUE)
			}
            
            M.loo <- M.d.loo
            
            maxval <- maxcpp1(IM.superset.loo, drug.num - 1, col.num)
            minval <- mincpp1(IM.subset.loo, drug.num - 1, col.num)
            subset.min.loo <- minval$min
            min.idx.loo <- minval$min.idx
            superset.max.loo <- maxval$max
            max.idx.loo <- maxval$max.idx
            
            
            cell.max.avg <- is.nan(M.d.loo) & is.finite(superset.max.loo)
            cell.max.avg <- which(cell.max.avg == TRUE)
            
            cell.min.avg <- is.nan(M.d.loo) & is.finite(subset.min.loo)
            cell.min.avg <- which(cell.min.avg == TRUE)
            
            # does the removed drug need max averaging
            max.eq.identical <- which(cell.max.avg == idx.identical[i])
            # does the removed drug need min averaging
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
                
                # find common item with drug.sub.cell and supersets_small
                cell.common <- which(drug.sub.cell & supersets.small)
                if (length(cell.common) != 0) {
                  # max averaging
				  if (weighted) {
	                  total <- sum(1/IM.sup.weight.loo[cell.common, cell.max.idx]) + 1/IM.sup.weight.loo[drug.maxsens.idx, cell.max.idx]
	  				  superset.max.loo[cell.max.idx] <- 1/IM.sup.weight.loo[drug.maxsens.idx, cell.max.idx]/total*sens.loo[drug.maxsens.idx] + sum(1/IM.sup.weight.loo[cell.common, cell.max.idx]/total*sens.loo[cell.common])
				  } else {
	                  k <- 1
	                  for (j in cell.common) {
	                    superset.max.loo[cell.max.idx] <- (superset.max.loo[cell.max.idx] * k + sens.loo[j])/(k + 1)
	                    k <- k + 1
	                  }
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
                
                # find the subsets of S(index,:) in S that has higher sensitivity
                subsets.high <- IM.superset.loo[, minsens.dec.idx] > subset.min.loo[cell.min.idx]
                # find common item with drug.sub.cell and supersets_small
                cell.common <- which(drug.sub.cell & subsets.high)
                if (length(cell.common) != 0) {
                  # min averaging
				  if (weighted) {
	                  total <- sum(1/IM.sub.weight.loo[cell.common, cell.min.idx]) + 1/IM.sub.weight.loo[index, cell.min.idx]
					  subset.min.loo[cell.min.idx] <- 1/IM.sub.weight.loo[index, cell.min.idx]/total * sens.loo[index] + sum(1/IM.sub.weight.loo[cell.common, cell.min.idx]/total * sens.loo[cell.common])
				  } else {
	                  k <- 1
	                  for (j in cell.common) {
	                    subset.min.loo[cell.min.idx] <- (subset.min.loo[cell.min.idx] * k + sens.loo[j])/(k + 1)
	                    k <- k + 1
	                  }
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
                
                # find common item with drug.sub.cell and supersets_small
                cell.common <- which(drug.sub.cell & supersets.small)
                if (length(cell.common) != 0) {
                  # max averaging
				  if (weighted) {
	                  total <- sum(1/IM.sup.weight.loo[cell.common, cell.max.idx]) + 1/IM.sup.weight.loo[drug.maxsens.idx, cell.max.idx]
	  				  superset.max.loo[cell.max.idx] <- 1/IM.sup.weight.loo[drug.maxsens.idx, cell.max.idx]/total*sens.loo[drug.maxsens.idx] + sum(1/IM.sup.weight.loo[cell.common, cell.max.idx]/total*sens.loo[cell.common])
				  } else {
	                  k <- 1
	                  for (j in cell.common) {
	                    superset.max.loo[cell.max.idx] <- (superset.max.loo[cell.max.idx] * k + sens.loo[j])/(k + 1)
	                    k <- k + 1
	                  }
				  }
                }
                
                cell.min.idx <- cell.min.avg[min.eq.identical]
                
                drug.sub.cell <- !is.infinite(IM.subset.loo[, cell.min.idx])
                
                index <- min.idx.loo[cell.min.idx]
                minsens.dec.idx <- drug.idx.loo[index]
                
                # find the subsets of S(index,:) in S that has higher sensitivity
                subsets.high <- IM.superset.loo[, minsens.dec.idx] > subset.min.loo[cell.min.idx]
                # find common item with drug.sub.cell and supersets_small
                cell.common <- which(drug.sub.cell & subsets.high)
                if (length(cell.common) != 0) {
                  # min averaging
				  if (weighted) {
	                  total <- sum(1/IM.sub.weight.loo[cell.common, cell.min.idx]) + 1/IM.sub.weight.loo[index, cell.min.idx]
					  subset.min.loo[cell.min.idx] <- 1/IM.sub.weight.loo[index, cell.min.idx]/total * sens.loo[index] + sum(1/IM.sub.weight.loo[cell.common, cell.min.idx]/total * sens.loo[cell.common])
				  } else {
	                  k <- 1
	                  for (j in cell.common) {
	                    subset.min.loo[cell.min.idx] <- (subset.min.loo[cell.min.idx] * k + sens.loo[j])/(k + 1)
	                    k <- k + 1
	                  }
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
    return(list(efficacy.mat = M, error = err.pred, prediction = pred))
} 
