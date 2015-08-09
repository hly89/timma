#' Model selection with sffs for the binary drug-target interaction data
#' 
#' A function to select the most predictive targets with sffs for the binary drug-target interaction data using
#' orignal maximization and minimization rules
#' 
#' @param data drug-target interaction data which is a matrix with drugs as row indexes and targets
#' as column indexes.
#' @param sens a drug sensitivity vector.
#' @param sp an integer to specify the starting point for sequential forward floating search (sffs) search 
#' algorithm to navigate the target set space. By default, sp = 1.
#' @param k.max an integer to sepcify the maximum number of targets that can be selected by the sffs
#' algorithm. By default, k.max = 8. In practice it should not be over than 10 as the number of target combinations will increase exponentially.
#' @param loo a logical value indicating whether to use the leave-one-out cross-validation in the model
#' selection process. By default, loo = True. 
#' @param averaging a parameter to specify which one of the averaging algorithms will be applied 
#' in the model construction. By default, averaging = "one.sided", which is the original model 
#' construction algorithm. When averaging = "two.sided", a modified averaging algorithm will be used.
#' These two variants only differ for the case where the minimization and maximization rules are not
#' simultaneously satisfied. For example, for a queried target set if the supersets but not the subsets
#' can be found in the training data, the one.sided algorithm will take the prediction from the averages
#' on the supersets sensitivities using the minimization rule. The two.sided algorithm, however, will 
#' lower the predicted sensitivity by averaging it with 0, which is the theoretical lower boundary of 
#' the sensitivities that could be obtained in the subsets. 
#' @param weighted a parameter to specify if the similarity between the queried target set and 
#' its subsets/supersets is considered as a weight factor in the averaging. When weighted =T RUE,
#' the similarity is considered as a weight factor such that those related target sets will be
#' weighted more in the final predictions.
#' @param verbosity a boolean value to decide if the information should be displayed. If it is TRUE, the information
#' will be displayed while the model is running. Otherwise, the information will not be displayed. By default, it is
#' FALSE.
#' @return A list containing the following components:
#' \item{timma}{a list contains: the predicted efficacy matrix, prediction error and predicted drug sensitivity}
#' \item{k_sel}{the indexes for selected targets}
#' 
#' @author Jing Tang \email{jing.tang@@helsinki.fi} 
#' @references Tang J, Karhinen L, Xu T, Szwajda A, Yadav B, Wennerberg K, Aittokallio T. 
#' Target inhibition networks: predicting selective combinations of druggable targets to block cancer 
#' survival pathways. PLOS Computational Biology 2013; 9: e1003226.
#' @examples 
#' \dontrun{
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' results <- RunSffs(tyner_interaction_binary, tyner_sensitivity[, 1])
#' }

RunSffs <- function(data, sens, sp = 1, k.max = 8, loo = TRUE, verbosity = FALSE, averaging="one.sided", class = 2, weighted = FALSE, filtering = FALSE, fixed = c()) {	
    # sffs for binary timma model
	
	# if the averaging method parameter is correct
	if (!(averaging %in% c("one.sided", "two.sided"))) {
		stop("Unkown averaging method!")
	}
	
	# if fixed is a target name list, map to targets index
	if (class(fixed) == "character") {
		fixed = match(fixed,colnames(data))
	}
	
    if (!(sp %in% fixed) && class(fixed)=="numeric") {
		# change starting point to the fixed targets index
    	sp = fixed
    } 
        
	
	if (filtering) {
		list.corr <- cor(data, sens)
		list.keep <- which(list.corr > 0)
		data <- data[, list.keep]
		list.len <- length(list.keep)
		# get the new starting point
		# if old sp list is in the filtered target list
		sp.new <- match(sp, list.keep)
		# remove NA from the matching results
		sp.new <- sort(sp.new, na.last = NA)
	    if (length(sp.new) != 0) {
	        target.list <- rep(0, list.len)
	        target.list[sp.new] <- 1
	    } else {
	        target.list <- c(1, rep(0, list.len - 1))
	        cat("The original starting point has beed filtered. \n")
	    }
	} else {
	    # binary vector to indicate if a target is in the cancer-specific target sets
	    target.list <- rep(0, ncol(data))
	    # set the sp target indicator to be 1
	    target.list[sp] <- 1
	}
    
    # get the number of the drugs
    drug.num <- nrow(data)
    # get the number of the targets
    target.num <- ncol(data)
    
    
    
    # the criterion function J(X_k)
    J <- rep(0, target.num)
    # the set X_k
    X <- matrix(0, target.num, target.num)
    
    run <- 0
    step <- 0
    change <- 1
    err.global.best <- 100
    
    while (change) {
        run <- run + 1
        target.list.tmp <- target.list
        
        while (sum(target.list.tmp) <= k.max) {
            err <- rep(Inf, target.num)
            step <- sum(target.list.tmp)
            data.temp <- data[, which(target.list.tmp == 1)]
			timma.results <- RunTimma(data.temp, sens, loo, averaging, class, weighted)
            
            # the score for X_k, the lower score is better
            J[step] <- mean(timma.results$error)
            # if target.list.tmp has been visited
            if (all(X[, step] == target.list.tmp)) {
                if(verbosity){
                  cat("current target list has been visted! \n")
                }
                
                tmp <- which(J == min(J[which(J > 0)]))
                target.list.tmp <- X[, tmp[1]]
                change <- 0
                break
            } else {
                X[, step] <- target.list.tmp
            }
            if(verbosity){
              cat("Number of selected targets: ", step, "MAE = ", J[step], "\n")
            }
            
            
            if (step > 1 && J[step] == J[step - 1]) {
                # if no improvement
                target.list <- rep(0, target.num)
                target.list[ind1] <- 1
                if (J[step] < err.global.best) {
                  err.global.best = J[step]
                  change <- 1
                  if(verbosity){
                    cat("Best MAE:", err.global.best, "\n")
                  }
                  
                } else {
                  change <- 0
                }
                break
            }
            
            # step 1: Inclusion
			if (class == 2) {
				space <- GetSearchSpace(drug.num, target.list.tmp, data, sens)
			}
            
            for (i in 1:target.num) {
                if (target.list.tmp[i] == 0) {
                  # including new target i
                  err[i] <- 0
				  if (class == 2) {
					  err1 <- GetPredErr(data[, i], space, sens, loo, averaging)
				  } else {
					  set.add.one <- target.list.tmp
					  set.add.one[i] <- 1
					  err1 <- RunTimma(data[, which(set.add.one == 1)], sens, loo, averaging, class, weighted)$error
				  }				  
                  err[i] <- mean(err1)
                }
            }
            
            # select X_{k+1} the most significant feature with respect to X_k
            dummy <- min(err, na.rm = TRUE)
            ind1 <- which.min(err)
            target.list.tmp[ind1] <- 1
            if(verbosity){
              cat("Including target #", ind1, "MAE = ", dummy, "\n")
            }
            
            
            # step2: Conditional exclusion find the least significant feature in the set X_{k+1}
            err <- rep(Inf, target.num)
            for (i in 1:target.num) {
                target.list.current <- target.list.tmp  # the current target set
                if (target.list.tmp[i] == 1) {
                  # if removing the i-th target
                  err[i] <- 0
                  target.list.current[i] <- 0
                  data.temp <- data[, which(target.list.current == 1)]  # which(target.list.current==1)
				  timma.results <- RunTimma(data.temp, sens, loo, averaging, class, weighted)
                  
                  err[i] <- mean(timma.results$error)
                }
            }
            err.worst <- min(err)
            worst.idx <- which.min(err)
            if (err[ind1] != err.worst && err.worst < J[step] && !(worst.idx %in% fixed)) {
                # the new added target is not the least significant feature
                target.list.tmp[worst.idx] <- 0  # exclude the least significant feature
                if(verbosity){
                  cat("Excluding target #", worst.idx, "MAE = ", err.worst, "\n")
                }
                
            }
            # step 3: Continuation of conditional exclusion
            while (sum(target.list.tmp) > 2) {
                err <- rep(Inf, target.num)
                for (i in 1:target.num) {
                  target.list.current <- target.list.tmp
                  if (target.list.tmp[i] == 1) {
                    err[i] <- 0
                    target.list.current[i] <- 0
					          timma.results <- RunTimma(data.temp, sens, loo, averaging, class, weighted)
                    
                    err[i] <- mean(timma.results$error)
                    
                  }
                }
                dummy <- min(err, na.rm = TRUE)
                worst.idx <- which.min(err)
                if (dummy < J[length(which(target.list.tmp == 1)) - 1] && !(worst.idx %in% fixed)) {
                  # if better result found !!J[step-1]!!
                  target.list.tmp[worst.idx] <- 0
                  if(verbosity){
                    cat("Continuing excluding target #", worst.idx, "MAE = ", dummy, "\n")
                  }
                  
                } else {
                  break
                  if(verbosity){
                    cat("break \n")
                  }
                  
                }
            }
            change <- 0
        }
    }
    
    temp <- which(J == min(J[which(J > 0)]))
    if (length(temp) > 1) {
        target.list.tmp <- X[, temp[2]]
    } else {
        target.list.tmp <- X[, temp[1]]
    }
    
    target.selected <- which(target.list.tmp == 1)
    #data.filtered <- unique(data[, target.selected], MARGIN = 2)
    data.filtered <- as.matrix(data[, target.selected])
	
	
	if (class == 2) {
		timma <- RunTimmaFullBin(data.filtered, sens, loo, averaging)
	} else {
		timma <- RunTimma(data.filtered, sens, loo, averaging, class, weighted)
	}
    
    #target.selected <- match(dimnames(data.filtered)[[2]],dimnames(data)[[2]])
    target.selected <- dimnames(data)[[2]][target.selected]
    return(list(timma = timma, target.selected = target.selected, data.selected = data.filtered))
} 
