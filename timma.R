#' Main function for the timma package
#' 
#' Target inhibition inference using maximization and minimization averaging
#'
#' @param x a drug-target interaction matrix. Row names are drug names and column names are target names.
#' @param y a normalized drug sensitivity vector.
#' @param sp an integer to specify the starting point for the sffs search algorithm. The number cannot be larger than the total number of targets in the drug-target interaction data. 
#' By default, the starting point is the first target, namely, sp = 1.
#' @param k.max an integer to specify the maximal number of targets that can be selected by the sffs algorithm. In practice it is advised to keep it under 10 
#' as the number of sensitivities to be predicted will increase exponentially. By default, k.max = 5.
#' @param filtering a logical parameter to determine whether the targets should be filtered before the model selection. 
#' By default, the value is FALSE, meaning that all the available targets will be considered in the model selection. 
#' If the value is TRUE, those targets that are negatively correlated with the drug sensitivity data will be removed.
#' @param class an integer to specify the number of classes in the drug-target interaction data. For a binary drug-target 
#' interaction data, class = 2. For a multi-class drug-target interaction data, class should be the number of classes.
#' @param averaging a parameter to specify which one of the averaging algorithms will be applied in the model construction. 
#' By default, averaging = "one.sided", which is the original model construction algorithm. When averaging = "two.sided", 
#' a modified averaging algorithm will be used. These two variants only differ for the case where the minimization and 
#' maximization rules are not simultaneously satisfied. For example, for a queried target set if the supersets but not the 
#' subsets can be found in the training data, the one.sided algorithm will take the prediction from the averages on the 
#' supersets sensitivities using the minimization rule. The two.sided algorithm, however, will lower the predicted 
#' sensitivity by averaging it with 0, which is the theoretical lower boundary of the sensitivities that could be 
#' obtained in the subsets.
#' @param weighted When averaging = "weighted", the similarity between the queried target set and 
#' its subsets/supersets is considered as a weight factor in the averaging, such that those related target 
#' sets will be more weighted in the final predictions.
#' @param verbosity a boolean value to decide if the information should be displayed. If it is TRUE, the information
#' will be displayed while the model is running. Otherwise, the information will not be displayed. By default, it is
#' FALSE.
#' @param use When use = "observed", the true drug sensitivity data will be used for drawing target inhibition network.
#' When use = "predicted", the predicted drug sensitivity data will be used for drawing target inhibition network.
#' @return an R image of the input and output data. 
#' @author Jing Tang \email{jing.tang@@helsinki.fi} 
#' @references Tang J, Karhinen L, Xu T, Szwajda A, Yadav B, Wennerberg K, Aittokallio T. 
#' Target inhibition networks: predicting selective combinations of druggable targets to block cancer 
#' survival pathways. PLOS Computational Biology 2013; 9: e1003226.
#' 
#' @examples 
#' \dontrun{
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' median_sensitivity<-tyner_sensitivity[, 1]
#' results<-timma(tyner_interaction_binary, median_sensitivity)
#' }

timma <- function(x, y, sp = 1, k.max = 5, filtering = FALSE, class = 2, averaging = "one.sided", weighted = FALSE, verbosity = FALSE, use = "observed", fixed = c()) {
  
  # get drug names
  drug.names <- dimnames(x)[[1]]
  
  # get kinase names
  kinase.names <- dimnames(x)[[2]]
  
  # preprocessing the data: merge the same colomns and delete all zeros colomns
  list.del <- c()
  name.merged <- c() # the same length as the number of columns

  list.ones <- rep(1, nrow(x))              
  for (i in 1:ncol(x)){
    if (i %in% list.del) {
      # if the i-th colomn is in the delete list
      name.merged <- c(name.merged, kinase.names[i])
    } else {
      # detect the col with all zeros
      if (list.ones %*% x[,i] == 0) {
        list.del <- c(list.del, i)
        name.merged <- c(name.merged, kinase.names[i])
      } else {
        col.same <- FindSameCol(x, x[,i])
        name.merged <- c(name.merged, paste(kinase.names[which(col.same == 1)], collapse = ";"))
        list.del <- c(list.del,which(col.same==1)[-1])
      }
    }
    
  }

  x<-x[,-list.del]
  name.merged <- name.merged[-list.del]
  dimnames(x)[[2]] <- name.merged
  kinase.names <- name.merged
  
  cat("----------------Start Running TIMMA----------------------------- \n") 
  float<-RunSffs(x, y, sp, k.max, loo = TRUE, verbosity, averaging, class, weighted, filtering, fixed)  
  cat("----------------Complete running TIMMA model-------------------- \n")
  
  err <- float$timma$error
  error <- mean(err)
  RMSE <- sqrt(mean(err^2))
  RMSE.baseline <- sqrt(mean((y - median(y))^2))
  
  # write selected targets file
  #x <- data.frame(x)
  target.selected <- float$target.selected
  #target.selected.names <- FindSameSet(x, target.selected, kinase.names)
  #data.selected <- x[, target.selected]
  data.selected <- data.frame(float$data.selected)
  target.selected.names <- dimnames(float$data.selected)[[2]]
  
  if (k.max > 1) {  
    target.selected.res <- cbind(data.selected, y, float$timma$prediction)
    dimnames(data.selected)[[2]] <- target.selected.names
    dimnames(target.selected.res)[[2]] <- c(target.selected.names, "sensitivity", "LOO sensitivity")
    write.table(target.selected.res, file = "selectedTargets.csv", sep = ",", col.names = NA)
  } else {
    data.selected <- rbind(target.selected.names, data.selected)
  }
  cat("----------------Saving selectedTargets.csv---------------------- \n")
  
  # write Timma
  if (class == 2) {
    graycode <- GetGrayCodeMat(length(target.selected))
    graycode.names <- GetGrayCodeNames(length(target.selected), target.selected.names, graycode$graycode.row, graycode$graycode.col)
    nr <- graycode.names$graycode.names.row
    nc <- t(graycode.names$graycode.names.col)
    timma.row <- nrow(nr) + nrow(nc)
    timma.col <- ncol(nr) + ncol(nc)
    timma <- array("", dim = c(timma.row, timma.col))
    timma[(nrow(nc) + 1):timma.row, 1:ncol(nr)] <- nr
    timma[1:nrow(nc), (ncol(nr) + 1):timma.col] <- nc
    timma[(nrow(nc) + 1):timma.row, (ncol(nr) + 1):timma.col] <- float$timma$efficacy.mat
    write.table(timma, file = "predictedSensitivities.csv", sep = ",", col.names = FALSE, row.names = FALSE)
    cat("----------------Saving predictedSensitivities.csv--------------- \n")
    
    target.comb.rank <- RankTargetComb(data.selected, timma)
    write.table(target.comb.rank, file = "predictedTargetScoring.csv",sep=",",row.names = FALSE)
    cat("----------------Saving predictedTargetScoring.csv--------------------- \n")
    
    drug.comb.rank <- RankDrugComb(data.selected, timma, y)
    write.table(drug.comb.rank, file = "predictedDrugScoring.csv", sep = ",", row.names = FALSE)
    cat("----------------Saving predictedDrugScoring.csv--------------------- \n")
    
    # write the R image
    save(x,data.selected, timma, y, file="result.RData")
    cat("----------------Saving result.RData---------------------- \n")
    
    #prediction.loo <- float$timma$prediction
    if(use == "observed"){
      prediction.loo <- y
    }else if(use == "predicted"){
      prediction.loo <- float$timma$prediction
    }else{
      stop("Parameter use must be either observed or predicted")
    }
    
    one<-which(prediction.loo>0.5)
    zero<-which(prediction.loo<=0.5)
    SENS<-prediction.loo
    SENS[one]<-1
    SENS[zero]<-0
   
    draw.data<-cbind(data.selected, SENS)
    DrawNetwork(draw.data)
    cat("----------------Saving targetInhibitionNetwork.pdf-------------- \n")
    cat("----------------Saving targetInhibitionNetwork.nnf-------------- \n")
  }
  
  dir.current <- getwd()
  cat("Analysis finished. All the results are saved in", dir.current)

  return(list(data = x, data.selected = data.selected, prediction = float, sensitivity = y))
} 
