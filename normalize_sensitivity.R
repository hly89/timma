#' Normalize the drug sensitivity data
#' 
#' A function to normalize the drug sensitivity data to [0,1]
#' 
#' @param IC50 a vector contains the drug sensitivity in the form of IC50.
#' @param method a string to specify the method used to normalize the sensitivity data. If it is
#' "minMax", the sensitivity is scaled by (Max_IC50-IC50)/(Max_IC50-Min_IC50). If it is "logistic",
#' it is scaled by 1/(1+exp(-1/IC50)). If it is "hyperbolic", it is scaled by tanh(1/IC50).
#' @return A vector contains the normalized drug sensitivity data.
#' 
#' @author Jing Tang \email{jing.tang@@helsinki.fi} 
#' @examples 
#' data(tyner_sensitivity)
#' normalizedSensitivity<-NormalizeSensitivity(tyner_sensitivity[,1])
#' 
NormalizeSensitivity<-function(IC50, method="minMax"){
  IC50.normalized <- IC50
  if(method == "minMax"){
    # get the max IC50
    IC50.max <- max(IC50, na.rm=TRUE)
    # get the min IC50
    IC50.min <- min(IC50, na.rm=TRUE)
    for(i in 1:length(IC50)){
      IC50.normalized[i] <- (IC50.max-IC50[i])/(IC50.max-IC50.min)
    }
  }else if(method == "logistic"){
    for(i in 1:length(IC50)){
      IC50.normalized[i] <- 1/(1+exp(-1/IC50[i]))
    }
  }else if(method == "hyperbolic"){
    for(i in 1:length(IC50)){
      IC50.normalized[i] <- tanh(1/IC50[i])
    }
  }else{
    stop("method must be one of minmax/logistic/hyperbolic.")
  }
  return(IC50.normalized)
}