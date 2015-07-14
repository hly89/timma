#' Generate search space
#' 
#' A function to generate the search space for sffs
#' 
#' @param drug.num an integer to specify the number of drugs
#' @param target.selected a vector to specify the selected target set
#' @param data drug-target interaction data
#' @param sens the drug sensitivity data
#' 
#' @return a list of the following components:
#' \item{space.identical}{search space of identical sets}
#' \item{space.superset}{search space of supersets}
#' \item{space.subset}{search space of subsets}
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' num<-length(tyner_sensitivity[,1])
#' target.selected<-rep(0, dim(tyner_interaction_binary)[2])
#' target.selected[1]<-1
#' space<-GetSearchSpace(num, target.selected, tyner_interaction_binary, tyner_sensitivity[,1])

GetSearchSpace <- function(drug.num, target.selected, data, sens) {
    # parameter 1: drug_num, the number of drugs 
    # parameter 2: target.selected, the current kinase set 
    # parameter 3: data, the complete drug-target-profile data
    # parameter 4: sens, the actual efficacy
    
    # extend one col by 0
    data.extend.zero <- cbind(data[, which(target.selected == 1)], rep(0, drug.num))
    
    
    # extend one col by 1
    data.extend.one <- cbind(data[, which(target.selected == 1)], rep(1, drug.num))
    data.extend <- unique(rbind(data.extend.one, data.extend.zero))
    data.extend.dec <- apply(data.extend, 1, function(x) strtoi(paste(x, collapse = ""), base = 2))
    
    
    dec <- apply(data.extend.zero, 1, function(x) strtoi(paste(x, collapse = ""), base = 2))
    # for identical
    col.num <- length(data.extend.dec)
    identical.idx <- sapply(dec, function(x) which(data.extend.dec == x))
    
    
    
    space.identical <- array(NA, dim = c(drug.num, col.num, 2))
    space.subset <- array(Inf, dim = c(drug.num, col.num, 2))
    space.superset <- array(-Inf, dim = c(drug.num, col.num, 2))
    
    for (i in 1:drug.num) {
        # get the decimal
        space.identical[i, identical.idx[i], 1] <- sens[i]
        
        
        # get the binary set: superset and subset
        bin.set <- GetBinSet(data.extend.zero[i, ])
        # ismember function R version: match
        subset.idx <- data.extend.dec %in% bin.set$subset
        space.subset[i, subset.idx, 1] <- sens[i]
        
        superset.idx <- data.extend.dec %in% bin.set$superset
        space.superset[i, superset.idx, 1] <- sens[i]
        
    }
    
    
    dec <- apply(data.extend.one, 1, function(x) strtoi(paste(x, collapse = ""), base = 2))
    # for identical
    col.num <- length(data.extend.dec)
    identical.idx <- sapply(dec, function(x) which(data.extend.dec == x))
    
    for (i in 1:drug.num) {
        # get the decimal
        space.identical[i, identical.idx[i], 2] <- sens[i]
        
        
        # get the binary set: superset and subset
        bin.set <- GetBinSet(data.extend.one[i, ])
        # ismember function R version: match
        subset.idx <- data.extend.dec %in% bin.set$subset
        space.subset[i, subset.idx, 2] <- sens[i]
        
        superset.idx <- data.extend.dec %in% bin.set$superset
        space.superset[i, superset.idx, 2] <- sens[i]
        
    }
    
    return(list(space.identical = space.identical, space.superset = space.superset, space.subset = space.subset))
} 
