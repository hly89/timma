#' Search for supersets and subsets
#' 
#' A function for searching the supersets and subsets of the binary drug-target interaction data.
#' 
#' @param profile_data the binary drug-target interaction matrix with row indexes as drugs and column
#' indexes as targets.
#' @return A list contains the following components:
#' \item{superset}{all the possible supersets of the input drug-target interaction data}
#' \item{subset}{all the possible subsets of the input drug-target interaction data}
#' 
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' data(tyner_interaction_binary)
#' sets<-GetBinSet(tyner_interaction_binary[1, 1:3])
#' 
GetBinSet <- function(data) {
    # return all possible supersets and subsets of input drug target profile data as binary set 
    # e.g. data=[1 0 0 1], subset=[0 0 0 1], superset=[1 1 1 1] return multiple paramters
    
    
    # [1 1 1 1]: 2^4-1 [1 0 0 1]: 2^4-1-2^2-2^1
    profile.len <- length(data)
    zeros.dec <- 2^(profile.len - which(data == 0))
    ones.dec <- 2^(profile.len - which(data == 1))
    data.dec <- 2^profile.len - 1 - sum(zeros.dec)
    
    # gray code for zeros 
    zeros.len <- length(zeros.dec)
	# check zeros.dec is empty or not
    if (zeros.len == 0) {
        superset = data.dec
        
    } else {
        graycode.zeros <- ConvertDecToBin(GetGrayCodeVec(zeros.len), zeros.len)
        graycode.zeros.len <- nrow(graycode.zeros)
        zeros.dec.bs <- t(matrix(rep(zeros.dec, graycode.zeros.len), ncol = graycode.zeros.len))
        superset <- rowSums(zeros.dec.bs * graycode.zeros) + data.dec
    }
    
    
    
    # gray code for ones 
    ones.len <- length(ones.dec)
	# check ones.dec is empty or not
    if (ones.len == 0) {
        # set subset NULL
        subset <- vector("numeric")
    } else {
        graycode.ones <- ConvertDecToBin(GetGrayCodeVec(ones.len), ones.len)
        # remove first line
        graycode.ones <- graycode.ones[-1, ]
        graycode.ones.len <- nrow(graycode.ones)
        if (is.null(graycode.ones.len)) 
            graycode.ones.len <- 1
        ones.dec.bs <- t(matrix(rep(ones.dec, graycode.ones.len), ncol = graycode.ones.len))
        subset <- data.dec - rowSums(ones.dec.bs * graycode.ones)
    }
    
    return(list(superset = superset, subset = subset))
}

 
