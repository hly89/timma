#' Generate gray code
#' 
#' A function to generate gray code
#' 
#' @param bit.num an integer to specify the number of bits.
#' @return a vector of the decimal gray code. 
#' @references Dah jyh Guan. (Scientific Note) Generalized Gray Codes with Applications. 1998
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' code<-GetGrayCodeVec(3)

GetGrayCodeVec <- function(bit.num) {
    graycode <- rep(0, 2^bit.num)
    
    # n=1 [0 1]
    graycode[2] <- 1
    if (bit.num == 1) {
        return(graycode)
    } else {
        t <- 2
        # for bit.num lager than 2
        for (i in 2:bit.num) {
            t2 <- t + t
            graycode[(t + 1):t2] <- t + rev(graycode[1:t])
            t <- t2
        }
        return(graycode)
    }
} 
