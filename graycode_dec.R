#' Graycode Function
#' 
#' A function to generate decimal graycode
#' 
#' @param a the number of targets
#' 
#' @return A list contains the following components:
#' \item{rows}{the number of rows}
#' \item{cols}{the number of columns}
#' \item{dec}{the decimal graycode results}
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @references Dah jyh Guan. (Scientific Note) Generalized Gray Codes with Applications. 1998
#' @examples 
#' code<-GetGrayCodeDec(5)
GetGrayCodeDec <- function(a) {
    # parameter1 a: the number of targets
    
    rows <- 2^floor(a/2)
    cols <- 2^(a - floor(a/2))
    # get the decimal
    dec <- GetGrayCodeVec(a)
    # reshape the index to a matrix
    dim(dec) <- c(cols, rows)
    dec <- t(dec)
    # flip every second row of the matrix
    dec.nrow <- nrow(dec)
    if (dec.nrow == 1) {
        return(list(rows, cols, dec))
    } else if (dec.nrow == 2) {
        dec[2, ] <- rev(dec[2, ])
        return(list(rows, cols, dec))
    } else {
        dec[seq(2, dec.nrow, 2), ] <- t(apply(dec[seq(2, dec.nrow, 2), ], 1, rev))  #does not work for nrow<=2
        return(list(rows, cols, dec))
    }
    
} 
