#' Gray code function for matrix indexes
#' 
#' A function to generate gray code used for matrix row and column names
#' 
#' @param bit.num an integer to specify the number of bits
#' 
#' @return a list of the following components:
#' \item{graycode.row}{binary gray code as row names of the predicted sensitivity matrix}
#' \item{graycode.col}{binary gray code as column names of the predicted sensitivity matrix}
#' \item{dec.row}{decimal gray code as row names of the predicted sensitivity matrix}
#' \item{dec.col}{decimal gray code as column names of the predicted sensitivity matrix}
#' 
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' names<-GetGrayCodeMat(3)

GetGrayCodeMat <- function(bit.num) {
    if (bit.num == 1) {
        graycode.row <- c(0, 1)
        graycode.col <- c(0, 1)
        dec.row <- graycode.row
        dec.col <- graycode.col
    } else {
        rows <- floor(bit.num/2)
        cols <- bit.num - rows
        dec.row <- GetGrayCodeVec(rows)
        dec.col <- GetGrayCodeVec(cols)
        graycode.row <- ConvertDecToBin(dec.row, rows)
        graycode.col <- ConvertDecToBin(dec.col, cols)
    }
    return(list(graycode.row = graycode.row, graycode.col = graycode.col, dec.row = dec.row, dec.col = dec.col))
} 
