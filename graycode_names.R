#' Names for the predicted sensitivity matrix
#' 
#' A function to make the target names in the format of gray code for the predected sensitivity matrix
#' 
#' @param target.num an integer to specify the number of targets
#' @param names a vector of the names of the targets
#' @param graycode.row the gray code as row indexes. It can be returned by \code{\link{graycode3}}.
#' @param graycode.col the gray code as column indexes. It can be returned by \code{\link{graycode3}}.
#' 
#' @return a list of the following components:
#' \item{graycode.names.row}{the gray code format target names as row names.}
#' \item{graycode.names.col}{the gray code format target names as column names.}
#' 
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' \dontrun{
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' k_select<-sffsBinary(tyner_interaction_binary, tyner_sensitivity[, 1])$k_sel
#' gc_timma<-GetGrayCodeMat(length(k_select))
#' select_kinase_names<-dimnames(tyner_interaction_binary)[[2]][k_select]
#' gc_names<-GetGrayCodeNames(length(k_select), select_kinase_names, gc_timma$graycode.row, gc_timma$graycode.col)
#' }
GetGrayCodeNames <- function(target.num, names, graycode.row, graycode.col) {
    # parameter target.num: the number of how many selected kinases    
    # parameter names: the names for selected kinases
    row.dim <- dim(graycode.row)
    col.dim <- dim(graycode.col)
    
    # get the kinase names for rows
    row.names <- names[1:row.dim[2]]
    # get the kinase names for columns
    col.names <- names[row.dim[2] + 1:target.num]
    
    graycode.names.row <- array(NA, dim = row.dim)
    graycode.names.col <- array(NA, dim = col.dim)
    for (i in 1:row.dim[1]) {
        for (j in 1:row.dim[2]) {
            if (graycode.row[i, j] == 0) {
                graycode.names.row[i, j] <- "-"
            } else {
                graycode.names.row[i, j] <- row.names[j]
            }
        }
    }
    for (i in 1:col.dim[1]) {
        for (j in 1:col.dim[2]) {
            if (graycode.col[i, j] == 0) {
                graycode.names.col[i, j] <- "-"
            } else {
                graycode.names.col[i, j] <- col.names[j]
            }
        }
    }
    return(list(graycode.names.row = graycode.names.row, graycode.names.col = graycode.names.col))
} 
