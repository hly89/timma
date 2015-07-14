#' Binary set for multiclass data
#' 
#' A function to get the supersets and subsets for multiclass data
#' 
#' @param input a vector of multiclass data
#' @param data a matrix of multiclass data as training data
#' @return a list of the following components:
#' \item{superset}{the supersets of the input data from the training data}
#' \item{subset}{the subsets of the input data from the training data}
#' 
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' data(tyner_interaction_multiclass)
#' sets<-GetBinSetMulti(tyner_interaction_multiclass[1,], tyner_interaction_multiclass)
GetBinSetMulti <- function(input, data, weighted = FALSE){
  #example:
  #input<-c(0,1,2,3)
  #data<-sample(0:5,40, T)
  #data<-array(data,dim=c(10,4))
  input1<-matrix(rep(input, nrow(data)), ncol=ncol(data), byrow=TRUE)
  res<-input1-data
  drugs<-c(1:nrow(data))
  val.neg<-which(res<0, arr.ind=TRUE)
  row.neg<-unique(val.neg[,1])
 
  zeros<-rowSums(abs(res))
  identical<-which(zeros==0)
  # subset no negative values
  sub<-setdiff(drugs, row.neg)
  sub<-setdiff(sub, intersect(sub, identical))
  if (weighted) {
	  # weight for every subset to the input data
	  if(length(sub)!=0){
	    weight.sub<-rowSums(matrix(res[sub,], ncol=ncol(data)))
	  }else{
		  weight.sub<-0
	  }
  }
  
  
  val.pos<-which(res>0, arr.ind=TRUE)
  row.pos<-unique(val.pos[,1])
  # superset no positive values
  sup<-setdiff(drugs, row.pos)
  if (weighted) {
	  # weight for every superset to the input data
	  if(length(sup)!=0){
	    weight.sup<-abs(rowSums(matrix(res[sup,], ncol=ncol(data))))
	  }else{weight.sup<-0}
  }
  
  if (weighted) {
  	return(list(subset=sub, superset=sup, subw=weight.sub, supw=weight.sup))
  } else {
  	return(list(subset=sub, superset=sup))
  }

  
}