#' Generate the list of ranked target combinations
#' 
#' A function to provide a list of target combiantions ranked by their predicted synergy scores
#' 
#' @param data.selected the drug-target interaction profile for the selected targets
#' @param efficacy.mat the predicted efficacy matrix
#' 
#' @return a matrix containing the list of target combinations
#' @author Jing Tang \email{jing.tang@@helsinki.fi} 
#' @references Tang J, Karhinen L, Xu T, Szwajda A, Yadav B, Wennerberg K, Aittokallio T. 
#' Target inhibition networks: predicting selective combinations of druggable targets to block cancer 
#' survival pathways. PLOS Computational Biology 2013; 9: e1003226.
#' @examples
#' \dontrun{
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' float<-sffsBinary(tyner_interaction_binary, tyner_sensitivity[, 1], max_k = 8)
#' k_select<-float$k_sel
#' x<-data.frame(tyner_interaction_binary)
#' kinase_names <- dimnames(x)[[2]]
#' select_kinase_names <- findSameSet(x, k_select, kinase_names)
#' gc_timma <- graycode3(length(k_select))
#' gc_names <- graycodeNames(length(k_select), select_kinase_names, gc_timma$gc_row, gc_timma$gc_col)
#' nr <- gc_names$nr
#' nc <- t(gc_names$nc)
#' timma_row <- nrow(nr) + nrow(nc)
#' timma_col <- ncol(nr) + ncol(nc)
#' timma <- array("", dim = c(timma_row, timma_col))
#' timma[(nrow(nc) + 1):timma_row, 1:ncol(nr)] <- nr
#' timma[1:nrow(nc), (ncol(nr) + 1):timma_col] <- nc
#' timma[(nrow(nc) + 1):timma_row, (ncol(nr) + 1):timma_col] <- float$timma$dummy
#' data.selected<-data.frame(tyner_interaction_binary)[, k_select]
#' target_combo_rank<-RankTargetComb(data.selected, timma)
#' }


RankTargetComb <- function(data.selected, efficacy.mat) {
  data.selected <- cbind(dimnames(data.selected)[[1]], data.selected)
  drug.num = dim(data.selected)[1]
  target.num = dim(data.selected)[2] - 1
  
  timma <- efficacy.mat
  len1=(target.num-floor(target.num/2))
  len2=2^(floor(target.num/2))
  len3=floor(target.num/2)
  len4=2^(target.num-floor(target.num/2))
  
  if(len3==1) target.row = matrix(timma[c((1 + len1):(len2 + len1)),1],ncol=1) else target.row = apply(timma[c((1 + len1):(len2 + len1)), c(1:len3)], 2, as.character)
  if(len1==1) target.col = matrix(timma[1, c((1 + len3):(len4 + len3))],nrow=1) else target.col = apply(timma[c(1:len1), c((1 + len3):(len4 + len3))], 2, as.character)
  
  efficacy.mat = apply(timma[c((1 + len1):(len2 + len1)), c((1 + len3):(len4 + len3))], 2, as.numeric)
  # efficacy.mat = efficacy.mat/max(efficacy.mat)  # normalize the efficacy into range [0 1]
  efficacy.vec = matrix(efficacy.mat, length(c(efficacy.mat)), 1)
  
  # identical targets separated by '/'.  Replace ';' with '/' in target1.sens.csv Replace ';' with '/' in
  # selectedKinasse1.csv empty entries denoted by NA. Replace '-' with NA in target1.sens.csv
  target.row[which(target.row == "-")] = NA
  target.col[which(target.col == "-")] = NA
  target.row = gsub(";", "/", target.row)
  target.col = gsub(";", "/", target.col)
  
  names.mat = matrix(0,dim(target.row)[1], dim(target.col)[2])  # target combinations
  target.num.comb = matrix(0,dim(target.row)[1], dim(target.col)[2])  # number of target nodes in the combination
  for (i in 1:dim(target.row)[1]) {
    for (j in 1:dim(target.col)[2]) {
      str1 = paste(paste(target.row[i, ], collapse = " "), paste(target.col[, j], collapse = " "), collapse = "", 
                   sep = " ")
      # str2 = gsub('NA','',str1)
      str2 = gsub("\\<NA\\>", "", str1)  # exact match
      str3 = gsub("^\\s+|\\s+$", "", str2)  # remove leading and trailing space
      str4 = strsplit(str3, "\\s+")
      names.mat[i, j] = paste(unlist(str4), collapse = ";")
      target.num.comb[i, j] = length(str4[[1]])
    }
  }
  names.vec = matrix(names.mat, length(c(names.mat)), 1)
  number.vec = matrix(target.num.comb, length(c(target.num.comb)), 1)
  target.info = data.frame(cbind(names.vec, efficacy.vec, number.vec))
  colnames(target.info) = c("Gene", "timma", "Number")
  target.info$timma = as.numeric(as.character(target.info$timma))
  
  # deal with only single and pairwise targets
  data.single = target.info[which(target.info$Number == 1), ]
  data.pairwise = target.info[which(target.info$Number == 2), ]
  
  
  inhibition.single = matrix(0,dim(data.pairwise)[1], 5)
  colnames(inhibition.single) = c("target1.sens", "target2.sens", "synergy.add", "synergy.multi", "synergy.highest")
  for (i in 1:dim(data.pairwise)[1]) {
    pair = data.pairwise$Gene[i]
    pair.name = unlist(strsplit(as.character(pair), ";"))
    index = lapply(pair.name, function(x) which(data.single$Gene==x))
    target1.sens = data.single$timma[unlist(index[1])]
    target2.sens = data.single$timma[unlist(index[2])]
    synergy.add = data.pairwise$timma[i] - target1.sens - target2.sens
    synergy.multi = data.pairwise$timma[i] - target1.sens * target2.sens
    # synergy.multi2 = synergy.add+target1.sens*target2.sens
    synergy.highest = data.pairwise$timma[i] - max(target1.sens, target2.sens)
    inhibition.single[i, ] = c(target1.sens, target2.sens, synergy.add, synergy.multi, synergy.highest)
  }
  data.pairwise = cbind(data.pairwise, inhibition.single)
  colnames(data.pairwise) = c("Target.pair","Sensitivity","Size","Target1.sen","Target2.sen","Synergy.add","Synergy.multi","Synergy.highest")
  
  return(data.pairwise)
}