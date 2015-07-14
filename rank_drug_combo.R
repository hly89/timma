#' Generate the list of ranked drug combinations
#' 
#' A function to provide a list of drug combinations ranked by their synergy scores
#' 
#' @param data.selected the selected drug-target interaction data
#' @param predicted_matrix the predicted efficacy matrix
#' @param sens the drug sensitivity vector.
#' 
#' @return a matrix contains the information about the list of drug combinations ranked by their synergy scores. 
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
#' drug_combo_rank<-RankDrugComb(data.selected, timma, tyner_sensitivity[, 1])
#' }

RankDrugComb <- function(data.selected, efficacy.mat, sens) {
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
    
    # identical targets separated by '/'.  Replace ';' with '/' in timma1.csv Replace ';' with '/' in
    # selectedKinasse1.csv empty entries denoted by NA. Replace '-' with NA in timma1.csv
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
    colnames(inhibition.single) = c("timma1", "timma2", "timma.add", "timma.multi", "timma.highest")
    for (i in 1:dim(data.pairwise)[1]) {
        pair = data.pairwise$Gene[i]
        pair.name = unlist(strsplit(as.character(pair), ";"))
        # index = lapply(pair.name, function(x) grep(x, data.single$Gene, fixed = TRUE))
        index = lapply(pair.name, function(x) which(data.single$Gene==x))
        timma1 = data.single$timma[unlist(index[1])]
        timma2 = data.single$timma[unlist(index[2])]
        timma.add = data.pairwise$timma[i] - timma1 - timma2
        timma.multi = data.pairwise$timma[i] - timma1 * timma2
        # timma.multi2 = timma.add+timma1*timma2
        timma.highest = data.pairwise$timma[i] - max(timma1, timma2)
        inhibition.single[i, ] = c(timma1, timma2, timma.add, timma.multi, timma.highest)
    }
    data.pairwise = cbind(data.pairwise, inhibition.single)
    
    # find the drugs that are corresponding to the targets
    names(data.selected)[1] = "Drugs"
    colnames(data.selected) = gsub(";", "/", colnames(data.selected), fixed = TRUE)
    drug.table = list(0)
    for (i in 1:target.num) {
        target = colnames(data.selected)[i + 1]
        index = which(data.selected[, i + 1] == 1)
        drugs = data.selected$Drugs[index]
        # drug.table[[i]] = expand.grid(target,drugs)
        drug.table[[i]] = expand.grid(target, index)
    }
    drug.table = do.call("rbind", drug.table)
    
    
    index = drug.table[, 2]
    drug.table = cbind(drug.table, data.selected$Drugs[index], sens[index])
    colnames(drug.table) = c("Gene", "Index", "Drug", "Sens")
    
    drug.score = list(0)
    for (i in 1:dim(data.pairwise)[1]) {
        pair = data.pairwise$Gene[i]
        pair.str = strsplit(as.character(pair), ";")[[1]]
        gene1 = pair.str[1]
        gene2 = pair.str[2]
        # drug1.idx = grep(gene1, drug.table$Gene, fixed = TRUE)  # indices of drugs that are targeting gene1
        drug1.idx = which(drug.table$Gene==gene1)
        # drug2.idx = grep(gene2, drug.table$Gene, fixed = TRUE)  # indices of drugs that are targeting gene2
        drug2.idx = which(drug.table$Gene==gene2)
        # drugs that are targeting gene1
        drug1 = drug.table$Drug[drug1.idx]
        # drugs that are targeting gene2
        drug2 = drug.table$Drug[drug2.idx]
        n = expand.grid(drug1, drug2)
        n = t(apply(n, 1, sort))  # sort drug names alphabetically
        order.new = t(apply(expand.grid(drug1, drug2), 1, order))
        
        tmp = mat.or.vec(dim(n)[1], 8)
        for (j in 1:dim(n)[1]) {
            tmp[j, 1] = as.character(n[j, 1])
            tmp[j, 2] = as.character(n[j, 2])
            tmp[j, 3] = drug.table$Sens[which(drug.table$Drug == n[j, 1])[1]]
            tmp[j, 4] = drug.table$Sens[which(drug.table$Drug == n[j, 2])[1]]
            tmp[j, 5] = data.pairwise[i, 3 + order.new[j, 1]]
            tmp[j, 6] = data.pairwise[i, 3 + order.new[j, 2]]
            tmp[j, 7] = pair.str[order.new[j, 1]]
            tmp[j, 8] = pair.str[order.new[j, 2]]
        }
        
        tmp = cbind(tmp, data.pairwise[rep(i, each = dim(n)[1]), ])
        
        drug.score[[i]] = tmp
    }
    
    drug.score.mat = do.call("rbind", drug.score)
    colnames(drug.score.mat)[1:8] = c("Drug1", "Drug2", "Drug1.sens", "Drug2.sens", "timma1.s", "timma2.s", 
        "Target1", "Target2")
    drug.score.mat$Drug1.sens = as.numeric(as.character(drug.score.mat$Drug1.sens))
    drug.score.mat$Drug2.sens = as.numeric(as.character(drug.score.mat$Drug2.sens))
    drug.score.mat$Drug1 = as.character(drug.score.mat$Drug1)
    drug.score.mat$Drug2 = as.character(drug.score.mat$Drug2)
    drug.score.mat$timma1.s = as.numeric(as.character(drug.score.mat$timma1.s))  # sorted timma1
    drug.score.mat$timma2.s = as.numeric(as.character(drug.score.mat$timma2.s))  # sorted timma2
    
    drug.combinations = c(0)
    for (i in 1:dim(drug.score.mat)[1]) {
        drug.combinations[i] = paste(sort(c(drug.score.mat$Drug1[i], drug.score.mat$Drug2[i])), collapse = " ")
    }
    drug.score.mat$drug.combinations = as.factor(drug.combinations)
    comp1 = aggregate(drug.score.mat[, c(1:2)], list(drug.combinations), unique)
    comp2 = aggregate(drug.score.mat[, c(3:6, 10, 12:16)], list(drug.combinations), mean)
    comp3 = aggregate(drug.score.mat[, 7], list(drug.combinations), function(i) paste(unique(i), collapse = ","))
    comp4 = aggregate(drug.score.mat[, 8], list(drug.combinations), function(i) paste(unique(i), collapse = ","))
    colnames(comp3)[2] = "Target1"
    colnames(comp4)[2] = "Target2"
    drugcomb.score = cbind(comp1, comp2, comp3, comp4)
    drugcomb.score = drugcomb.score[,c(2,3,9,5,6,12,13,14,16,18)]
    colnames(drugcomb.score) = c("Drug1", "Drug2", "Predicted Sensitivity", "Drug1 Sensitivity", "Drug2 Sensitivity", "Synergy.add", "Synergy.multi", "Synergy.highest", 
    "Target1", "Target2")
    return(drugcomb.score)
}
 
