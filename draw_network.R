#' Draw graph function
#' 
#' A function to draw the target inhibition network.
#' 
#' @param data a data frame combining drug-target interaction data with drug sensitivity. The column names
#' must be upper case.
#' @return An image in both pdf and nnf format of the estimated target inhibition network.
#' @author Jing Tang \email{jing.tang@@helsinki.fi}
#' @references Tang J, Karhinen L, Xu T, Szwajda A, Yadav B, Wennerberg K, Aittokallio T. 
#' Target inhibition networks: predicting selective combinations of druggable targets to block cancer 
#' survival pathways. PLOS Computational Biology 2013; 9: e1003226.
#' @examples 
#' \dontrun{
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' y <- tyner_sensitivity[, 1]
#' results <- RunSffs(tyner_interaction_binary, y)
#' x < -data.frame(results$data.selected)
#' #binarize the sensitivity data
#' one <- which(y > 0.5)
#' zero <- which(y <= 0.5)
#' SENS <- y
#' SENS[one] <- 1
#' SENS[zero] <- 0
#' data <- cbind(x, SENS)
#' DrawNetwork(data)
#' }

DrawNetwork <- function(data) {
    
    # column names must be upper case and contains only letters!
    col.names.actual <- dimnames(data)[[2]]
    col.names <- unlist(lapply(seq(1:length(col.names.actual)),function(i) paste(LETTERS[i],LETTERS[i],sep="")))
    colnames(data) <- col.names
    dimnames(data)[[2]][length(col.names)] <- "SENS"
    
    boolean <- eqmcc(data, outcome = "SENS",row.dom = F,omit = 1)
    #boolean$essential
    components.num <- length(boolean$essential)
    expressions <- sapply(boolean$essential, function(x) strsplit(x, "\\*"))
    
    expressions.compact <- list()
    # only treat the true conditions
    is.upper <- "[A-Z]"
    num <- 1
    for (i in 1:components.num) {
        
        result <- grepl(pattern = is.upper, x = expressions[[i]])
        if(any(result)){
          expressions.compact[[num]] <- unlist(lapply(expressions[[i]][which(result == TRUE)],
                                                    function(x) unlist(strsplit(x, ".", fixed=TRUE))[1]))
          num <- num + 1
        }
        
    }
    
    # remove duplicated compact expressions
    expressions.compact <- expressions.compact[which(duplicated(expressions.compact)==F)]
    
    # change the labels back to the target names
    expressions.compact <- lapply(expressions.compact,function(i) col.names.actual[grep(paste(i,collapse="|"),col.names)])
    
    components.num <- length(expressions.compact)
    components.height <- lapply(expressions.compact, length)
    
    
    seg <- 1 # the width for one component
    seg.in <- 0.2 # the width of the small line within the components
    seg.out <-0.2 # the width between the components
    scale <- mean(unlist(lapply(expressions.compact,function(x) max(nchar(x))/3)))
    
    figure.height <- max(unlist(components.height))
    figure.width <- components.num * (seg*scale + seg.out)
    
    # drawing
    dummy <- 0
    pdf(file="targetInhibitionNetwork.pdf", width=figure.width, height=figure.height,pointsize=12)
    par(mar=c(0,0,0,0))
	# no axes, no labels
    plot(dummy, dummy, type = "n", axes = FALSE, ann = FALSE, xlim = c(0, figure.width), ylim = c(0, figure.height)) 
    start.line <- seg.in
    
    for (i in 1:components.num) {
        leg <- components.height[[i]]
        target.len <- c()
        for(j in 1:leg){
          target.len <- c(target.len, nchar(expressions.compact[[i]][j]))
        }
        len.max <- max(target.len) # the length of the longest target names in the component

        y0 <- figure.height/2 - (leg - 1)/2 # y for the lowest target
        y0 <- seq(y0, y0 + leg - 1) # y for the other targets

        x0 <- rep(start.line, leg)
        x1 <- rep(start.line+seg.in, leg)
        y1 <- y0
        x2 <- x1 + (seg - 2*seg.in)*len.max/3
        x3 <- x2 + seg.in
        
        segments(x0, y0, x1, y1)
        segments(x2, y0, x3, y1)
        lines(x0[c(1, 1)], y0[c(1, leg)], type = "l")
        lines(x3[c(1, 1)], y0[c(1, leg)], type = "l")
        lines(c(x3[1], x3[1] + seg.out), c(figure.height/2, figure.height/2), type = "l")
        
        start.line <- x3[1] + seg.out
        
        # write names
        text(x1, y0, labels = expressions.compact[[i]], pos = 4)
    }
    dev.off()
    
    # two terminal version, output to cytoscape nnf file
    cat(c(),file = "targetInhibitionNetwork.nnf")
    for (i in 1:components.num){
      write(paste(c('TargetInhibitionNetwork',paste('M',i,sep="")),collapse="\t"),file="targetInhibitionNetwork.nnf",sep="\n",append=T)
    }
    
    Terminals <- paste("T",1:2, sep = "") # the number of terminals for the SIF format
    
    for (i in 1:components.num){
      for (j in expressions.compact[[i]]){
      lines1 <-paste(c(paste('M',i,sep = ""),Terminals[1],"pp",j), collapse = "\t")
      lines2 <-paste(c(paste('M',i,sep = ""),j,"pp",Terminals[2]), collapse = "\t")
      write(lines1,file="targetInhibitionNetwork.nnf",sep = "\n",append = T)
      write(lines2,file="targetInhibitionNetwork.nnf",sep = "\n",append = T)
      }
    } 
    cat()
} 
