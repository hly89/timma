#' Heat diffusion functions
#'
#' A function to implement heatS and probS algorithms for undirected/directed networks.
#' @param graph an igraph object with gene.exp as node attributes.
#' @param heat a vector of initial heat for each node.
#' @param t an integer to specify the number of iterations.
#' @return a plot of the heat of each node at each iteration and a list of the following components:
#' \item steady.state a vector of the heat of each node at steady state if it can be determined.
#' \item heat.change a matrix of the heat of each node at each iteration.
#' \item fig the plot object.
#' @author Liye He \email{liye.he@@helsinki.fi}
#' @references T. Zhou, Z. Kuscsik, J.G. Liu, M. Medo, J.R. Wakeling, and Y.C. Zhang.
#' Solving the apparent diversity-accuracy dilemma of recommender systems. Proceedings of the National Academy of Sciences 
#' 107(10):4511--4515 (2010)
#' @examples
#' library(igraph) # v1.0.1 or later
#' library(ggplot2)
#' n <- 100
#' g <- sample_pa(n,directed = FALSE) # scale-free network
#' gene.exp <- runif(n, 0.1, 100)
#' g <- set.vertex.attribute(g, name = "gene.exp", value = gene.exp)
#' heat <- rep(0, n)
#' heat[sample(1:n, 1)] <- 1
#' results <- heatDiffusionWeighted(g, heat, 100)


heatDiffusionWeighted <- function(graph, heat, t = 50) {
  # graph is an igraph object
  nodes.num <- vcount(graph)
  heat.info <- matrix(0, nrow = nodes.num, ncol = (t + 1))
  heat.info[,1] <- heat
  
  gene.exp <- V(graph)$gene.exp
  gene.exp.table <- t(matrix(rep(gene.exp, nodes.num), nodes.num, nodes.num))
  net <- as.matrix(get.adjacency(graph))
  if (is_directed(graph)) {
      net.tmp <- t(t(gene.exp.table)*net)/gene.exp
      #trans.mat <- net.tmp / rowSums(net.tmp)
   
    } else {
      net.tmp <- gene.exp.table * net
      net.tmp <- net.tmp / gene.exp
  }
  trans.mat <- net.tmp / rowSums(net.tmp)
  nan.idx <- which(is.nan(trans.mat) == TRUE)
  trans.mat[nan.idx] <- 0
  # get the steady-state vector
  n <- ncol(trans.mat)
  a <- t(trans.mat-diag(n))
  a <- rbind(a, rep(1,n))
  # check if a is a square matrix
  a.tmp <- unique(a)
  mat.singular <- FALSE
  if (ncol(a.tmp) == nrow(a.tmp)) { # a square matrix
    a.det <- det(a.tmp)
    # check if singular matrix
    if (a.det == 0) {
      # singluar matrix
      cat("No steady state vector can be determined! \n")
      mat.singular <- TRUE
    } else {
      # nonsingluar matrix
      b <- c(rep(0,n), 1)
      mu<-qr.solve(a,b)
    }
  } else { # not a square matrix
    b <- c(rep(0,n), 1)
    mu<-qr.solve(a,b)
  }
  
  
  for (i in seq_len(t)) {
    heat.info[ , i+1] <- as.vector(heat.info[ , i] %*% trans.mat)
  }
  
  
  
  heat.change.mat <- matrix(0, nrow = nodes.num*(t + 1), ncol = 3)
  node.list <- c()
  for (i in 1:nodes.num) {
    node.list <- c(node.list, rep(i, t+1))
  }
  heat.change.mat[ , 1] <- node.list
  heat.change.mat[ , 2] <- as.vector(t(heat.info))
  heat.change.mat[ , 3] <- rep(c(1:(t + 1)), nodes.num)
  colnames(heat.change.mat) <- c("Node", "Heat", "Time")
  
  
  heat.change.df <- data.frame(heat.change.mat)
  heat.change.df$Node<-as.factor(heat.change.df$Node)
  fig <- ggplot(data = heat.change.df, aes(x = Time, y = Heat, group = Node, color = Node)) + geom_line() + geom_point()
  ggsave("heatW.png")
  
  if (mat.singular) {
    return(list(heat.change = heat.info, fig = fig))
  } else {
    return(list(steady.state = sum(heat)*mu, heat.change = heat.info, fig = fig))
  }
  
}





