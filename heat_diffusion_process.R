#' Heat diffusion process function
#'
#' A function to implement the heat diffusion process algorithm
#' @param graph an igraph object.
#' @param heat a vector of initial heat for each node.
#' @param t an integer to specify the number of iterations.
#' @param alpha the thermal conductivity.
#' @return a plot of the heat of each node at each iteration and a list of the following components:
#' \item steady.state a vector of the heat of each node at steady state if the steady state can be determined.
#' \item heat.change a matrix of the heat of each node at each iteration.
#' \item fig the plot object.
#' @author Liye He \email{liye.he@@helsinki.fi}
#' @references Ma, Hao and Yang, Haixuan and Lyu, Michael R. and King, Irwin. 
#' Mining Social Networks Using Heat Diffusion Processes for Marketing Candidates Selection.
#' Proceedings of the 17th ACM Conference on Information and Knowledge Management, 2008, p233-242.
#' @examples
#' library(igraph) # v1.0.1 or later
#' library(ggplot2)
#' library(Matrix) # v1.2-2
#' # g <- sample_growing(10, 1, directed = TRUE, citation = TRUE)
#' n <- 100
#' g <- sample_pa(n, directed = FALSE) # scale-free network
#' net <- as.matrix(get.adjacency(g))
#' # heat <- c(1,0,0,0,1,0,1,0,1,0)
#' heat <- rep(0,n)
#' heat[sample(1:n, 1)] <- 1
#' results <- heatDiffusionProcess(g, heat, 50)

heatDiffusionProcess <- function(graph, heat, t = 50, alpha = 1) {
	# graph is an igraph object
	
  # get the adjacency matrix from an igraph object
  net <- as.matrix(get.adjacency(graph))
  # nrow(net) : num of nodes
  nodes.num <- nrow(net)
  if (is_directed(graph)) {
    # get the transition matrix for directed network
    node.outdegree <- degree(graph, mode = "out")
    tau <- node.outdegree
    tau[which(tau != 0)] <- 1
    node.outdegree[which(node.outdegree ==  0)] <- 1
    mat.tmp <- matrix(rep(node.outdegree, nodes.num), nodes.num, nodes.num)
    net2 <- t(net)/t(mat.tmp)
    diag(net2) <- -tau
  } else {
    # get the transition matrix for undirected network
    net2 <- net - diag(rowSums(net))
  }
  trans.mat <- as.matrix(expm(alpha*net2))
  # get the steady-state vector
  n <- ncol(trans.mat)
  a <- t(trans.mat-diag(n))
  a <- as.matrix(a)
  a <- rbind(a, rep(1,n))
  # check if a is a square matrix
  a.tmp <- unique(a)
  mat.singular <- FALSE
  if (ncol(a.tmp) == nrow(a.tmp)) { # a square matrix
    a.det <- det(a.tmp)
    # check if singular matrix
    if (a.det == 0) {
      # singluar matrix
      mat.singular <- TRUE
      cat("No steady state vector can be determined! \n")
    } else {
      # nonsingluar matrix
      b <- c(rep(0,n), 1)
      mu<-qr.solve(a,b)
    }
  } else { # not a square matrix
    b <- c(rep(0,n), 1)
    mu<-qr.solve(a,b)
  }
  
  

	idx <- seq(0.1, t*0.1, 0.1)
	
	heat.info <- matrix(0, nrow = nodes.num, ncol = t+1)
	heat.info[ , 1] <- heat
	for (i in seq_len(t)) {
		heat.info[ , i+1] <- as.vector(heat %*% as.matrix(expm(alpha*idx[i]*net2)))
	}
	# plot the heat changes for each nodes
	heat.change.mat <- matrix(0, nrow = nodes.num*(t+1), ncol = 3)
    node.list <- c()
    for (i in 1:nrow(net)) {
      node.list <- c(node.list, rep(i, t+1))
    }
    heat.change.mat[ , 1] <- node.list
    heat.change.mat[ , 2] <- as.vector(t(heat.info))
    heat.change.mat[ , 3] <- rep(c(1:(t + 1)), nodes.num)
    colnames(heat.change.mat) <- c("Node", "Heat", "Time")  
    heat.change.df <- data.frame(heat.change.mat)
    heat.change.df$Node<-as.factor(heat.change.df$Node)
	  fig = ggplot(data = heat.change.df, aes(x = Time, y = Heat, group = Node, color = Node)) + geom_line() + geom_point()
	  ggsave(filename = "heatDP.png")
	  if (mat.singular) {
	    return(list(heat.change = heat.info, fig = fig))
	  } else {
	    return(list(steady.state = sum(heat)*mu, heat.change = heat.info, fig = fig))
	  }
	  
}