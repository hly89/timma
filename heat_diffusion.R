#' Heat diffusion functions
#'
#' A function to implement heatS and probS algorithms for undirected/directed networks.
#' @param graph an igraph object
#' @param heat a vector of initial heat for each node.
#' @param t an integer to specify the number of iterations
#' @param method a string to specify the heat diffusion algorithm. In heatS algorithm, each node receives the mean amouont of the heat possessed by its neighboring nodes. In probS algorithm, heat possessed by each node is evenly distributed to its neighboring nodes. In probS, the total heat remains constant while heat is not so.
#' @return a plot of the heat of each node at each iteration and a list of the following components:
#' \item steady.state a vector of the heat of each node at steady state if the steady state can be determined.
#' \item heat.change a matrix of the heat of each node at each iteration.
#' @author Liye He \email{liye.he@@helsinki.fi}
#' @references T. Zhou, Z. Kuscsik, J.G. Liu, M. Medo, J.R. Wakeling, and Y.C. Zhang.
#' Solving the apparent diversity-accuracy dilemma of recommender systems. Proceedings of the National Academy of Sciences 
#' 107(10):4511--4515 (2010)
#' @examples
#' library(igraph) # v1.0.1 or later
#' library(ggplot2)
#' # g <- sample_growing(10, 1, directed = FALSE, citation = TRUE)
#' n <- 100
#' g <- sample_pa(n,directed=T) # scale-free network
#' # heat <- c(1,0,0,0,1,0,1,0,1,0)
#' heat <- rep(0,n)
#' heat[sample(1:n, 1)] = 1
#' results <- heatDiffusion(g, heat, 100, method="heatS")


heatDiffusion <- function(graph, heat, t = 50, method = "probS") {
  # graph is an igraph object
  
  if (!(method %in% c("probS", "heatS"))) {
     stop("Parameter method can only be either probS or heatS")
  }
  nodes.num <- vcount(graph)
  heat.info <- matrix(0, nrow = nodes.num, ncol = (t + 1))
  heat.info[,1] <- heat
  
  net <- as.matrix(get.adjacency(graph))
  if (method == "probS") {
    if (is_directed(graph)) {
      trans.mat <- net/degree(graph, mode = "out")
      # get nan idx
      nan.idx <- which(is.nan(trans.mat) == TRUE)
      trans.mat[nan.idx] <- 0
    } else {
      trans.mat <- net/degree(graph)
    }
  	
  } else {
  	if (is_directed(graph)) {
  	  trans.mat <- t(t(net)/degree(graph, mode = "in"))
  	  # get nan idx
  	  nan.idx <- which(is.nan(trans.mat) == TRUE)
  	  trans.mat[nan.idx] <- 0
  	} else {
  	  trans.mat <- t(net/degree(graph))
  	}
  }
  
  # get the steady-state vector
  n <- ncol(trans.mat)
  a <- t(trans.mat-diag(n))
  a.det <- det(a)
  if (a.det == 0) {
    cat("No steady state vector can be determined. \n")
  } else {
    a <- rbind(a, rep(1,n))
    b <- c(rep(0,n), 1)
    mu<-qr.solve(a,b)
  }
  

  for (i in seq_len(t)) {
    heat.info[ , i+1] <- heat.info[ , i] %*% trans.mat
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
  fig = ggplot(data = heat.change.df, aes(x = Time, y = Heat, group = Node, color = Node)) + geom_line() + geom_point()
  if (method == "probS") {
  	ggsave(filename = "probS.png")
  } else {
  	ggsave(filename = "heatS.png")
  }
  
  if (a.det == 0) {
    return(list(heat.change = heat.info, fig = fig))
  } else {
    return(list(steady.state = sum(heat)*mu, heat.change = heat.info, fig = fig))
  }
  
}





