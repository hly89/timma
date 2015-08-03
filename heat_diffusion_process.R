#' Heat diffusion process function
#'
#' A function to implement the heat diffusion process algorithm
#' @param net an adjacency matrix of each node in the undirected network.
#' @param heat a vector of initial heat for each node.
#' @param t an integer to specify the number of iterations.
#' @param alpha the thermal conductivity.
#' @return a plot of the heat of each node at each iteration and a list of the following components:
#' \item steady.state a vector of the heat of each node at steady state.
#' \item heat.change a matrix of the heat of each node at each iteration.
#' @author Liye He \email{liye.he@@helsinki.fi}
#' @references Ma, Hao and Yang, Haixuan and Lyu, Michael R. and King, Irwin. 
#' Mining Social Networks Using Heat Diffusion Processes for Marketing Candidates Selection.
#' Proceedings of the 17th ACM Conference on Information and Knowledge Management, 2008, p233-242.
#' @examples
#' library(igraph) # v1.0.1 or later
#' library(ggplot2)
#' library(Matrix)
#' # g <- sample_growing(10, 1, directed = FALSE, citation = TRUE)
#' n <- 100
#' g <- sample_pa(n,directed=F) # scale-free network
#' net <- as.matrix(get.adjacency(g))
#' # heat <- c(1,0,0,0,1,0,1,0,1,0)
#' heat <- rep(0,n)
#' heat[sample(1:n, 1)] = 1
#' results <- heatDiffusionProcess(net, heat, 50)

heatDiffusionProcess <- function(net, heat, t = 50, alpha = 1) {
	# net should be adjaceny matrix
	
	# get the transition matrix
	net2 <- net - diag(rowSums(net))
	trans.mat <- expm(net2)
	# get the steady-state vector
	n <- ncol(trans.mat)
	a <- t(trans.mat-diag(n))
	a <- rbind(a, rep(1,n))
	b <- c(rep(0,n), 1)
	mu<-qr.solve(a,b)

	idx <- seq(0.1, t*0.1, 0.1)
	# nrow(net) : num of nodes
	nodes.num <- nrow(net2)
	heat.info <- matrix(0, nrow = nodes.num, ncol = t+1)
	heat.info[ , 1] <- heat
	for (i in seq_len(t)) {
		heat.info[ , i+1] <- as.vector(heat %*% expm(alpha*idx[i]*net2))
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
	  return(list(steady.state = sum(heat)*mu, heat.change = heat.info, fig = fig))
}