heatDiffusionProcess <- function(net, heat, t) {
	# net should be adjaceny matrix
	
	# get the transition matrix
	net <- net - diag(rowSums(net))
	trans.mat <- expm(net)
	# get the steady-state vector
	n <- ncol(trans.mat)
	a <- t(trans.mat-diag(n))
	a <- rbind(a, rep(1,n))
	b <- c(rep(0,n), 1)
	mu<-qr.solve(a,b)

	idx <- seq(0.1, t*0.1, 0.1)
	# nrow(net) : num of nodes
	heat.info <- matrix(0, nrow = nrow(net), ncol = t+1)
	heat.info[ , 1] <- heat
	for (i in seq_len(t)) {
		heat.info[ , i+1] <- as.vector(heat.info[ , i] %*% expm(idx[i]*net))
	}
	# plot the heat changes for each nodes
	heat.change.mat <- matrix(0, nrow = 10*(t+1), ncol = 3)
    node.list <- c()
    for (i in 1:nrow(net)) {
      node.list <- c(node.list, rep(i, t+1))
    }
    heat.change.mat[ , 1] <- node.list
    heat.change.mat[ , 2] <- as.vector(t(heat.info))
    heat.change.mat[ , 3] <- rep(c(1:51), 10)
    colnames(heat.change.mat) <- c("Node", "Heat", "Time")  
    heat.change.df <- data.frame(heat.change.mat)
    heat.change.df$Node<-as.factor(heat.change.df$Node)
	ggplot(data = heat.change.df, aes(x = Time, y = Heat, group = Node, color = Node)) + geom_line() + geom_point()
	ggsave(filename = "heatDP.png")
	return(list(steady.state = mu, heat.change = heat.info))
}