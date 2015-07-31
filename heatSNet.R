heatSNet <- function(net, heat, t) {
  # nodes num
  nodes.num <- nrow(net)
  heat.info <- matrix(0, nrow = nodes.num, ncol = t+1)
  heat.info[,1] <- heat
  net.degree <- rowSums(net)
  # the transition matrix
  trans.mat <- t(net/net.degree)
  # get the steady-state vector
  n <- ncol(trans.mat)
  a <- t(trans.mat-diag(n))
  a <- rbind(a, rep(1,n))
  b <- c(rep(0,n), 1)
  mu<-qr.solve(a,b)
  
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
  ggplot(data = heat.change.df, aes(x = Time, y = Heat, group = Node, color = Node)) + geom_line() + geom_point()
  ggsave(filename = "heatS.png")
  return(list(steady.state = mu, heat.change = heat.info))
}