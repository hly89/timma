#------------------
# Visualization in Cytoscape 2.8/3.0
#------------------
VisualizeCytoscape <- function (graph) {
	# the edge list in EntrezID
	edge.list.entrez <- get.edgelist(graph)
	# the edge list in Symbol
	nodes.entrez <- V(graph)$EntrezID
	nodes.symbol <- V(graph)$Symbol
	edge.list.symbol <- apply(edge.list.entrez, c(1,2), function(x) nodes.symbol[which(nodes.entrez == x)])
	
	# generate edge color
	edges.type <- E(graph)$kind
	edges.color <- vector(mode = "character", length = length(edges.type))
	edges.color[grep("Interaction", edges.type)] <- "0,255,0" # green edges for interaction
	edges.color[grep("Activation", edges.type)] <- "255,0,0" # red edges for activation
	edges.color[grep("Inhibition", edges.type)] <- "0,0,255" # blue edges for inhibition
	
	# generate edge width
	edges.width <- E(graph)$weight
	edges.width <- paste("w:", edges.width, sep = "")
	# generate edge arrow
	edges.arrow <- vector(mode = "character", length = length(edges.type))
	edges.arrow[grep("Interaction", edges.type)] <- "NONE" 
	edges.arrow[grep("Activation", edges.type)] <- "ARROW" 
	edges.arrow[grep("Inhibition", edges.type)] <- "T" 
	# generate the edge attributes file
	edges.attr.table <- cbind(edge.list.symbol, edges.color, edges.width, edges.arrow)
	colnames(edges.attr.table) <- c("source", "target", "edge.color", "edge.label", "edge.targetArrowShape")
	
	write.table(edges.attr.table, file = "edge.attributes.cytoscape.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
	
	# generate node attributes file: a table with columns as attributes
	# node color
	nodes.num <- vcount(graph)
	if (is.null(V(graph)$color)) {
		nodes.color <- rep("84,84,84", nodes.num)
	} else {
		nodes.color <- V(graph)$color
	}
	
	# node shape 
	shape.style <- c('ellipse','rectangle','triangle','diamond','hexagon','octagon','parallelogram','roundrect','vee')
	variable.type.list <- c("timma_parental","timma_tam","timma_parental timma_tam","cnv","exp","mutation","others")
	nodes.type <- c()
	nodes.alteration <- V(graph)$Alterations
	for (i in seq_len(nodes.num)) {
		tmp <- which(variable.type.list == nodes.alteration[i])
		if (length(tmp) == 0) {
			nodes.type[i] <- shape.style[7]
		} else {
			nodes.type[i] <- shape.style[tmp]
		}
	}
	# generate the node attribute table
	node.attr.table <- cbind(nodes.symbol, nodes.type, nodes.alteration, nodes.color, V(graph)$Value)
	colnames(node.attr.table) <- c("name", "node.shape", "type", "node.fillColor", "value")
	
	write.table(node.attr.table, file = "node.attributes.cytoscape.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
}



