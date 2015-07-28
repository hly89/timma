#' Visualize network in Cytoscape 2.x/3.x
#'
#' A function to visualize the network in Cytoscape 2.x/3.x
#' 
#' @param graph an iGraph object. 
#' @param mode whether to visualize the network automatically or manually. If mode is "automatic", the network is visualized automatically in Cytoscape, which requires CytoscapeRPC plugin installed beforehand. If mode is "manual", the nodes and edges attributes files are generated for loading into Cytoscape.
#'
#' @author Liye He \email{liye.he@@helsinki.fi}
#' @examples
#' \dontrun{
#' data(graph)
#' VisualizeNetCYT(graph)
#'}
VisualizeNetCYT <- function(graph, mode = "automatic") {
	#------------------
	# Visualization in Cytoscape 2.8/3.0
	#------------------
	if (!(mode %in% c("automatic", "manual"))) {
		stop ("The mode parameter can only be either automatic or manual!")
	}
    # the required attributes for graph: node.fillColor, node.shape
	
	# no node.fillColor attributes
	if (is.null(V(graph)$node.fillColor)) {
		graph <- set.vertex.attribute(graph, "node.fillColor", value = "255,0,255") # default color
	} 
	
	# node shape 
	if (is.null(V(graph)$node.shape)) {
		graph <- set.vertex.attribute(graph, "node.shape", value = "ellipse") # default shape
	}
	
	
	
	# edge color
	if (is.null(E(graph)$edge.color)) {
		graph <- set.edge.attribute(graph, "edge.color", value = "0,255,0")
	}
	
	
	
	# generate edge arrow
	edges.type <- E(graph)$kind
	edges.arrow <- vector(mode = "character", length = length(edges.type))
	edges.arrow[grep("Interaction", edges.type)] <- "NONE" 
	edges.arrow[grep("Activation", edges.type)] <- "ARROW" 
	edges.arrow[grep("Inhibition", edges.type)] <- "T"
	
	graph <- set.edge.attribute(graph, "edge.targetArrowShape", value = edges.arrow)
	
	# generate edge line width
	graph <- set.edge.attribute(graph, "edge.lineWidth", value = E(graph)$weight)
	
	
	if (mode == "automatic") {
		# convert to graphNEL 
		graph.cyt <- igraph.to.graphNEL(graph)
	
		# create attributes for graphNEL
		graph.cyt <- initNodeAttribute(graph.cyt, 'Symbol', 'char', "symbol")	
		graph.cyt <- initNodeAttribute(graph.cyt, 'EntrezID', 'numeric', 0)
		graph.cyt <- initNodeAttribute(graph.cyt, 'Alterations', 'char', 'exp')
		graph.cyt <- initNodeAttribute(graph.cyt, 'node.fillColor', 'char', '84,84,84')
		graph.cyt <- initNodeAttribute(graph.cyt, 'node.shape', 'char', 'ellipse')
		graph.cyt <- initNodeAttribute(graph.cyt, 'Value', 'numeric', 0)
		graph.cyt <- initEdgeAttribute(graph.cyt, 'weight', 'numeric', 0)
	    graph.cyt <- initEdgeAttribute(graph.cyt, 'kind', 'char', 'interaction')
	    graph.cyt <- initEdgeAttribute(graph.cyt, 'edge.color', 'char', '255,0,0')
	    graph.cyt <- initEdgeAttribute(graph.cyt, 'edge.targetArrowShape', 'char', 'NONE')
	    graph.cyt <- initEdgeAttribute(graph.cyt, 'edge.lineWidth', 'char', '1')
		# start the clock
		ptm <- proc.time()
		os.type <- Sys.info()["sysname"]
		if (os.type == "Darwin") { # mac os
			system("open -a /Applications/Cytoscape_v2.8.2/Cytoscape.app")
			# check if the port is open
			port.info <- suppressWarnings(system("lsof -i :9000", intern = TRUE))
			while (length(port.info) == 0) {
				port.info <- suppressWarnings(system("lsof -i :9000", intern = TRUE))
				running.time <- proc.time() - ptm
				if (running.time[3] > 180) {
					stop ("Taking too much time to open Cytoscape, please open Cytoscape manually and run again!")
				}
			}
		    cyt.window <- new.CytoscapeWindow("Network", graph = graph.cyt, overwriteWindow = TRUE)
			displayGraph(cyt.window)
		    cy <- CytoscapeConnection ()
		    layout.name <- getLayoutNames(cy)[15]
		    layoutNetwork(cyt.window, layout.name)
		    setDefaultNodeSize(cyt.window, 70)
			redraw(cyt.window)
		} else {
		    cyt.window <- new.CytoscapeWindow("Network", graph = graph.cyt, overwriteWindow = TRUE)
			displayGraph(cyt.window)
		    cy <- CytoscapeConnection ()
		    layout.name <- getLayoutNames(cy)[15]
		    layoutNetwork(cyt.window, layout.name)
		    setDefaultNodeSize(cyt.window, 70)
			redraw(cyt.window)
		}
	    
		
	    
	} else {
		# the edge list in Symbol
		edge.list.symbol <- get.edgelist(graph)
		# the edge list in Symbol
		nodes.entrez <- V(graph)$EntrezID
		nodes.symbol <- V(graph)$Symbol
		edges.attr.table <- cbind(edge.list.symbol, edges.color, E(graph)$weight, edges.arrow)
		colnames(edges.attr.table) <- c("source", "target", "edge.color", "edge.Width", "edge.targetArrowShape")
		write.table(edges.attr.table, file = "edge.attributes.cytoscape.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
		# generate the node attribute table
		node.attr.table <- cbind(nodes.symbol, nodes.type, nodes.alteration, V(graph)$node.fillColor, V(graph)$Value)
		colnames(node.attr.table) <- c("name", "node.shape", "type", "node.fillColor", "value")
		write.table(node.attr.table, file = "node.attributes.cytoscape.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
	}
	

}
