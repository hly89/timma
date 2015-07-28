# timma

requirement: Rcpp, QCA, iGraph, RCytoscape package 

usage:

1. loading R functions into environment <br />
   file.sources = list.files(pattern="*.R") <br />
   sapply(file.sources,source,.GlobalEnv)
2. loading Rcpp functions into environment <br />
   library(Rcpp) <br />
   sourceCpp("max_min.cpp")


how to visualize in Cytoscape

1. pass through iGraph object to VisualizeCytoscape <br />
   VisualizeCytoscape(graph)
2. after step 1, you get the following two files for Cytoscape <br />
   edge.attributes.cytoscape.txt and node.attributes.cytoscape.txt
3. create network in Cytoscape
   File ->

	
   
   
  