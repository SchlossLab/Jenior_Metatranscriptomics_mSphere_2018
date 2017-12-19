
# Load dependencies
if ('igraph' %in% installed.packages()[,"Package"] == FALSE){install.packages('igraph', quiet=TRUE)}
library('igraph', verbose=FALSE, character.only=TRUE)

# Function to analyse connectedness of metabolic network
analyze_network <- function(network_file) {
  
  # Read in metabolic network data
  network <- read.delim(network_file, header=FALSE, sep='\t')
  rm(network_file)
  
  # Format directed graph
  raw_graph <- graph.data.frame(network, directed=TRUE)
  simple_graph <- simplify(raw_graph)
  rm(network, raw_graph)
  
  # Calculate summary statistice
  edges <- length(as.vector(E(simple_graph)))
  nodes <- length(as.vector(V(simple_graph)))
  kos <- length(as.vector(grep('K', V(simple_graph)$name, value=TRUE)))
  compounds <- length(as.vector(grep('C', V(simple_graph)$name, value=TRUE)))
  comps <- components(simple_graph)$no
  med_path <- median(distances(simple_graph))
  between <- median(betweenness(simple_graph))
  close <- median(closeness(simple_graph, vids=V(simple_graph), mode='total'))
  
  # Generate output
  output <- c(edges,nodes,kos,compounds,comps,med_path,between,close)
  names(output) <- c('total_edges','total_nodes','KOs','metabolites','component_graphs','median_shortest_paths','median_betweenness','median_closeness')
  return(output)
}

#-------------------------------------------------------------------------------------------------------------------------------------#

# cef
akkermansia_cef <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/cefoperazone/infected/akkermansia_630.bipartite.files/graph.tsv')
bacteroides_cef <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/cefoperazone/infected/bacteroides_630.bipartite.files/graph.tsv')
blautia_cef <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/cefoperazone/infected/blautia_630.bipartite.files/graph.tsv')
clostridium_cef <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/cefoperazone/infected/clostridium_630.bipartite.files/graph.tsv')
escherichia_cef <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/cefoperazone/infected/escherichia_630.bipartite.files/graph.tsv')
lactobacillus_cef <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/cefoperazone/infected/lactobacillus_630.bipartite.files/graph.tsv')

# strep
bacteroides_strep <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/streptomycin/infected/bacteroides_630.bipartite.files/graph.tsv')
bifidobacterium_strep <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/streptomycin/infected/bifidobacterium_630.bipartite.files/graph.tsv')
clostridium_strep <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/streptomycin/infected/clostridium_630.bipartite.files/graph.tsv')
escherichia_strep <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/streptomycin/infected/escherichia_630.bipartite.files/graph.tsv')
lactobacillus_strep <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/streptomycin/infected/lactobacillus_630.bipartite.files/graph.tsv')
olsenella_strep <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/streptomycin/infected/olsenella_630.bipartite.files/graph.tsv')
parabacteroides_strep <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/streptomycin/infected/parabacteroides_630.bipartite.files/graph.tsv')

# clinda
bacteroides_clinda <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/clindamycin/infected/bacteroides_630.bipartite.files/graph.tsv')
escherichia_clinda <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/clindamycin/infected/escherichia_630.bipartite.files/graph.tsv')
lachnospiraceae_clinda <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/clindamycin/infected/lachnospiraceae_630.bipartite.files/graph.tsv')
lactobacillus_clinda <- analyze_network('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolic_models/clindamycin/infected/lactobacillus_630.bipartite.files/graph.tsv')



# resistant

#-------------------------------------------------------------------------------------------------------------------------------------#



