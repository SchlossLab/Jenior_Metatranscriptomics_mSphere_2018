
# Install packages
#install.packages('vegan')
library(vegan)

# Define variables
cef_input <- '/Users/schloss/Desktop/cdf_bipartite/cefoperazone_630.bipartite.files/cefoperazone_630.input_score.txt'
strep_input <- '/Users/schloss/Desktop/cdf_bipartite/clindamycin_630.bipartite.files/clindamycin_630.input_score.txt'
clinda_input <- '/Users/schloss/Desktop/cdf_bipartite/streptomycin_630.bipartite.files/streptomycin_630.input_score.txt'
gf_input <- '/Users/schloss/Desktop/cdf_bipartite/germfree.bipartite.files/germfree.input_score.txt'

#--------------------------------------------------------------------------------------------------------------#

# Define formatting function
format_compound <- function(file_name, sample_name){
  scores <- read.delim(file_name, header = TRUE, row.names = 1)
  colnames(scores) <- c('kegg_code', sample_name)
  scores <- scores[order(-scores[,2]), ] 
  scores$kegg_code <- NULL
  return(scores)
}

#--------------------------------------------------------------------------------------------------------------#

gene_expression <- merge(cef_expression_all, clinda_expression_all, by='row.names')
rownames(gene_expression) <- gene_expression$Row.names
gene_expression$Row.names <- NULL
gene_expression <- merge(gene_expression, strep_expression_all, by='row.names')
rownames(gene_expression) <- gene_expression$Row.names
gene_expression$Row.names <- NULL
gene_expression <- merge(gene_expression, gf_expression_all, by='row.names')
rownames(gene_expression) <- gene_expression$Row.names
gene_expression$Row.names <- NULL
colnames(gene_expression) <- c('Cefoperazone', 'Clindamycin', 'Streptomycin', 'Germfree')

#--------------------------------------------------------------------------------------------------------------#

# Prep all of the data
cef_input_scores <- format_compound(cef_input, 'Cefoperazone')
clinda_input_scores <- format_compound(clinda_input, 'Clindamycin')
strep_input_scores <- format_compound(strep_input, 'Streptomycin')
gf_input_scores <- format_compound(gf_input, 'Germfree')

# Merge compound tables
input_scores <- merge(cef_input_scores , clinda_input_scores, by = 'row.names')
rownames(input_scores) <- input_scores$Row.names
input_scores$Row.names <- NULL
input_scores <- merge(input_scores , strep_input_scores, by = 'row.names')
rownames(input_scores) <- input_scores$Row.names
input_scores$Row.names <- NULL
input_scores <- merge(input_scores , gf_input_scores, by = 'row.names')
rownames(input_scores) <- input_scores$Row.names
input_scores$Row.names <- NULL
colnames(input_scores) <- c('Cefoperazone', 'Clindamycin', 'Streptomycin', 'Germfree')

# Remove categories where all of the values are 0 and sort
filter_input_scores <- input_scores[rowSums(input_scores[, -1]) > 0, ]
filter_input_scores <- filter_input_scores[order(rowSums(filter_input_scores), decreasing = T), ]
filter_input_scores <- as.matrix(filter_input_scores)

# Generate distance matrix
dist_matrix <- vegdist(filter_input_scores, method='bray')
nmds <- metaMDS(dist_matrix)

# Plot 
stressplot(nmds)
plot(nmds)

#--------------------------------------------------------------------------------------------------------------#

# PCA ordination

# Input compounds
pca <- rda(filter_input_scores)
biplot(pca, scal=2, xlab='PC1 (88.9%)', ylab='PC2 (5.5%)')
text(c(-2,-2.5,-3,-2.2), c(1.25,0.35,-0.35,-0.65), c('Germfree','Cefoperazone','Clindamycin','Streptomycin'))

# Pathways
pca2 <- rda(filter_all_expression)
biplot(pca2, xlab='PC1 (88.6%)', ylab='PC2 (10.6%)')

# Genes
pdf(file='/Users/schloss/Desktop/genes.biplot.pdf', width=10, height=10)
pca3 <- rda(gene_expression)
biplot(pca3, scal=2, xlab='PC1 (82.1%)', ylab='PC2 (11.3%)')
text(c(3,0.5,1,1.3), c(0.1,-2,-1.7,-1.3), c('Germfree','Cefoperazone','Clindamycin','Streptomycin'))
dev.off()

#--------------------------------------------------------------------------------------------------------------#


