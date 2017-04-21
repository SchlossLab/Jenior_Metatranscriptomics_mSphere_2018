
# Set up environment
rm(list=ls())
gc()

# Load in functions
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define files

# Normalized Metatranscriptomes
noabx_normalized_reads <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/conv_normalized.tsv'
cef_normalized_reads <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/cef_normalized.tsv'
clinda_normalized_reads <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/clinda_normalized.tsv'
strep_normalized_reads <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/strep_normalized.tsv'

# KEGG pathway annotations for genes
noabx_pathways <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/conv_pathways.tsv'
cef_pathways <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/cef_pathways.tsv'
clinda_pathways <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/clinda_pathways.tsv'
strep_pathways <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/strep_pathways.tsv'

# KEGG taxonomy IDs
kegg_tax <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/kegg_taxonomy.tsv'

# Taxonomy colors
tax_colors <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/taxonomy_color.tsv'

# Output plot
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/figures/figure_2.pdf'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data

# Normalized Metatranscriptomes
noabx_normalized_reads <- read.delim(noabx_normalized_reads, sep='\t', header=TRUE, row.names=5)
noabx_normalized_reads[is.na(noabx_normalized_reads)] <- 'none'
cef_normalized_reads <- read.delim(cef_normalized_reads, sep='\t', header=TRUE, row.names=6)
cef_normalized_reads[is.na(cef_normalized_reads)] <- 'none'
clinda_normalized_reads <- read.delim(clinda_normalized_reads, sep='\t', header=TRUE, row.names=6)
clinda_normalized_reads[is.na(clinda_normalized_reads)] <- 'none'
strep_normalized_reads <- read.delim(strep_normalized_reads, sep='\t', header=TRUE, row.names=6)
strep_normalized_reads[is.na(strep_normalized_reads)] <- 'none'

# Pathway mapping data
noabx_pathways <- read.delim(noabx_pathways, sep='\t', header=TRUE, row.names=2)
cef_pathways <- read.delim(cef_pathways, sep='\t', header=TRUE, row.names=3)
clinda_pathways <- read.delim(clinda_pathways, sep='\t', header=TRUE, row.names=3)
strep_pathways <- read.delim(strep_pathways, sep='\t', header=TRUE, row.names=3)

# KEGG organism file
kegg_tax <- read.delim(kegg_tax, sep='\t', header=TRUE)
kegg_tax[] <- lapply(kegg_tax, as.character)

# Taxonomy colors
tax_colors <- read.delim(tax_colors, sep='\t', header=TRUE)
tax_colors[] <- lapply(tax_colors, as.character)

#-------------------------------------------------------------------------------------------------------------------------#

# Screen for those genes that have a gene annotation
noabx_annotated <- subset(noabx_normalized_reads, gene != 'none')
cef_annotated <- subset(cef_normalized_reads, gene != 'none')
clinda_annotated <- subset(clinda_normalized_reads, gene != 'none')
strep_annotated <- subset(strep_normalized_reads, gene != 'none')
rm(noabx_normalized_reads, cef_normalized_reads, clinda_normalized_reads, strep_normalized_reads)

# Find patterns of highest expression in each community


noabx_annotated <- noabx_annotated[order(-cyl),] 




#--------------------------------#

# Format pathway mapping ionformation

# Renames columns
colnames(noabx_pathways) <- c('noabx_mock')
colnames(cef_pathways) <- c('cef_630','cef_mock')
colnames(clinda_pathways) <- c('clinda_630','clinda_mock')
colnames(strep_pathways) <- c('strep_630','strep_mock')


#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()





