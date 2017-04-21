
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

# Output plot
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/figures/figure_2.pdf'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data
# Normalized Metatranscriptomes
noabx_normalized_reads <- read.delim(noabx_normalized_reads, sep='\t', header=TRUE, row.names=5, stringsAsFactors=FALSE)
cef_normalized_reads <- read.delim(cef_normalized_reads, sep='\t', header=TRUE, row.names=6, stringsAsFactors=FALSE)
clinda_normalized_reads <- read.delim(clinda_normalized_reads, sep='\t', header=TRUE, row.names=6, stringsAsFactors=FALSE)
strep_normalized_reads <- read.delim(strep_normalized_reads, sep='\t', header=TRUE, row.names=6, stringsAsFactors=FALSE)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data
# Remove excess columns
noabx_normalized_reads$ko <- NULL
noabx_normalized_reads$pathway <- NULL
cef_normalized_reads$ko <- NULL
cef_normalized_reads$pathway <- NULL
clinda_normalized_reads$ko <- NULL
clinda_normalized_reads$pathway <- NULL
strep_normalized_reads$ko <- NULL
strep_normalized_reads$pathway <- NULL

# Screen for those genes that have a gene annotation
noabx_annotated <- subset(noabx_normalized_reads, gene != '')
noabx_annotated <- noabx_annotated[!rownames(noabx_annotated) %in% rownames(noabx_annotated[grep('unknown_\\d', noabx_annotated$gene),]), ]
cef_annotated <- subset(cef_normalized_reads, gene != '')
cef_annotated <- cef_annotated[!rownames(cef_annotated) %in% rownames(cef_annotated[grep('unknown_\\d', cef_annotated$gene),]), ]
clinda_annotated <- subset(clinda_normalized_reads, gene != '')
clinda_annotated <- clinda_annotated[!rownames(clinda_annotated) %in% rownames(clinda_annotated[grep('unknown_\\d', clinda_annotated$gene),]), ]
strep_annotated <- subset(strep_normalized_reads, gene != '')
strep_annotated <- strep_annotated[!rownames(strep_annotated) %in% rownames(strep_annotated[grep('unknown_\\d', strep_annotated$gene),]), ]
rm(noabx_normalized_reads, cef_normalized_reads, clinda_normalized_reads, strep_normalized_reads)

# Reverse log2 transformation temporarily
noabx_annotated$conv_metaT_reads <- 2 ^ noabx_annotated$conv_metaT_reads
cef_annotated$cef_630_metaT_reads <- 2 ^ cef_annotated$cef_630_metaT_reads
cef_annotated$cef_mock_metaT_reads <- 2 ^ cef_annotated$cef_mock_metaT_reads
clinda_annotated$clinda_630_metaT_reads <- 2 ^ clinda_annotated$clinda_630_metaT_reads
clinda_annotated$clinda_mock_metaT_reads <- 2 ^ clinda_annotated$clinda_mock_metaT_reads
strep_annotated$strep_630_metaT_reads <- 2 ^ strep_annotated$strep_630_metaT_reads
strep_annotated$strep_mock_metaT_reads <- 2 ^ strep_annotated$strep_mock_metaT_reads

# Subset to treatment groups
noabx_annotated$conv_metaT_reads <- as.numeric(noabx_annotated$conv_metaT_reads)
cef_630_annotated <- cef_annotated
cef_630_annotated$cef_630_metaT_reads <- as.numeric(cef_630_annotated$cef_630_metaT_reads)
cef_630_annotated$cef_mock_metaT_reads <- NULL
cef_mock_annotated <- cef_annotated
cef_mock_annotated$cef_mock_metaT_reads <- as.numeric(cef_mock_annotated$cef_mock_metaT_reads)
cef_mock_annotated$cef_630_metaT_reads <- NULL
clinda_630_annotated <- clinda_annotated
clinda_630_annotated$clinda_630_metaT_reads <- as.numeric(clinda_630_annotated$clinda_630_metaT_reads)
clinda_630_annotated$clinda_mock_metaT_reads <- NULL
clinda_mock_annotated <- clinda_annotated
clinda_mock_annotated$clinda_mock_metaT_reads <- as.numeric(clinda_mock_annotated$clinda_mock_metaT_reads)
clinda_mock_annotated$clinda_630_metaT_reads <- NULL
strep_630_annotated <- strep_annotated
strep_630_annotated$strep_630_metaT_reads <- as.numeric(strep_630_annotated$strep_630_metaT_reads)
strep_630_annotated$strep_mock_metaT_reads <- NULL
strep_mock_annotated <- strep_annotated
strep_mock_annotated$strep_mock_metaT_reads <- as.numeric(strep_mock_annotated$strep_mock_metaT_reads)
strep_mock_annotated$strep_630_metaT_reads <- NULL
rm(cef_annotated, clinda_annotated, strep_annotated)

# Aggregate identical genes, regardless of organism - retransform
noabx_annotated <- aggregate(noabx_annotated$conv_metaT_reads, by=list(noabx_annotated$gene), FUN=sum)
colnames(noabx_annotated) <- c('gene', 'read_abundance')
noabx_annotated$read_abundance <- log2(noabx_annotated$read_abundance)
cef_630_annotated <- aggregate(cef_630_annotated$cef_630_metaT_reads, by=list(cef_630_annotated$gene), FUN=sum)
colnames(cef_630_annotated) <- c('gene', 'read_abundance')
cef_630_annotated$read_abundance <- log2(cef_630_annotated$read_abundance)
cef_mock_annotated <- aggregate(cef_mock_annotated$cef_mock_metaT_reads, by=list(cef_mock_annotated$gene), FUN=sum)
colnames(cef_mock_annotated) <- c('gene', 'read_abundance')
cef_mock_annotated$read_abundance <- log2(cef_mock_annotated$read_abundance)
clinda_630_annotated <- aggregate(clinda_630_annotated$clinda_630_metaT_reads, by=list(clinda_630_annotated$gene), FUN=sum)
colnames(clinda_630_annotated) <- c('gene', 'read_abundance')
clinda_630_annotated$read_abundance <- log2(clinda_630_annotated$read_abundance)
clinda_mock_annotated <- aggregate(clinda_mock_annotated$clinda_mock_metaT_reads, by=list(clinda_mock_annotated$gene), FUN=sum)
colnames(clinda_mock_annotated) <- c('gene', 'read_abundance')
clinda_mock_annotated$read_abundance <- log2(clinda_mock_annotated$read_abundance)
strep_630_annotated <- aggregate(strep_630_annotated$strep_630_metaT_reads, by=list(strep_630_annotated$gene), FUN=sum)
colnames(strep_630_annotated) <- c('gene', 'read_abundance')
strep_630_annotated$read_abundance <- log2(strep_630_annotated$read_abundance)
strep_mock_annotated <- aggregate(strep_mock_annotated$strep_mock_metaT_reads, by=list(strep_mock_annotated$gene), FUN=sum)
colnames(strep_mock_annotated) <- c('gene', 'read_abundance')
strep_mock_annotated$read_abundance <- log2(strep_mock_annotated$read_abundance)

# Get the top expressed genes from each
noabx_top <- noabx_annotated[order(-noabx_annotated$read_abundance),][c(1:100),] 
cef_mock_top <- cef_630_annotated[order(-cef_630_annotated$read_abundance),][c(1:100),] 
cef_630_top <- cef_mock_annotated[order(-cef_mock_annotated$read_abundance),][c(1:100),] 
clinda_630_top <- clinda_630_annotated[order(-clinda_630_annotated$read_abundance),][c(1:100),] 
clinda_mock_top <- clinda_mock_annotated[order(-clinda_mock_annotated$read_abundance),][c(1:100),] 
strep_630_top <- strep_630_annotated[order(-strep_630_annotated$read_abundance),][c(1:100),] 
strep_mock_top <- strep_mock_annotated[order(-strep_mock_annotated$read_abundance),][c(1:100),] 
rm(noabx_annotated, cef_630_annotated, cef_mock_annotated, clinda_630_annotated, clinda_mock_annotated, strep_630_annotated, strep_mock_annotated)

#-------------------------------------------------------------------------------------------------------------------------------------#






#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()





