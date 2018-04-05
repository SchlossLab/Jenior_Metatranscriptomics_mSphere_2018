
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')

# Define files
# Normalized Metatranscriptomes
cef_normalized_reads <- 'data/read_mapping/cef_normalized_metaT.tsv'
clinda_normalized_reads <- 'data/read_mapping/clinda_normalized_metaT.tsv'
strep_normalized_reads <- 'data/read_mapping/strep_normalized_metaT.tsv'
noabx_normalized_reads <- 'data/read_mapping/noabx_normalized_metaT.tsv'

# Output plot
plot_file <- 'results/supplement/figures/figure_S5.pdf'
legend_file <- 'results/supplement/figures/figure_S5legend.pdf'

#--------------------------------------------------------------------------------------------------#

# Read in data
# Normalized Metatranscriptomes
cef_normalized_reads <- read.delim(cef_normalized_reads, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=7)
clinda_normalized_reads <- read.delim(clinda_normalized_reads, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=7)
strep_normalized_reads <- read.delim(strep_normalized_reads, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=7)
noabx_normalized_reads <- read.delim(noabx_normalized_reads, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=6)

#--------------------------------------------------------------------------------------------------#

# Format data

# Remove C. difficile genes
cdiff_omit <- c('cdf','pdc','cdc','cdl','pdf')
cef_normalized_reads <- cef_normalized_reads[!cef_normalized_reads$organism %in% cdiff_omit,]
clinda_normalized_reads <- clinda_normalized_reads[!clinda_normalized_reads$organism %in% cdiff_omit,]
strep_normalized_reads <- strep_normalized_reads[!strep_normalized_reads$organism %in% cdiff_omit,]
noabx_normalized_reads <- noabx_normalized_reads[!noabx_normalized_reads$organism %in% cdiff_omit,]
rm(cdiff_omit)

# Screen for those genes that have a gene annotation
cef_annotated <- cef_normalized_reads[!rownames(cef_normalized_reads) %in% rownames(cef_normalized_reads[grep('unknown_\\d', rownames(cef_normalized_reads)),]), ]
clinda_annotated <- clinda_normalized_reads[!rownames(clinda_normalized_reads) %in% rownames(clinda_normalized_reads[grep('unknown_\\d', rownames(clinda_normalized_reads)),]), ]
strep_annotated <- strep_normalized_reads[!rownames(strep_normalized_reads) %in% rownames(strep_normalized_reads[grep('unknown_\\d', rownames(strep_normalized_reads)),]), ]
noabx_annotated <- noabx_normalized_reads[!rownames(noabx_normalized_reads) %in% rownames(noabx_normalized_reads[grep('unknown_\\d', rownames(noabx_normalized_reads)),]), ]
rm(cef_normalized_reads, clinda_normalized_reads, strep_normalized_reads, noabx_normalized_reads)

# Screen for genes with no reads
cef_annotated[,2:3] <- round(cef_annotated[,2:3])
cef_annotated <- subset(cef_annotated, rowSums(cef_annotated[,2:3]) != 0)
clinda_annotated[,2:3] <- round(clinda_annotated[,2:3])
clinda_annotated <- subset(clinda_annotated, rowSums(clinda_annotated[,2:3]) != 0)
strep_annotated[,2:3] <- round(strep_annotated[,2:3])
strep_annotated <- subset(strep_annotated, rowSums(strep_annotated[,2:3]) != 0)
noabx_annotated[,2] <- round(noabx_annotated[,2])
noabx_annotated <- subset(noabx_annotated, noabx_annotated[,2] != 0)

# Screen out ribosomal genes + hypothetical, putative, and uncharacterized annotations
cef_annotated <- subset(cef_annotated, !grepl('Ribosomal_RNA*', cef_annotated$description))
cef_annotated <- subset(cef_annotated, !grepl('ribosomal_RNA*', cef_annotated$description))
cef_annotated <- subset(cef_annotated, !grepl('*ribosomal_RNA*', cef_annotated$description))
cef_annotated <- subset(cef_annotated, !grepl('*ribosomal_protein*', cef_annotated$description))
cef_annotated <- subset(cef_annotated, !grepl('16S_rRNA_*', cef_annotated$description))
cef_annotated <- subset(cef_annotated, !grepl('23S_rRNA_*', cef_annotated$description))
cef_annotated <- subset(cef_annotated, !grepl('uncharacterized_*', cef_annotated$description))
cef_annotated <- subset(cef_annotated, !grepl('putative_*', cef_annotated$description))
cef_annotated <- subset(cef_annotated, !grepl('Predicted_*', cef_annotated$description))
cef_annotated <- subset(cef_annotated, description != 'hypothetical_protein')
clinda_annotated <- subset(clinda_annotated, !grepl('Ribosomal_RNA*', clinda_annotated$description))
clinda_annotated <- subset(clinda_annotated, !grepl('ribosomal_RNA*', clinda_annotated$description))
clinda_annotated <- subset(clinda_annotated, !grepl('*ribosomal_RNA*', clinda_annotated$description))
clinda_annotated <- subset(clinda_annotated, !grepl('*ribosomal_protein*', clinda_annotated$description))
clinda_annotated <- subset(clinda_annotated, !grepl('16S_rRNA_*', clinda_annotated$description))
clinda_annotated <- subset(clinda_annotated, !grepl('23S_rRNA_*', clinda_annotated$description))
clinda_annotated <- subset(clinda_annotated, !grepl('uncharacterized_*', clinda_annotated$description))
clinda_annotated <- subset(clinda_annotated, !grepl('putative_*', clinda_annotated$description))
clinda_annotated <- subset(clinda_annotated, !grepl('Predicted_*', clinda_annotated$description))
clinda_annotated <- subset(clinda_annotated, description != 'hypothetical_protein')
strep_annotated <- subset(strep_annotated, !grepl('Ribosomal_RNA*', strep_annotated$description))
strep_annotated <- subset(strep_annotated, !grepl('ribosomal_RNA*', strep_annotated$description))
strep_annotated <- subset(strep_annotated, !grepl('*ribosomal_RNA*', strep_annotated$description))
strep_annotated <- subset(strep_annotated, !grepl('*ribosomal_protein*', strep_annotated$description))
strep_annotated <- subset(strep_annotated, !grepl('16S_rRNA_*', strep_annotated$description))
strep_annotated <- subset(strep_annotated, !grepl('23S_rRNA_*', strep_annotated$description))
strep_annotated <- subset(strep_annotated, !grepl('uncharacterized_*', strep_annotated$description))
strep_annotated <- subset(strep_annotated, !grepl('putative_*', strep_annotated$description))
strep_annotated <- subset(strep_annotated, !grepl('Predicted_*', strep_annotated$description))
strep_annotated <- subset(strep_annotated, description != 'hypothetical_protein')
noabx_annotated <- subset(noabx_annotated, !grepl('Ribosomal_RNA*', noabx_annotated$description))
noabx_annotated <- subset(noabx_annotated, !grepl('ribosomal_RNA*', noabx_annotated$description))
noabx_annotated <- subset(noabx_annotated, !grepl('*ribosomal_RNA*', noabx_annotated$description))
noabx_annotated <- subset(noabx_annotated, !grepl('*ribosomal_protein*', noabx_annotated$description))
noabx_annotated <- subset(noabx_annotated, !grepl('16S_rRNA_*', noabx_annotated$description))
noabx_annotated <- subset(noabx_annotated, !grepl('23S_rRNA_*', noabx_annotated$description))
noabx_annotated <- subset(noabx_annotated, !grepl('uncharacterized_*', noabx_annotated$description))
noabx_annotated <- subset(noabx_annotated, !grepl('putative_*', noabx_annotated$description))
noabx_annotated <- subset(noabx_annotated, !grepl('Predicted_*', noabx_annotated$description))
noabx_annotated <- subset(noabx_annotated, description != 'hypothetical_protein')

# Simplify column names
colnames(cef_annotated) <- c('gene', 'cef_630_reads', 'cef_mock_reads','organism','description','ko','pathways')
colnames(clinda_annotated) <- c('gene', 'clinda_630_reads', 'clinda_mock_reads','organism','description','ko','pathways')
colnames(strep_annotated) <- c('gene', 'strep_630_reads', 'strep_mock_reads','organism','description','ko','pathways')
colnames(noabx_annotated) <- c('gene', 'noabx_mock_reads','organism','description','ko','pathways')

#--------------------------------------------------------------------------------------------------#

# Aggregate genes with same annotation
strep_genes <- aggregate(cbind(strep_annotated$strep_630_reads, strep_annotated$strep_mock_reads), by=list(strep_annotated$description), FUN=sum)
colnames(strep_genes) <- c('description', 'strep_630_reads', 'strep_mock_reads')
cef_genes <- aggregate(cbind(cef_annotated$cef_630_reads, cef_annotated$cef_mock_reads), by=list(cef_annotated$description), FUN=sum)
colnames(cef_genes) <- c('description', 'cef_630_reads', 'cef_mock_reads')
clinda_genes <- aggregate(cbind(clinda_annotated$clinda_630_reads, clinda_annotated$clinda_mock_reads), by=list(clinda_annotated$description), FUN=sum)
colnames(clinda_genes) <- c('description', 'clinda_630_reads', 'clinda_mock_reads')
noabx_genes <- aggregate(noabx_annotated$noabx_mock_reads, by=list(noabx_annotated$description), FUN=sum)
colnames(noabx_genes) <- c('description', 'noabx_mock_reads')

# Calculate and take top differences from top 15%
cef_10 <- round(nrow(cef_genes) * 0.15)
clinda_10 <- round(nrow(clinda_genes) * 0.15)
strep_10 <- round(nrow(strep_genes) * 0.15)
cef_genes$infect_diff <- abs(cef_genes$cef_630_reads - cef_genes$cef_mock_reads)
cef_genes <- cef_genes[order(-cef_genes$infect_diff),]
top_cef_genes <- cef_genes$description[1:cef_10]
cef_genes$infect_diff <- NULL
clinda_genes$infect_diff <- abs(clinda_genes$clinda_630_reads - clinda_genes$clinda_mock_reads)
clinda_genes <- clinda_genes[order(-clinda_genes$infect_diff),]
top_clinda_genes <- clinda_genes$description[1:clinda_10]
clinda_genes$infect_diff <- NULL
strep_genes$infect_diff <- abs(strep_genes$strep_630_reads - strep_genes$strep_mock_reads)
strep_genes <- strep_genes[order(-strep_genes$infect_diff),]
top_strep_genes <- strep_genes$description[1:strep_10]
strep_genes$infect_diff <- NULL

# Calculate overlap of top differences
top_cef_clinda_genes <- intersect(top_cef_genes, top_clinda_genes)
top_cef_strep_genes <- intersect(top_cef_genes, top_strep_genes)
top_clinda_strep_genes <- intersect(top_clinda_genes, top_strep_genes)
top_genes <- intersect(top_cef_genes, intersect(top_clinda_genes, top_strep_genes))
top_cef_genes <- subset(top_cef_genes, !top_cef_genes %in% top_cef_clinda_genes)
top_cef_genes <- subset(top_cef_genes, !top_cef_genes %in% top_cef_strep_genes)
top_cef_genes <- subset(top_cef_genes, !top_cef_genes %in% top_genes)
top_clinda_genes <- subset(top_clinda_genes, !top_clinda_genes %in% top_cef_clinda_genes)
top_clinda_genes <- subset(top_clinda_genes, !top_clinda_genes %in% top_clinda_strep_genes)
top_clinda_genes <- subset(top_clinda_genes, !top_clinda_genes %in% top_genes)
top_strep_genes <- subset(top_strep_genes, !top_strep_genes %in% top_clinda_strep_genes)
top_strep_genes <- subset(top_strep_genes, !top_strep_genes %in% top_cef_strep_genes)
top_strep_genes <- subset(top_strep_genes, !top_strep_genes %in% top_genes)
rm(top_cef_clinda_genes, top_cef_strep_genes, top_clinda_strep_genes, top_genes)

# Get individual gene sets
top_strep_genes <- subset(strep_genes, strep_genes$description %in% top_strep_genes)
top_cef_genes <- subset(cef_genes, cef_genes$description %in% top_cef_genes)
top_clinda_genes <- subset(clinda_genes, clinda_genes$description %in% top_clinda_genes)
rm(strep_genes, cef_genes, clinda_genes)

# Merge with untreated group
top_strep_genes <- merge(top_strep_genes, noabx_genes, by='description')
rownames(top_strep_genes) <- top_strep_genes$description
top_strep_genes$description <- NULL
top_cef_genes <- merge(top_cef_genes, noabx_genes, by='description')
rownames(top_cef_genes) <- top_cef_genes$description
top_cef_genes$description <- NULL
top_clinda_genes <- merge(top_clinda_genes, noabx_genes, by='description')
rownames(top_clinda_genes) <- top_clinda_genes$description
top_clinda_genes$description <- NULL
rm(noabx_genes)

# Pick top 10 differences between avg. abx-treated and untreated
top_strep_genes$infect_mean <- rowMeans(top_strep_genes[,2:3])
top_strep_genes$infect_diff <- abs(top_strep_genes$infect_mean - top_strep_genes$noabx_mock_reads)
top_strep_genes <- top_strep_genes[order(-top_strep_genes$infect_diff),]
top_strep_genes <- top_strep_genes[1:11,]
top_strep_genes <- subset(top_strep_genes, rownames(top_strep_genes) != 'type_I')
top_strep_genes$infect_diff <- NULL
top_strep_genes$infect_mean <- NULL
top_cef_genes$infect_mean <- rowMeans(top_cef_genes[,2:3])
top_cef_genes$infect_diff <- abs(top_cef_genes$infect_mean - top_cef_genes$noabx_mock_reads)
top_cef_genes <- top_cef_genes[order(-top_cef_genes$infect_diff),]
top_cef_genes <- top_cef_genes[1:10,]
top_cef_genes$infect_diff <- NULL
top_cef_genes$infect_mean <- NULL
top_clinda_genes$infect_mean <- rowMeans(top_clinda_genes[,2:3])
top_clinda_genes$infect_diff <- abs(top_clinda_genes$infect_mean - top_clinda_genes$noabx_mock_reads)
top_clinda_genes <- top_clinda_genes[order(-top_clinda_genes$infect_diff),]
top_clinda_genes <- top_clinda_genes[1:11,]
top_clinda_genes <- subset(top_clinda_genes, rownames(top_clinda_genes) != 'homology_to_Homo_sapiens')
top_clinda_genes$infect_diff <- NULL
top_clinda_genes$infect_mean <- NULL

# Transform data and format names
top_strep_genes <- log2(top_strep_genes + 1)
rownames(top_strep_genes) <- gsub('_', ' ', rownames(top_strep_genes))
rownames(top_strep_genes) <- gsub(' \\(EC\\:2\\.7\\.9\\.1)', '', rownames(top_strep_genes))
rownames(top_strep_genes) <- gsub('binding-protein-dependent transport system inner membrane protein', 'transport system inner membrane protein', rownames(top_strep_genes))
rownames(top_strep_genes) <- capitalize(rownames(top_strep_genes))
top_cef_genes <- log2(top_cef_genes + 1)
rownames(top_cef_genes) <- gsub('_', ' ', rownames(top_cef_genes))
rownames(top_cef_genes) <- gsub(' \\(EF\\-1A\\/EF\\-Tu\\) \\(EC\\:3\\.6\\.5\\.3\\)', '', rownames(top_cef_genes))
rownames(top_cef_genes) <- gsub('two-component system sensor histidine kinase/response regulator', 'two-component system histidine kinase regulator', rownames(top_cef_genes))
rownames(top_cef_genes) <- capitalize(rownames(top_cef_genes))
top_clinda_genes <- log2(top_clinda_genes + 1)
rownames(top_clinda_genes) <- gsub('_', ' ', rownames(top_clinda_genes))
rownames(top_clinda_genes) <- gsub(' \\(EC\\:3\\.2\\.1\\.23)', '', rownames(top_clinda_genes))
rownames(top_clinda_genes) <- capitalize(rownames(top_clinda_genes))

# Reverse row order for plotting
top_strep_genes <- top_strep_genes[nrow(top_strep_genes):1, ]
top_cef_genes <- top_cef_genes[nrow(top_cef_genes):1, ]
top_clinda_genes <- top_clinda_genes[nrow(top_clinda_genes):1, ]

# Transpose
top_strep_genes <- as.matrix(t(top_strep_genes))
top_cef_genes <- as.matrix(t(top_cef_genes))
top_clinda_genes <- as.matrix(t(top_clinda_genes))

#-------------#

# Aggregate pathways
cef_pathways <- subset(cef_annotated, !is.na(cef_annotated$pathways))
cef_pathways <- aggregate(cbind(cef_pathways$cef_630_reads, cef_pathways$cef_mock_reads), by=list(cef_pathways$pathways), FUN=sum)
rownames(cef_pathways) <- cef_pathways$Group.1
cef_pathways$Group.1 <- NULL
cef_pathways <- subset(cef_pathways, !rownames(cef_pathways) %in% c('Ribosome:Translation','Metabolic_pathways'))
colnames(cef_pathways) <- c('cef_630_reads', 'cef_mock_reads')
clinda_pathways <- subset(clinda_annotated, !is.na(clinda_annotated$pathways))
clinda_pathways <- aggregate(cbind(clinda_pathways$clinda_630_reads, clinda_pathways$clinda_mock_reads), by=list(clinda_pathways$pathways), FUN=sum)
rownames(clinda_pathways) <- clinda_pathways$Group.1
clinda_pathways$Group.1 <- NULL
clinda_pathways <- subset(clinda_pathways, !rownames(clinda_pathways) %in% c('Ribosome:Translation','Metabolic_pathways'))
colnames(clinda_pathways) <- c('clinda_630_reads', 'clinda_mock_reads')
strep_pathways <- subset(strep_annotated, !is.na(strep_annotated$pathways))
strep_pathways <- aggregate(cbind(strep_pathways$strep_630_reads, strep_pathways$strep_mock_reads), by=list(strep_pathways$pathways), FUN=sum)
rownames(strep_pathways) <- strep_pathways$Group.1
strep_pathways$Group.1 <- NULL
strep_pathways <- subset(strep_pathways, !rownames(strep_pathways) %in% c('Ribosome:Translation','Metabolic_pathways'))
colnames(strep_pathways) <- c('strep_630_reads', 'strep_mock_reads')
noabx_pathways <- subset(noabx_annotated, !is.na(noabx_annotated$pathways))
noabx_pathways <- aggregate(noabx_pathways$noabx_mock_reads, by=list(noabx_pathways$pathways), FUN=sum)
rownames(noabx_pathways) <- noabx_pathways$Group.1
noabx_pathways$Group.1 <- NULL
noabx_pathways <- subset(noabx_pathways, !rownames(noabx_pathways) %in% c('Ribosome:Translation','Metabolic_pathways'))
colnames(noabx_pathways) <- 'noabx_mock_reads'
rm(cef_annotated, clinda_annotated, strep_annotated, noabx_annotated)

# Calculate and take top differences
cef_10 <- round(nrow(cef_pathways) * 0.15)
clinda_10 <- round(nrow(clinda_pathways) * 0.15)
strep_10 <- round(nrow(strep_pathways) * 0.15)
cef_pathways$diff <- abs(cef_pathways$cef_630_reads - cef_pathways$cef_mock_reads)
cef_pathways <- cef_pathways[order(-cef_pathways$diff),]
top_cef_paths <- rownames(cef_pathways[1:cef_10,])
cef_pathways$diff <- NULL
clinda_pathways$diff <- abs(clinda_pathways$clinda_630_reads - clinda_pathways$clinda_mock_reads)
clinda_pathways <- clinda_pathways[order(-clinda_pathways$diff),]
top_clinda_paths <- rownames(clinda_pathways[1:clinda_10,])
clinda_pathways$diff <- NULL
strep_pathways$diff <- abs(strep_pathways$strep_630_reads - strep_pathways$strep_mock_reads)
strep_pathways <- strep_pathways[order(-strep_pathways$diff),]
top_strep_paths <- rownames(strep_pathways[1:strep_10,])
strep_pathways$diff <- NULL

# Calculate overlap of top differences
top_cef_clinda_paths <- intersect(top_cef_paths, top_clinda_paths)
top_cef_strep_paths <- intersect(top_cef_paths, top_strep_paths)
top_clinda_strep_paths <- intersect(top_clinda_paths, top_strep_paths)
top_paths <- intersect(top_cef_paths, intersect(top_clinda_paths, top_strep_paths))
top_cef_paths <- subset(top_cef_paths, !top_cef_paths %in% top_cef_clinda_paths)
top_cef_paths <- subset(top_cef_paths, !top_cef_paths %in% top_cef_strep_paths)
top_cef_paths <- subset(top_cef_paths, !top_cef_paths %in% top_paths)
top_clinda_paths <- subset(top_clinda_paths, !top_clinda_paths %in% top_cef_clinda_paths)
top_clinda_paths <- subset(top_clinda_paths, !top_clinda_paths %in% top_clinda_strep_paths)
top_clinda_paths <- subset(top_clinda_paths, !top_clinda_paths %in% top_paths)
top_strep_paths <- subset(top_strep_paths, !top_strep_paths %in% top_clinda_strep_paths)
top_strep_paths <- subset(top_strep_paths, !top_strep_paths %in% top_cef_strep_paths)
top_strep_paths <- subset(top_strep_paths, !top_strep_paths %in% top_paths)
rm(top_cef_clinda_paths,top_cef_strep_paths,top_clinda_strep_paths,top_paths)

# Get individual gene sets
top_strep_paths <- subset(strep_pathways, rownames(strep_pathways) %in% top_strep_paths)
top_cef_paths <- subset(cef_pathways, rownames(cef_pathways) %in% top_cef_paths)
top_clinda_paths <- subset(clinda_pathways, rownames(clinda_pathways) %in% top_clinda_paths)
rm(strep_pathways, cef_pathways, clinda_pathways)

# Merge with untreated group
top_strep_paths <- merge(top_strep_paths, noabx_pathways, by='row.names')
rownames(top_strep_paths) <- top_strep_paths$Row.names
top_strep_paths$Row.names <- NULL
top_cef_paths <- merge(top_cef_paths, noabx_pathways, by='row.names')
rownames(top_cef_paths) <- top_cef_paths$Row.names
top_cef_paths$Row.names <- NULL
top_clinda_paths <- merge(top_clinda_paths, noabx_pathways, by='row.names')
rownames(top_clinda_paths) <- top_clinda_paths$Row.names
top_clinda_paths$Row.names <- NULL
rm(noabx_pathways)

# Transform data and format names
top_strep_paths <- log2(top_strep_paths + 1)
rownames(top_strep_paths) <- strsplit(rownames(top_strep_paths), split=':')[[1]][1]
rownames(top_strep_paths) <- gsub('_', ' ', rownames(top_strep_paths))
rownames(top_strep_paths) <- capitalize(rownames(top_strep_paths))
top_cef_paths <- log2(top_cef_paths + 1)
rownames(top_cef_paths) <- strsplit(rownames(top_cef_paths), split=':')[[1]][1]
rownames(top_cef_paths) <- gsub('_', ' ', rownames(top_cef_paths))
rownames(top_cef_paths) <- capitalize(rownames(top_cef_paths))
path_names <- as.vector(rownames(top_cef_paths))
for (x in c(1:length(path_names))) {
  path_names[x] <- strsplit(path_names[x], split=':')[[1]][1]
}
rownames(top_cef_paths) <- path_names
rm(path_names)
rownames(top_cef_paths) <- gsub('_', ' ', rownames(top_cef_paths))
rownames(top_cef_paths) <- capitalize(rownames(top_cef_paths))
top_clinda_paths <- log2(top_clinda_paths + 1)
path_names <- as.vector(rownames(top_clinda_paths))
for (x in c(1:length(path_names))) {
  path_names[x] <- strsplit(path_names[x], split=':')[[1]][1]
}
rownames(top_clinda_paths) <- path_names
rm(path_names)
rownames(top_clinda_paths) <- gsub('_', ' ', rownames(top_clinda_paths))
rownames(top_clinda_paths) <- capitalize(rownames(top_clinda_paths))

# Reverse row order for plotting
top_cef_paths <- top_cef_paths[nrow(top_cef_paths):1, ]
top_clinda_paths <- top_clinda_paths[nrow(top_clinda_paths):1, ]

# Transpose
top_strep_paths <- as.matrix(t(top_strep_paths))
top_cef_paths <- as.matrix(t(top_cef_paths))
top_clinda_paths <- as.matrix(t(top_clinda_paths))

# Combine for plot
top_paths <- cbind(top_clinda_paths, top_cef_paths, top_strep_paths)
# Strep = 1
# Cef = 7 - 
# Clinda = 2
rm(top_strep_paths, top_cef_paths, top_clinda_paths)

#--------------------------------------------------------------------------------------------------#

# Generate bar plots
pdf(file=plot_file, width=17, height=8)
layout(matrix(c(1,2,
                3,4),
              nrow=2, ncol=2, byrow=TRUE))

# Genes
par(mar=c(3,21,2,1), mgp=c(2.5, 0.75, 0), las=1, xaxs='i', yaxs='i')
# Strep
barplot(top_strep_genes, xaxt='n', xlim=c(0,14), ylim=c(0,41), beside=TRUE, horiz=TRUE, 
        xlab='', ylab='', angle=20, density=c(NA,20,NA),  
        col=c(strep_col,strep_col,noabx_col))
abline(h=seq(4.5,40,4), lty=5, col='black')
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2))
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.5, cex=0.9)
mtext('Streptomycin-pretreated (Genes)', side=3)
mtext('A', side=2, padj=-8, adj=16, font=2, cex=1.6)
# Cef
barplot(top_cef_genes, xaxt='n', xlim=c(0,14), ylim=c(0,41), beside=TRUE, horiz=TRUE, 
        xlab='', ylab='', angle=20, density=c(NA,20,NA),  
        col=c(cef_col,cef_col,noabx_col))
abline(h=seq(4.5,40,4), lty=5, col='black')
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2))
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.5, cex=0.9)
mtext('Cefoperazone-pretreated (Genes)', side=3)
mtext('B', side=2, padj=-8, adj=16, font=2, cex=1.6)
# Clinda
barplot(top_clinda_genes, xaxt='n', xlim=c(0,14), ylim=c(0,41), beside=TRUE, horiz=TRUE, 
        xlab='', ylab='', angle=20, density=c(NA,20,NA),  
        col=c(clinda_col,clinda_col,noabx_col))
abline(h=seq(4.5,40,4), lty=5, col='black')
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2))
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.5, cex=0.9)
mtext('Clindamycin-pretreated (Genes)', side=3)
mtext('C', side=2, padj=-8, adj=16, font=2, cex=1.6)

# Pathways
barplot(top_paths, xaxt='n', xlim=c(0,14), ylim=c(0,41), beside=TRUE, horiz=TRUE, 
        xlab='', ylab='', angle=20, density=c(NA,20,NA), 
        col=c(clinda_col,clinda_col,noabx_col,
              clinda_col,clinda_col,noabx_col,
              cef_col,cef_col,noabx_col,
              cef_col,cef_col,noabx_col,
              cef_col,cef_col,noabx_col,
              cef_col,cef_col,noabx_col,
              cef_col,cef_col,noabx_col,
              cef_col,cef_col,noabx_col,
              cef_col,cef_col,noabx_col,
              strep_col,strep_col,noabx_col))
abline(h=36.5, lwd=2, col='black')
abline(h=8.5, lwd=2, col='black')
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2), cex.axis=1.2)
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.5, cex=0.9)
mtext('All Pretreatments (Pathways)', side=3)
mtext('D', side=2, padj=-8, adj=16, font=2, cex=1.6)

dev.off()

# Legends
pdf(file=legend_file, width=10, height=5)
par(mar=c(0,0,0,0))
plot(0, xaxt='n', yaxt='n', bty='n', pch='', ylab='', xlab='', xlim=c(-10,10), ylim=c(-5,5))
legend(x=0, y=2, legend=c('Mock-infected',as.expression(bquote(paste(italic('C. difficile'),'-infected')))), 
       fill='black', density=c(NA, 20), pt.cex=2.3, cex=1.5)
legend(x=0, y=0, legend=c('No antibiotics','Streptomycin-pretreated','Cefoperazone-pretreated','Clindamycin-pretreated'), 
       pt.bg=c(noabx_col,strep_col,cef_col,clinda_col), pch=22, pt.cex=2.3, cex=1.5)
dev.off()

#--------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only=TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()

