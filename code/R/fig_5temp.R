
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

# Output plot
plot_file <- 'results/figures/figure_5temp.pdf'

#--------------------------------------------------------------------------------------------------#

# Read in data
# Normalized Metatranscriptomes
cef_normalized_reads <- read.delim(cef_normalized_reads, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=7)
clinda_normalized_reads <- read.delim(clinda_normalized_reads, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=7)
strep_normalized_reads <- read.delim(strep_normalized_reads, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=7)

#--------------------------------------------------------------------------------------------------#

# Format data

# Remove C. difficile genes
cdiff_omit <- c('cdf','pdc','cdc','cdl','pdf')
cef_normalized_reads <- cef_normalized_reads[!cef_normalized_reads$organism %in% cdiff_omit,]
clinda_normalized_reads <- clinda_normalized_reads[!clinda_normalized_reads$organism %in% cdiff_omit,]
strep_normalized_reads <- strep_normalized_reads[!strep_normalized_reads$organism %in% cdiff_omit,]
rm(cdiff_omit)

# Screen for those genes that have a gene annotation
cef_annotated <- cef_normalized_reads[!rownames(cef_normalized_reads) %in% rownames(cef_normalized_reads[grep('unknown_\\d', rownames(cef_normalized_reads)),]), ]
clinda_annotated <- clinda_normalized_reads[!rownames(clinda_normalized_reads) %in% rownames(clinda_normalized_reads[grep('unknown_\\d', rownames(clinda_normalized_reads)),]), ]
strep_annotated <- strep_normalized_reads[!rownames(strep_normalized_reads) %in% rownames(strep_normalized_reads[grep('unknown_\\d', rownames(strep_normalized_reads)),]), ]
rm(cef_normalized_reads, clinda_normalized_reads, strep_normalized_reads)

# Screen for genes with no reads
cef_annotated[,2:3] <- round(cef_annotated[,2:3])
cef_annotated <- subset(cef_annotated, rowSums(cef_annotated[,2:3]) != 0)
clinda_annotated[,2:3] <- round(clinda_annotated[,2:3])
clinda_annotated <- subset(clinda_annotated, rowSums(clinda_annotated[,2:3]) != 0)
strep_annotated[,2:3] <- round(strep_annotated[,2:3])
strep_annotated <- subset(strep_annotated, rowSums(strep_annotated[,2:3]) != 0)

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

# Simplify column names
colnames(cef_annotated) <- c('gene', 'cef_630_reads', 'cef_mock_reads','organism','description','ko','pathways')
colnames(clinda_annotated) <- c('gene', 'clinda_630_reads', 'clinda_mock_reads','organism','description','ko','pathways')
colnames(strep_annotated) <- c('gene', 'strep_630_reads', 'strep_mock_reads','organism','description','ko','pathways')

#--------------------------------------------------------------------------------------------------#

# Write supplementary table panels
#table_s3a <- strep_annotated
#table_s3a$diff <- NULL
#colnames(table_s3a) <- c('KEGG_gene_hit','Normalized_cDNA_Reads(Infected)','Normalized_cDNA_Reads(Mock)')
#write.table(table_s3a, file='results/supplement/tables/Table_s3a.tsv', sep='\t', row.names=FALSE, quote=FALSE)
#table_s3b <- cef_annotated
#table_s3b$diff <- NULL
#colnames(table_s3b) <- c('KEGG_gene_hit','Normalized_cDNA_Reads(Infected)','Normalized_cDNA_Reads(Mock)')
#write.table(table_s3b, file='results/supplement/tables/Table_s3b.tsv', sep='\t', row.names=FALSE, quote=FALSE)
#table_s3c <- clinda_annotated
#table_s3c$diff <- NULL
#colnames(table_s3c) <- c('KEGG_gene_hit','Normalized_cDNA_Reads(Infected)','Normalized_cDNA_Reads(Mock)')
#write.table(table_s3c, file='results/supplement/tables/Table_s3c.tsv', sep='\t', row.names=FALSE, quote=FALSE)
#rm(table_s3a, table_s3b, table_s3c)

# Calculate and take top differences from top 10%
cef_10 <- round(nrow(cef_annotated) * 0.15)
clinda_10 <- round(nrow(clinda_annotated) * 0.15)
strep_10 <- round(nrow(strep_annotated) * 0.15)
cef_annotated$diff <- abs(cef_annotated$cef_630_reads - cef_annotated$cef_mock_reads)
cef_annotated <- cef_annotated[order(-cef_annotated$diff),]
top_cef_genes <- cef_annotated$gene[1:cef_10]
cef_annotated$diff <- NULL
clinda_annotated$diff <- abs(clinda_annotated$clinda_630_reads - clinda_annotated$clinda_mock_reads)
clinda_annotated <- clinda_annotated[order(-clinda_annotated$diff),]
top_clinda_genes <- clinda_annotated$gene[1:clinda_10]
clinda_annotated$diff <- NULL
strep_annotated$diff <- abs(strep_annotated$strep_630_reads - strep_annotated$strep_mock_reads)
strep_annotated <- strep_annotated[order(-strep_annotated$diff),]
top_strep_genes <- strep_annotated$gene[1:strep_10]
strep_annotated$diff <- NULL

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
top_cef_clinda_genes <- as.numeric(length(top_cef_clinda_genes))
top_cef_strep_genes <- as.numeric(length(top_cef_strep_genes))
top_clinda_strep_genes <- as.numeric(length(top_clinda_strep_genes))
top_all_genes <- as.numeric(length(top_genes))
top_cef_genes <- as.numeric(length(top_cef_genes))
top_clinda_genes <- as.numeric(length(top_clinda_genes))
top_strep_genes <- as.numeric(length(top_strep_genes))

# Merge treatment groups
final_annotated <- merge(strep_annotated, cef_annotated, by='gene', all=TRUE)
final_annotated <- merge(final_annotated, clinda_annotated, by='gene', all=TRUE)
final_annotated[is.na(final_annotated)] <- 0 # Undetectable genes set to 0
final_annotated <- subset(final_annotated, final_annotated$gene %in% top_genes)
rownames(final_annotated) <- final_annotated$description
final_annotated$gene <- NULL
final_annotated$organism.x <- NULL
final_annotated$description.x <- NULL
final_annotated$ko.x <- NULL
final_annotated$pathways.x <- NULL
final_annotated$organism.y <- NULL
final_annotated$description.y <- NULL
final_annotated$ko.y <- NULL
final_annotated$pathways.y <- NULL
final_annotated$organism <- NULL
final_annotated$description <- NULL
final_annotated$ko <- NULL
final_annotated$pathways <- NULL
final_annotated <- log2(final_annotated + 1)

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
top_cef_clinda_paths <- as.numeric(length(top_cef_clinda_paths))
top_cef_strep_paths <- as.numeric(length(top_cef_strep_paths))
top_clinda_strep_paths <- as.numeric(length(top_clinda_strep_paths))
top_all_paths <- as.numeric(length(top_paths))
top_cef_paths <- as.numeric(length(top_cef_paths))
top_clinda_paths <- as.numeric(length(top_clinda_paths))
top_strep_paths <- as.numeric(length(top_strep_paths))

# Merge treatment groups
final_pathways <- merge(cef_pathways, clinda_pathways, by='row.names', all=TRUE)
rownames(final_pathways) <- final_pathways$Row.names
final_pathways$Row.names <- NULL
final_pathways <- merge(final_pathways, strep_pathways, by='row.names', all=TRUE)
rownames(final_pathways) <- final_pathways$Row.names
final_pathways$Row.names <- NULL
final_pathways[is.na(final_pathways)] <- 0 # Undetectable genes set to 0
final_pathways <- subset(final_pathways, rownames(final_pathways) %in% top_paths)
final_pathways <- log2(final_pathways + 1)
rm(cef_annotated, clinda_annotated, strep_annotated, 
   cef_pathways, clinda_pathways, strep_pathways)

# Reformat names
# genes
rownames(final_annotated) <- gsub('_', ' ', rownames(final_annotated))
rownames(final_annotated) <- gsub('\\:cation symporter', '', rownames(final_annotated))
rownames(final_annotated) <- gsub('glycoside-pentoside-hexuronide', 'symporter', rownames(final_annotated))
rownames(final_annotated) <- capitalize(rownames(final_annotated))
# pathways
path_names <- as.vector(rownames(final_pathways))
for (x in c(1:length(path_names))) {
  path_names[x] <- strsplit(path_names[x], split=':')[[1]][1]
}
path_names <- gsub('_', ' ', path_names)
path_names <- gsub('Microbial metabolism in diverse environments', 'Metabolism in diverse environments', path_names)
path_names <- capitalize(path_names)
rownames(final_pathways) <- path_names
rm(path_names)

# Add pathway information to genes
rownames(final_annotated) <- c('(1) Dps family protein',
                               '(2,3,8) Formate acetyltransferase',            
                               '(4) Ribosome-associated factor Y',
                               '(5,8) Phosphopyruvate hydratase',         
                               '(6) HD superfamily hydrolase',
                               'ArsR family transcriptional regulator',
                               'GPH family symporter',
                               '(7) Protein grpE')

# Reverse row order for plotting
final_annotated <- final_annotated[nrow(final_annotated):1, ]
final_pathways <- final_pathways[nrow(final_pathways):1, ]

# Save variances
gene_var <- RowVar(final_annotated)
path_var <- RowVar(final_pathways)

#--------------------------------------------------------------------------------------------------#

# Generate bar plots

# Genes
pdf(file='results/figures/figure_5CD.pdf', width=12, height=8)
layout(matrix(c(1,2,
                1,2,
                3,3),
              nrow=3, ncol=2, byrow = TRUE))

par(mar=c(3,19.75,1,1), mgp=c(2.5, 0.75, 0), las=1, xaxs='i', yaxs='i')
barplot(as.matrix(t(final_annotated)), xaxt='n', xlim=c(0,14), ylim=c(0,57), beside=TRUE, horiz=TRUE, 
        xlab='', ylab='', angle=20, density=c(NA,20),  
        col=c(strep_col,strep_col,cef_col,cef_col,clinda_col,clinda_col), cex.names=1.4)
abline(h=seq(7.5,56,7), lty=5, col='black')
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2), cex.axis=1.2)
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.5, cex=0.9)
mtext('B', side=2, padj=-13, adj=15, font=2, cex=1.4)
mtext('A', side=2, padj=-13, adj=11, font=2, cex=1.4)

# Pathways
par(mar=c(3,19.75,1,1), mgp=c(2.5, 0.75, 0), las=1, xaxs='i', yaxs='i')
barplot(as.matrix(t(final_pathways)), xaxt='n', xlim=c(0,14), ylim=c(0,57), beside=TRUE, horiz=TRUE, 
        xlab='', ylab='', angle=20, density=c(NA,20), 
        col=c(strep_col,strep_col,cef_col,cef_col,clinda_col,clinda_col), cex.names=1.4)
abline(h=seq(7.5,56,7), lty=5, col='black')
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2), cex.axis=1.2)
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.5, cex=0.9)
mtext('D', side=2, padj=-13, adj=15, font=2, cex=1.4)
mtext('C', side=2, padj=-13, adj=11, font=2, cex=1.4)

# Legend plot
par(mar=c(0,0,0,0))
plot(0, xaxt='n', yaxt='n', bty='n', pch='', ylab='', xlab='', xlim=c(-10,10), ylim=c(-5,5))
legend(x=-9, y=3, legend=c('1. Metabolism of cofactors and vitamins',
                             '2. Pyruvate metabolism',
                             '3. Propanoate metabolism',
                             '4. RNA transport',
                             '5. Glycolysis / Gluconeogenesis',
                             '6. RNA degradation',
                             '7. Chaperones and folding catalysts',
                             '8. Microbial metabolism in diverse environments'),
       pt.cex=0, bty='n', cex=1.5)
legend(x=-1, y=3, legend=c(expression(italic('C. difficile')),'Mock-infected'), fill='black',
       density=c(NA, 20), pt.cex=2.3, cex=1.5)
legend(x=3, y=3.5, legend=c('Streptomycin-pretreated','Cefoperazone-pretreated','Clindamycin-pretreated'), 
       pt.bg=c(strep_col,cef_col,clinda_col), pch=22, pt.cex=2.3, cex=1.5)

dev.off()

#--------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()

