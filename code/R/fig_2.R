
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define files
# Normalized Metatranscriptomes
cef_normalized_reads <- 'data/read_mapping/cef_normalized_metaT.tsv'
clinda_normalized_reads <- 'data/read_mapping/clinda_normalized_metaT.tsv'
strep_normalized_reads <- 'data/read_mapping/strep_normalized_metaT.tsv'

# KEGG pathway annotations for genes
cef_pathways <- 'data/read_mapping/pathway_expression/cef_pathways.tsv'
clinda_pathways <- 'data/read_mapping/pathway_expression/clinda_pathways.tsv'
strep_pathways <- 'data/read_mapping/pathway_expression/strep_pathways.tsv'

# Gene look up table
genes <- 'data/gene_names.tsv'

# Output plot
plot_file <- 'results/figures/figure_2.pdf'

#--------------------------------------------------------------------------------------------------#

# Read in data
# Normalized Metatranscriptomes
cef_normalized_reads <- read.delim(cef_normalized_reads, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=7)
clinda_normalized_reads <- read.delim(clinda_normalized_reads, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=7)
strep_normalized_reads <- read.delim(strep_normalized_reads, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=7)

# Pooled pathway mappings
cef_pathways <- read.delim(cef_pathways, sep='\t', header=TRUE, row.names=3)
clinda_pathways <- read.delim(clinda_pathways, sep='\t', header=TRUE, row.names=3)
strep_pathways <- read.delim(strep_pathways, sep='\t', header=TRUE, row.names=3)

# Gene functions
genes <- read.delim(genes, sep='\t', header=TRUE, stringsAsFactors=FALSE)

#--------------------------------------------------------------------------------------------------#

# Format data

# Screen for those genes that have a gene annotation
cef_annotated <- cef_normalized_reads[!rownames(cef_normalized_reads) %in% rownames(cef_normalized_reads[grep('unknown_\\d', rownames(cef_normalized_reads)),]), ]
clinda_annotated <- clinda_normalized_reads[!rownames(clinda_normalized_reads) %in% rownames(clinda_normalized_reads[grep('unknown_\\d', rownames(clinda_normalized_reads)),]), ]
strep_annotated <- strep_normalized_reads[!rownames(strep_normalized_reads) %in% rownames(strep_normalized_reads[grep('unknown_\\d', rownames(strep_normalized_reads)),]), ]
rm(cef_normalized_reads, clinda_normalized_reads, strep_normalized_reads)

# Screen out ribosomal genes
cef_annotated <- subset(cef_annotated, !grepl('Ribosomal_RNA*', cef_annotated$description))
cef_annotated <- subset(cef_annotated, !grepl('ribosomal_RNA*', cef_annotated$description))
cef_annotated <- subset(cef_annotated, !grepl('*ribosomal_RNA*', cef_annotated$description))
clinda_annotated <- subset(clinda_annotated, !grepl('Ribosomal_RNA*', clinda_annotated$description))
clinda_annotated <- subset(clinda_annotated, !grepl('ribosomal_RNA*', clinda_annotated$description))
clinda_annotated <- subset(clinda_annotated, !grepl('*ribosomal_RNA*', clinda_annotated$description))
strep_annotated <- subset(strep_annotated, !grepl('Ribosomal_RNA*', strep_annotated$description))
strep_annotated <- subset(strep_annotated, !grepl('ribosomal_RNA*', strep_annotated$description))
strep_annotated <- subset(strep_annotated, !grepl('*ribosomal_RNA*', strep_annotated$description))

# Remove C. difficile genes
cdiff_omit <- c('cdf','pdc','cdc','cdl','pdf')
cef_annotated <- cef_annotated[!cef_annotated$organism %in% cdiff_omit,]
clinda_annotated <- clinda_annotated[!clinda_annotated$organism %in% cdiff_omit,]
strep_annotated <- strep_annotated[!strep_annotated$organism %in% cdiff_omit,]
rm(cdiff_omit)

# Remove hypothetical and uncharacterized annotations
cef_annotated <- subset(cef_annotated, description != 'hypothetical_protein')
cef_annotated <- subset(cef_annotated, description != 'uncharacterized_*')
clinda_annotated <- subset(clinda_annotated, description != 'hypothetical_protein')
clinda_annotated <- subset(clinda_annotated, description != 'uncharacterized_*')
strep_annotated <- subset(strep_annotated, description != 'hypothetical_protein')
strep_annotated <- subset(strep_annotated, description != 'uncharacterized_*')

# Save pathway information
cef_pathways <- cef_annotated[complete.cases(cef_annotated),]
clinda_pathways <- clinda_annotated[complete.cases(clinda_annotated),]
strep_pathways <- strep_annotated[complete.cases(strep_annotated),]

# Aggregate identical genes, regardless of organism
cef_annotated <- aggregate(cbind(cef_annotated$cef_630_metaT_reads, cef_annotated$cef_mock_metaT_reads), by=list(cef_annotated$description), FUN=sum)
colnames(cef_annotated) <- c('gene', 'cef_630_reads', 'cef_mock_reads')
clinda_annotated <- aggregate(cbind(clinda_annotated$clinda_630_metaT_reads, clinda_annotated$clinda_mock_metaT_reads), by=list(clinda_annotated$description), FUN=sum)
colnames(clinda_annotated) <- c('gene', 'clinda_630_reads', 'clinda_mock_reads')
strep_annotated <- aggregate(cbind(strep_annotated$strep_630_metaT_reads, strep_annotated$strep_mock_metaT_reads), by=list(strep_annotated$description), FUN=sum)
colnames(strep_annotated) <- c('gene', 'strep_630_reads', 'strep_mock_reads')

# Aggregate pathways
cef_pathways <- aggregate(cbind(cef_pathways$cef_630_metaT_reads, cef_pathways$cef_mock_metaT_reads), by=list(cef_pathways$pathways), FUN=sum)
colnames(cef_pathways) <- c('pathway', 'cef_630_reads', 'cef_mock_reads')
clinda_pathways <- aggregate(cbind(clinda_pathways$clinda_630_metaT_reads, clinda_pathways$clinda_mock_metaT_reads), by=list(clinda_pathways$pathways), FUN=sum)
colnames(clinda_pathways) <- c('pathway', 'clinda_630_reads', 'clinda_mock_reads')
strep_pathways <- aggregate(cbind(strep_pathways$strep_630_metaT_reads, strep_pathways$strep_mock_metaT_reads), by=list(strep_pathways$pathways), FUN=sum)
colnames(strep_pathways) <- c('pathway', 'strep_630_reads', 'strep_mock_reads')

#--------------------------------------------------------------------------------------------------#

# Calculate ratios of infected to mock
cef_annotated$infected_v_mock <- cef_annotated$cef_630_reads / cef_annotated$cef_mock_reads
cef_annotated$mock_v_infected <- -(cef_annotated$cef_mock_reads / cef_annotated$cef_630_reads)
clinda_annotated$infected_v_mock <- clinda_annotated$clinda_630_reads / clinda_annotated$clinda_mock_reads
clinda_annotated$mock_v_infected <- -(clinda_annotated$clinda_mock_reads / clinda_annotated$clinda_630_reads)
strep_annotated$infected_v_mock <- strep_annotated$strep_630_reads / strep_annotated$strep_mock_reads
strep_annotated$mock_v_infected <- -(strep_annotated$strep_mock_reads / strep_annotated$strep_630_reads)

# Calculate degree of difference between groups
cef_annotated$change <- cef_annotated$infected_v_mock - cef_annotated$mock_v_infected
clinda_annotated$change <- clinda_annotated$infected_v_mock - clinda_annotated$mock_v_infected
strep_annotated$change <- strep_annotated$infected_v_mock - strep_annotated$mock_v_infected




# Rank differences and subset to top 15
cef_annotated <- cef_annotated[order(-cef_annotated$change),]
cef_annotated_top <- cef_annotated[1:15,]
clinda_annotated <- clinda_annotated[order(-clinda_annotated$change),]
clinda_annotated_top <- clinda_annotated[1:15,]
strep_annotated <- strep_annotated[order(-strep_annotated$change),]
strep_annotated_top <- strep_annotated[1:15,]

# Log transform expression
cef_annotated_top[,2:3] <- log2(cef_annotated_top[,2:3])
clinda_annotated_top[,2:3] <- log2(clinda_annotated_top[,2:3])
strep_annotated_top[,2:3] <- log2(strep_annotated_top[,2:3])















# Reassociate with pathway annotations
cef_630_top_pathways <- merge(cef_630_top, all_pathways, by='gene', all.x=TRUE)
cef_mock_top_pathways <- merge(cef_mock_top, all_pathways, by='gene', all.x=TRUE)
clinda_630_top_pathways <- merge(clinda_630_top, all_pathways, by='gene', all.x=TRUE)
clinda_mock_top_pathways <- merge(clinda_mock_top, all_pathways, by='gene', all.x=TRUE)
strep_630_top_pathways <- merge(strep_630_top, all_pathways, by='gene', all.x=TRUE)
strep_mock_top_pathways <- merge(strep_mock_top, all_pathways, by='gene', all.x=TRUE)
rm(all_pathways)

# Add rownames
cef_630_top <- merge(cef_630_top, genes, by='gene', all.x=TRUE)
cef_630_top$gene <- NULL
cef_630_top$name <- gsub('_', ' ', cef_630_top$name)
rownames(cef_630_top) <- cef_630_top$name
cef_630_top$name <- NULL
cef_mock_top <- merge(cef_mock_top, genes, by='gene', all.x=TRUE)
cef_mock_top$gene <- NULL
cef_mock_top$name <- gsub('_', ' ', cef_mock_top$name)
rownames(cef_mock_top) <- cef_mock_top$name
cef_mock_top$name <- NULL
clinda_630_top <- merge(clinda_630_top, genes, by='gene', all.x=TRUE)
clinda_630_top$gene <- NULL
clinda_630_top$name <- gsub('_', ' ', clinda_630_top$name)
rownames(clinda_630_top) <- clinda_630_top$name
clinda_630_top$name <- NULL
clinda_mock_top <- merge(clinda_mock_top, genes, by='gene', all.x=TRUE)
clinda_mock_top$gene <- NULL
clinda_mock_top$name <- gsub('_', ' ', clinda_mock_top$name)
rownames(clinda_mock_top) <- clinda_mock_top$name
clinda_mock_top$name <- NULL
strep_630_top <- merge(strep_630_top, genes, by='gene', all.x=TRUE)
strep_630_top$gene <- NULL
strep_630_top$name <- gsub('_', ' ', strep_630_top$name)
rownames(strep_630_top) <- strep_630_top$name
strep_630_top$name <- NULL
strep_mock_top <- merge(strep_mock_top, genes, by='gene', all.x=TRUE)
strep_mock_top$gene <- NULL
strep_mock_top$name <- gsub('_', ' ', strep_mock_top$name)
rownames(strep_mock_top) <- strep_mock_top$name
strep_mock_top$name <- NULL
rm(genes)

# Reorder by expression in treatment group
cef_630_top <- cef_630_top[order(cef_630_top$cef_630_reads),]
cef_mock_top <- cef_mock_top[order(cef_mock_top$cef_mock_reads),]
clinda_630_top <- clinda_630_top[order(clinda_630_top$clinda_630_reads),]
clinda_mock_top <- clinda_mock_top[order(clinda_mock_top$clinda_mock_reads),]
strep_630_top <- strep_630_top[order(strep_630_top$strep_630_reads),]
strep_mock_top <- strep_mock_top[order(strep_mock_top$strep_mock_reads),]

# Convert to matrices for barplots
cef_630_top <- as.matrix(t(cef_630_top))
cef_mock_top <- as.matrix(t(cef_mock_top))
clinda_630_top <- as.matrix(t(clinda_630_top))
clinda_mock_top <- as.matrix(t(clinda_mock_top))
strep_630_top <- as.matrix(t(strep_630_top))
strep_mock_top <- as.matrix(t(strep_mock_top))

# Reverse the row order so infected plots first
cef_630_top <- cef_630_top[rev(rownames(cef_630_top)),]
cef_mock_top <- cef_mock_top[rev(rownames(cef_mock_top)),]
clinda_630_top <- clinda_630_top[rev(rownames(clinda_630_top)),]
clinda_mock_top <- clinda_mock_top[rev(rownames(clinda_mock_top)),]
strep_630_top <- strep_630_top[rev(rownames(strep_630_top)),]
strep_mock_top <- strep_mock_top[rev(rownames(strep_mock_top)),]

#-------------------------------------------------------------------------------------------------------------------------#

# Format pathway info

# Renames columns
colnames(cef_pathways) <- c('cef_630','cef_mock')
colnames(clinda_pathways) <- c('clinda_630','clinda_mock')
colnames(strep_pathways) <- c('strep_630','strep_mock')

# Calculate differences
cef_pathways$difference <- abs(cef_pathways$cef_630 - cef_pathways$cef_mock)
clinda_pathways$difference <- abs(clinda_pathways$clinda_630 - clinda_pathways$clinda_mock)
strep_pathways$difference <- abs(strep_pathways$strep_630 - strep_pathways$strep_mock)

# Calculate confidence interval
diffs <- c(cef_pathways$difference, clinda_pathways$difference, strep_pathways$difference)
confidence <- as.numeric(quantile(diffs, probs=0.05))
rm(diffs)

# Find the largest changes, take top 5
cef_pathways <- cef_pathways[order(-cef_pathways$difference),][1:6,]
clinda_pathways <- clinda_pathways[order(-clinda_pathways$difference),][1:6,]
strep_pathways <- strep_pathways[order(-strep_pathways$difference),][1:6,]

# Remove Global pathways (always top since its the largest)
cef_pathways <- cef_pathways[-1, ]
clinda_pathways <- clinda_pathways[-1, ]
strep_pathways <- strep_pathways[-1, ]

# Recalculate differences for plotting
cef_pathways$difference <- cef_pathways$cef_630 - cef_pathways$cef_mock
cef_pathways$cef_630 <- NULL
cef_pathways$cef_mock <- NULL
clinda_pathways$difference <- clinda_pathways$clinda_630 - clinda_pathways$clinda_mock
clinda_pathways$clinda_630 <- NULL
clinda_pathways$clinda_mock <- NULL
strep_pathways$difference <- strep_pathways$strep_630 - strep_pathways$strep_mock
strep_pathways$strep_630 <- NULL
strep_pathways$strep_mock <- NULL

# Transform data and format for plotting
cef_pathways <- -log2(abs(cef_pathways))
clinda_pathways <- log2(clinda_pathways)
strep_pathways <- log2(strep_pathways)

# Format pathway names
pathway_names <- c(gsub('_', ' ', rownames(strep_pathways)), gsub('_', ' ', rownames(cef_pathways)), gsub('_', ' ', rownames(clinda_pathways)))
cef_pathways <- cef_pathways$difference
clinda_pathways <- clinda_pathways$difference
strep_pathways <- strep_pathways$difference

# Add blank rows for proper barplotting
cef_pathways <- c(rep(0,6), cef_pathways)
clinda_pathways <- c(rep(0,12), clinda_pathways)

#--------------------------------------------------------------------------------------------------#

# Generate figure
pdf(file=plot_file, width=8, height=12)
layout(matrix(c(1,2,3,
                4,4,4),
              nrow=2, ncol=3, byrow = TRUE))

#------------------#

# Streptomycin
# Mock-infected
par(mar=c(4,16,3,1), mgp=c(2.5, 0.75, 0), las=1, xaxs='i')
barplot(strep_mock_top, xaxt='n', xlim=c(0,14), beside=TRUE, horiz=TRUE, 
        xlab='', ylab='', col=c(noabx_col, strep_col), cex.names=1) 
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2))
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.5, cex=0.8)
mtext('Mock-infected', side=3, padj=-0.3, cex=0.9)
mtext('a', side=2, padj=-11, adj=18, cex=1.2, font=2)

# 630-infected
par(mar=c(4,16,3,1), mgp=c(2.5, 0.75, 0), las=1, xaxs='i')
barplot(strep_630_top, xaxt='n', xlim=c(0,14), beside=TRUE, horiz=TRUE,
        xlab='', ylab='', col=c(noabx_col, strep_col), cex.names=1) 
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2))
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.5, cex=0.8)
mtext(expression(paste(italic('C. difficile'),' 630-infected')), side=3, padj=-0.3, cex=0.9)
mtext('b', side=2, padj=-11, adj=17, cex=1.2, font=2)

# Legend
par(mar=c(0,0,0,1))
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-5,5), ylim=c(-10,10))
legend('center', legend=c('Streptomycin-pretreated','No Antibiotics (No CDI)'), pt.bg=c(strep_col,noabx_col), 
       pch=22, pt.cex=2.4, cex=1.3)

#------------------#

# Cefoperazone
# Mock-infected
par(mar=c(4,16,3,1), mgp=c(2.5, 0.75, 0), las=1, xaxs='i')
barplot(cef_mock_top, xaxt='n', xlim=c(0,14), beside=TRUE, horiz=TRUE, 
        xlab='', ylab='', col=c(noabx_col, cef_col), cex.names=1) 
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2))
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.5, cex=0.8)
mtext('Mock-infected', side=3, padj=-0.3, cex=0.9)
mtext('c', side=2, padj=-11, adj=18, cex=1.2, font=2)

# 630-infected
par(mar=c(4,16,3,1), mgp=c(2.5, 0.75, 0), las=1, xaxs='i')
barplot(cef_630_top, xaxt='n', xlim=c(0,14), beside=TRUE, horiz=TRUE, 
        xlab='', ylab='', col=c(noabx_col, cef_col), cex.names=1) 
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2))
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.5, cex=0.8)
mtext(expression(paste(italic('C. difficile'),' 630-infected')), side=3, padj=-0.3, cex=0.9)
mtext('d', side=2, padj=-11, adj=17, cex=1.2, font=2)

# Legend
par(mar=c(0,0,0,1))
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-5,5), ylim=c(-10,10))
legend('center', legend=c('Cefoperazone-pretreated','No Antibiotics (No CDI)'), pt.bg=c(cef_col,noabx_col), 
       pch=22, pt.cex=2.4, cex=1.3)

#------------------#

# Clindamycin
# Mock-infected
par(mar=c(4,16,3,1), mgp=c(2.5, 0.75, 0), las=1, xaxs='i')
barplot(clinda_mock_top, xaxt='n', xlim=c(0,14), beside=TRUE, horiz=TRUE, 
        xlab='', ylab='', col=c(noabx_col, clinda_col), cex.names=1) 
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2))
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.5, cex=0.8)
mtext('Mock-infected', side=3, padj=-0.3, cex=0.9)
mtext('d', side=2, padj=-11, adj=18, cex=1.2, font=2)

# 630-infected
par(mar=c(4,16,3,1), mgp=c(2.5, 0.75, 0), las=1, xaxs='i')
barplot(clinda_630_top, xaxt='n', xlim=c(0,14), beside=TRUE, horiz=TRUE, 
        xlab='', ylab='', col=c(noabx_col, clinda_col), cex.names=1) 
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2))
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.5, cex=0.8)
mtext(expression(paste(italic('C. difficile'),' 630-infected')), side=3, padj=-0.3, cex=0.9)
mtext('f', side=2, padj=-11, adj=30, cex=1.2, font=2)

# Legend
par(mar=c(0,0,0,1))
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-5,5), ylim=c(-10,10))
legend('center', legend=c('Clindamycin-pretreated','No Antibiotics (No CDI)'), pt.bg=c(clinda_col,noabx_col), 
       pch=22, pt.cex=2.4, cex=1.3)

#-------------------#

# Overrepresented pathways
par(mar=c(15,4,1,2), las=1, mgp=c(1.6,0.7,0))
plot(0, type='n', xlab='', xaxt='n', yaxt='n', ylab=as.expression(bquote(paste(Delta,' Fold Abundance (',log[2],')'))), xlim=c(0.5,20), ylim=c(-12,12))
abline(h=0, lwd=1.5)
abline(h=c(-confidence,confidence), lwd=1.2, lty=5, col='gray30')
axis(side=2, at=seq(-12,12,3), labels=c(12,9,6,3,0,3,6,9,12))
text(x=c(2.6,2.39), y=c(12,-12), cex=0.9,
     labels=c(as.expression(bquote(paste('Greater in ',italic('C. difficile'),'-infected Metatranscriptome'))), 'Greater in Mock-infected Metatranscriptome'))
legend('topright', legend=c('Streptomycin-pretreated','Cefoperazone-pretreated','Clindamycin-pretreated'),
       pt.bg=c(strep_col, cef_col, clinda_col), pch=22, pt.cex=1.7, col='black', bty='n')
text(cex=1, x=c(seq(1.2,6,1.2),seq(8.4,13.2,1.2),seq(15.6,20.4,1.2)), y=-14, pathway_names, xpd=TRUE, srt=60, pos=2)
mtext('d', side=2, line=2, las=2, adj=1.5, padj=-6, cex=1.2, font=2)

# Add groups
barplot(strep_pathways, xlim=c(0.5,20), ylim=c(-12,12), col=adjustcolor(strep_col, alpha.f=0.75), yaxt='n', add=TRUE) # Streptomycin
barplot(cef_pathways, xlim=c(0.5,20), ylim=c(-12,12), col=adjustcolor(cef_col, alpha.f=0.75), yaxt='n', add=TRUE) # Cefoperazone
barplot(clinda_pathways, xlim=c(0.5,20), ylim=c(-12,12), col=adjustcolor(clinda_col, alpha.f=0.75), yaxt='n', add=TRUE) # Clindamycin

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
