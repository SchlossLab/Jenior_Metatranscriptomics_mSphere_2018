
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

# Calculate differences between infected to mock
cef_annotated$diff <- abs(cef_annotated$cef_630_reads - cef_annotated$cef_mock_reads)
clinda_annotated$diff <- abs(clinda_annotated$clinda_630_reads - clinda_annotated$clinda_mock_reads)
strep_annotated$diff <- abs(strep_annotated$strep_630_reads - strep_annotated$strep_mock_reads)
cef_pathways$diff <- abs(cef_pathways$cef_630_reads - cef_pathways$cef_mock_reads)
clinda_pathways$diff <- abs(clinda_pathways$clinda_630_reads - clinda_pathways$clinda_mock_reads)
strep_pathways$diff <- abs(strep_pathways$strep_630_reads - strep_pathways$strep_mock_reads)

# Rank differences and subset to top differences
cef_annotated <- cef_annotated[order(-cef_annotated$diff),][1:20,]
clinda_annotated <- clinda_annotated[order(-clinda_annotated$diff),][1:20,]
strep_annotated <- strep_annotated[order(-strep_annotated$diff),][1:20,]
cef_pathways <- cef_pathways[order(-cef_pathways$diff),][1:10,]
clinda_pathways <- clinda_pathways[order(-clinda_pathways$diff),][1:10,]
strep_pathways <- strep_pathways[order(-strep_pathways$diff),][1:10,]

# Log transform expression
cef_annotated[,c(2:4)] <- log2(cef_annotated[,c(2:4)])
clinda_annotated[,c(2:4)] <- log2(clinda_annotated[,c(2:4)])
strep_annotated[,c(2:4)] <- log2(strep_annotated[,c(2:4)])
cef_pathways[,c(2:4)] <- log2(cef_pathways[,c(2:4)])
clinda_pathways[,c(2:4)] <- log2(clinda_pathways[,c(2:4)])
strep_pathways[,c(2:4)] <- log2(strep_pathways[,c(2:4)])

# Format names for plotting
cef_annotated$gene <- gsub('_', ' ', cef_annotated$gene)
cef_annotated$gene <- gsub('mostly Fe transport', 'Iron transport', cef_annotated$gene)
cef_annotated$gene <- gsub('phosphopyruvate hydratase', 'Phosphopyruvate hydratase', cef_annotated$gene)
cef_annotated$gene <- gsub('outer membrane receptor protein', 'Outer membrane receptor protein', cef_annotated$gene)
cef_annotated$gene <- gsub('Predicted dehydrogenases and related proteins', 'Predicted dehydrogenases', cef_annotated$gene)
cef_annotated$gene <- gsub('SusC/RagA family TonB-linked outer membrane protein', 'SusC/RagA family TonB-linked protein', cef_annotated$gene)
cef_annotated$gene <- gsub('SusD family.', 'SusD family protein', cef_annotated$gene)
cef_annotated$gene <- gsub('methylmalonyl-CoA mutase \\(EC\\:5\\.4\\.99\\.2\\)', 'Methylmalonyl-CoA mutase', cef_annotated$gene)
cef_annotated$gene <- gsub('two-component system sensor histidine kinase/response regulator', 'Two-component histidine sensor regulator', cef_annotated$gene)
cef_annotated <- cef_annotated[cef_annotated$gene != '50S ribosomal protein L5',]
cef_annotated <- cef_annotated[cef_annotated$gene != 'putative',]
cef_annotated <- cef_annotated[c(1:15),]
rownames(cef_annotated) <- cef_annotated$gene
cef_annotated$gene <- NULL

clinda_annotated$gene <- gsub('_', ' ', clinda_annotated$gene)
clinda_annotated$gene <- gsub('fructose-bisphosphate aldolase', 'Fructose-bisphosphate aldolase', clinda_annotated$gene)
clinda_annotated$gene <- gsub('transcriptional regulator Spx', 'Transcriptional regulator Spx', clinda_annotated$gene)
clinda_annotated$gene <- gsub('glyceraldehyde 3-phosphate dehydrogenase \\(EC\\:1\\.2\\.1\\.12\\)', 'Glyceraldehyde 3-P dehydrogenase', clinda_annotated$gene)
clinda_annotated$gene <- gsub('Multi Drug ABC transporter transmembrane subunit family protein', 'Multi-Drug ABC transporter subunit', clinda_annotated$gene)
clinda_annotated$gene <- gsub('cell division initiation protein', 'Cell division initiation protein', clinda_annotated$gene)
clinda_annotated$gene <- gsub('cell division protein GpsB', 'Cell division protein GpsB', clinda_annotated$gene)
clinda_annotated$gene <- gsub('pyruvate formate-lyase activating enzyme \\(EC\\:1\\.97\\.1\\.4)', 'Pyruvate formate-lyase activating enzyme', clinda_annotated$gene)
clinda_annotated$gene <- gsub('phosphoenolpyruvate-protein phosphotransferase', 'Phosphoenolpyruvate phosphotransferase', clinda_annotated$gene)
clinda_annotated$gene <- gsub('phosphoglycerate mutase', 'Phosphoglycerate mutase', clinda_annotated$gene)
clinda_annotated$gene <- gsub('phosphopyruvate hydratase', 'Phosphopyruvate hydratase', clinda_annotated$gene)
clinda_annotated$gene <- gsub('ribosome-associated factor Y', 'Ribosome-associated factor Y', clinda_annotated$gene)
clinda_annotated <- clinda_annotated[clinda_annotated$gene != 'uncharacterized LOC100521496',]
clinda_annotated <- clinda_annotated[clinda_annotated$gene != '30S ribosomal protein S16',]
clinda_annotated <- clinda_annotated[c(1:15),]
rownames(clinda_annotated) <- clinda_annotated$gene
clinda_annotated$gene <- NULL

strep_annotated$gene <- gsub('_', ' ', strep_annotated$gene)
strep_annotated$gene <- gsub('glyceraldehyde 3-phosphate dehydrogenase \\(EC\\:1\\.2\\.1\\.12)', 'Glyceraldehyde 3-P dehydrogenase', strep_annotated$gene)
strep_annotated$gene <- gsub('alcohol dehydrogenase', 'Alcohol dehydrogenase', strep_annotated$gene)
strep_annotated$gene <- gsub('universal stress protein', 'Universal stress protein', strep_annotated$gene)
strep_annotated$gene <- gsub('fructose-bisphosphate aldolase', 'Fructose-bisphosphate aldolase', strep_annotated$gene)
strep_annotated$gene <- gsub('elongation factor Tu \\(EC\\:3\\.6\\.5\\.3)', 'Elongation factor Tu', strep_annotated$gene)
strep_annotated$gene <- gsub('alkyl hydroperoxide reductase', 'Alkyl hydroperoxide reductase', strep_annotated$gene)
strep_annotated$gene <- gsub('acetyltransferase', 'Acetyltransferase', strep_annotated$gene)
strep_annotated$gene <- gsub('phosphoenolpyruvate-protein phosphotransferase', 'Phosphoenolpyruvate phosphotransferase', strep_annotated$gene)
strep_annotated$gene <- gsub('phosphoglycerate kinase \\(EC\\:2\\.7\\.2\\.3)', 'Phosphoglycerate kinase', strep_annotated$gene)
strep_annotated$gene <- gsub('phosphoglycerate mutase', 'Phosphoglycerate mutase', strep_annotated$gene)
strep_annotated$gene <- gsub('transcriptional regulator', 'Transcriptional regulator', strep_annotated$gene)
strep_annotated <- strep_annotated[strep_annotated$gene != '50S ribosomal protein L12',]
strep_annotated <- strep_annotated[strep_annotated$gene != '50S ribosomal protein L10',]
strep_annotated <- strep_annotated[strep_annotated$gene != 'ribosomal protein',]
strep_annotated <- strep_annotated[c(1:15),]
rownames(strep_annotated) <- strep_annotated$gene
strep_annotated$gene <- NULL

cef_pathways$pathway <- gsub('_', ' ', cef_pathways$pathway)
cef_pathways <- cef_pathways[cef_pathways$pathway != 'Ribosome:Translation',]
cef_pathways <- cef_pathways[cef_pathways$pathway != 'Metabolic pathways',]
cef_pathways$pathway <- gsub(':', ': ', cef_pathways$pathway)
cef_pathways <- cef_pathways[c(1:5),]
rownames(cef_pathways) <- cef_pathways$pathway
cef_pathways$pathway <- NULL

clinda_pathways$pathway <- gsub('_', ' ', clinda_pathways$pathway)
clinda_pathways <- clinda_pathways[clinda_pathways$pathway != 'Ribosome:Translation',]
clinda_pathways <- clinda_pathways[clinda_pathways$pathway != 'Metabolic pathways',]
clinda_pathways$pathway <- gsub(':Membrane transport', '', clinda_pathways$pathway)
clinda_pathways$pathway <- gsub(' \\(PTS\\)', '', clinda_pathways$pathway)
clinda_pathways$pathway <- gsub(':Energy metabolism:Energy metabolism', '', clinda_pathways$pathway)
clinda_pathways <- clinda_pathways[c(1:5),]
rownames(clinda_pathways) <- clinda_pathways$pathway
clinda_pathways$pathway <- NULL

strep_pathways$pathway <- gsub('_', ' ', strep_pathways$pathway)
strep_pathways <- strep_pathways[strep_pathways$pathway != 'Ribosome:Translation',]
strep_pathways <- strep_pathways[strep_pathways$pathway != 'Metabolic pathways',]
strep_pathways$pathway <- gsub(':Translation', '', strep_pathways$pathway)
strep_pathways$pathway <- gsub(':Energy metabolism', '', strep_pathways$pathway)
strep_pathways$pathway <- gsub(':Carbohydrate metabolism', '', strep_pathways$pathway)
strep_pathways$pathway <- gsub(':Membrane transport', '', strep_pathways$pathway)
strep_pathways <- strep_pathways[c(1:5),]
rownames(strep_pathways) <- strep_pathways$pathway
strep_pathways$pathway <- NULL

# Reorder to plot correctly
cef_annotated <- cef_annotated[order(cef_annotated$diff),]
clinda_annotated <- clinda_annotated[order(clinda_annotated$diff),]
strep_annotated <- strep_annotated[order(strep_annotated$diff),]
cef_pathways <- cef_pathways[order(cef_pathways$diff),]
clinda_pathways <- clinda_pathways[order(clinda_pathways$diff),]
strep_pathways <- strep_pathways[order(strep_pathways$diff),]

# Convert to matrices for barplots
cef_annotated <- as.matrix(t(cef_annotated))
cef_pathways <- as.matrix(t(cef_pathways))
clinda_annotated <- as.matrix(t(clinda_annotated))
clinda_pathways <- as.matrix(t(clinda_pathways))
strep_annotated <- as.matrix(t(strep_annotated))
strep_pathways <- as.matrix(t(strep_pathways))

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
