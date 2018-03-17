
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

# Subset to pathway annotations for separate analysis
cef_pathways <- subset(cef_annotated, !is.na(cef_annotated$pathways))
cef_pathways$gene <- NULL
cef_pathways$organism <- NULL
cef_pathways$description <- NULL
cef_pathways$ko <- NULL
clinda_pathways <- subset(clinda_annotated, !is.na(clinda_annotated$pathways))
clinda_pathways$gene <- NULL
clinda_pathways$organism <- NULL
clinda_pathways$description <- NULL
clinda_pathways$ko <- NULL
strep_pathways <- subset(strep_annotated, !is.na(strep_annotated$pathways))
strep_pathways$gene <- NULL
strep_pathways$organism <- NULL
strep_pathways$description <- NULL
strep_pathways$ko <- NULL

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

# Aggregate identical genes, regardless of organism
cef_annotated <- aggregate(cbind(cef_annotated$cef_630_metaT_reads, cef_annotated$cef_mock_metaT_reads), by=list(cef_annotated$description), FUN=sum)
colnames(cef_annotated) <- c('gene', 'cef_630_reads', 'cef_mock_reads')
clinda_annotated <- aggregate(cbind(clinda_annotated$clinda_630_metaT_reads, clinda_annotated$clinda_mock_metaT_reads), by=list(clinda_annotated$description), FUN=sum)
colnames(clinda_annotated) <- c('gene', 'clinda_630_reads', 'clinda_mock_reads')
strep_annotated <- aggregate(cbind(strep_annotated$strep_630_metaT_reads, strep_annotated$strep_mock_metaT_reads), by=list(strep_annotated$description), FUN=sum)
colnames(strep_annotated) <- c('gene', 'strep_630_reads', 'strep_mock_reads')

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

# Calculate and take top differences
cef_annotated$diff <- abs(cef_annotated$cef_630_reads - cef_annotated$cef_mock_reads)
cef_annotated <- cef_annotated[order(-cef_annotated$diff),]
top_cef_genes <- cef_annotated$gene[1:25]
cef_annotated$diff <- NULL
clinda_annotated$diff <- abs(clinda_annotated$clinda_630_reads - clinda_annotated$clinda_mock_reads)
clinda_annotated <- clinda_annotated[order(-clinda_annotated$diff),]
top_clinda_genes <- clinda_annotated$gene[1:25]
clinda_annotated$diff <- NULL
strep_annotated$diff <- abs(strep_annotated$strep_630_reads - strep_annotated$strep_mock_reads)
strep_annotated <- strep_annotated[order(-strep_annotated$diff),]
top_strep_genes <- strep_annotated$gene[1:25]
strep_annotated$diff <- NULL
top_genes <- unique(c(top_cef_genes, top_clinda_genes, top_strep_genes))
rm(top_cef_genes, top_clinda_genes, top_strep_genes)

# Merge treatment groups
final_annotated <- merge(strep_annotated, cef_annotated, by='gene', all=TRUE)
final_annotated <- merge(final_annotated, clinda_annotated, by='gene', all=TRUE)
final_annotated[is.na(final_annotated)] <- 0 # Undetectable genes set to 0
rownames(final_annotated) <- final_annotated$gene
final_annotated$gene <- NULL
rm(cef_annotated, clinda_annotated, strep_annotated)

# Subset to top genes
final_annotated <- subset(final_annotated, rownames(final_annotated) %in% top_genes)
rm(top_genes)

# Reorder to the largest variance with largest variance
final_annotated$vars <- RowVar(final_annotated)
final_annotated <- final_annotated[order(-final_annotated$vars),]
final_annotated <- final_annotated[c(1:25),]
final_annotated <- final_annotated[order(final_annotated$vars),]
final_annotated$vars <- NULL

# Final transformation
final_annotated <- log2(final_annotated + 1)

# Reformat names
final_annotated$gene <- rownames(final_annotated)
final_annotated$gene <- gsub('_', ' ', final_annotated$gene)
final_annotated$gene <- gsub(' \\(EC\\:1\\.2\\.1\\.12\\)', '', final_annotated$gene)
final_annotated$gene <- gsub(' \\(EC\\:3\\.6\\.5\\.3\\)', '', final_annotated$gene)
final_annotated$gene <- gsub(' \\(EC\\:2\\.7\\.2\\.3\\)', '', final_annotated$gene)
final_annotated$gene <- gsub(' \\(EC\\:1\\.97\\.1\\.4\\)', '', final_annotated$gene)
final_annotated$gene <- gsub('phosphoenolpyruvate', 'PEP', final_annotated$gene)
final_annotated$gene <- gsub('Multi Drug ABC transporter transmembrane subunit family protein', 'Multi Drug ABC transporter', final_annotated$gene)
final_annotated$gene <- gsub('glyceraldehyde 3-phosphate dehydrogenase', 'GAPDH', final_annotated$gene)
final_annotated$gene <- capitalize(final_annotated$gene)
rownames(final_annotated) <- final_annotated$gene
final_annotated$gene <- NULL

#------------------#

# Pathway analysis
cef_pathways <- aggregate(cef_pathways[,1:2], by=list(cef_pathways$pathways), FUN=sum)
rownames(cef_pathways) <- cef_pathways$Group.1
cef_pathways$Group.1 <- NULL
clinda_pathways <- aggregate(clinda_pathways[,1:2], by=list(clinda_pathways$pathways), FUN=sum)
rownames(clinda_pathways) <- clinda_pathways$Group.1
clinda_pathways$Group.1 <- NULL
strep_pathways <- aggregate(strep_pathways[,1:2], by=list(strep_pathways$pathways), FUN=sum)
rownames(strep_pathways) <- strep_pathways$Group.1
strep_pathways$Group.1 <- NULL

# Remove uninformative pathways
cef_pathways <- subset(cef_pathways, !rownames(cef_pathways) %in% c('Ribosome:Translation','Metabolic_pathways','Carbon_metabolism'))
clinda_pathways <- subset(clinda_pathways, !rownames(clinda_pathways) %in% c('Ribosome:Translation','Metabolic_pathways','Carbon_metabolism'))
strep_pathways <- subset(strep_pathways, !rownames(strep_pathways) %in% c('Ribosome:Translation','Metabolic_pathways','Carbon_metabolism'))

# Calculate and take top differences
cef_pathways$diff <- abs(cef_pathways$cef_630_metaT_reads - cef_pathways$cef_mock_metaT_reads)
cef_pathways <- cef_pathways[order(-cef_pathways$diff),]
top_cef_pathways <- rownames(cef_pathways[1:15,])
cef_pathways$diff <- NULL
clinda_pathways$diff <- abs(clinda_pathways$clinda_630_metaT_reads - clinda_pathways$clinda_mock_metaT_reads)
clinda_pathways <- clinda_pathways[order(-clinda_pathways$diff),]
top_clinda_pathways <- rownames(clinda_pathways[1:15,])
clinda_pathways$diff <- NULL
strep_pathways$diff <- abs(strep_pathways$strep_630_metaT_reads - strep_pathways$strep_mock_metaT_reads)
strep_pathways <- strep_pathways[order(-strep_pathways$diff),]
top_strep_pathways <- rownames(strep_pathways[1:15,])
strep_pathways$diff <- NULL
top_pathways <- unique(c(top_cef_pathways, top_clinda_pathways, top_strep_pathways))
rm(top_cef_pathways, top_clinda_pathways, top_strep_pathways)

# Merge treatment groups
final_pathways <- merge(strep_pathways, cef_pathways, by='row.names', all=TRUE)
rownames(final_pathways) <- final_pathways$Row.names
final_pathways$Row.names <- NULL
final_pathways <- merge(final_pathways, clinda_pathways, by='row.names', all=TRUE)
rownames(final_pathways) <- final_pathways$Row.names
final_pathways$Row.names <- NULL
final_pathways[is.na(final_pathways)] <- 0 # Undetectable genes set to 0
rm(cef_pathways, clinda_pathways, strep_pathways)

# Subset to top 
final_pathways <- subset(final_pathways, rownames(final_pathways) %in% top_pathways)
rm(top_pathways)

# Reorder to the largest variance with largest variance
final_pathways$vars <- RowVar(final_pathways)
final_pathways <- final_pathways[order(-final_pathways$vars),]
final_pathways <- final_pathways[c(1:10),]
final_pathways <- final_pathways[order(final_pathways$vars),]
final_pathways$vars <- NULL

# Final transformation
final_pathways <- log2(final_pathways + 1)

# Reformat names
path_names <- as.vector(rownames(final_pathways))
for (x in c(1:length(path_names))) {
  path_names[x] <- strsplit(path_names[x], split=':')[[1]][1]
}
path_names <- gsub('_', ' ', path_names)
rownames(final_pathways) <- path_names
rm(path_names)

#--------------------------------------------------------------------------------------------------#

# Generate figure
pdf(file=plot_file, width=7, height=11)

#------------------#

# Pathways

par(mar=c(15,4,1,1), mgp=c(2.5, 0.75, 0), las=1, xaxs='i')

barplot(as.matrix(t(final_pathways)), xaxt='n', ylim=c(0,15), beside=TRUE,  
        xlab='', ylab='', angle=20, density=c(500,20), 
        col=c(strep_col,strep_col,cef_col,cef_col,clinda_col,clinda_col), cex.names=0.9)
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2))
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=2, padj=2.5, cex=0.9)



pdf(file=plot_file, width=7, height=12)

# Genes
par(mar=c(3,15.5,1,1), mgp=c(2.5, 0.75, 0), las=1, xaxs='i', yaxs='i')
barplot(as.matrix(t(final_annotated)), xaxt='n', xlim=c(0,14), ylim=c(0,176), beside=TRUE, horiz=TRUE, 
        xlab='', ylab='', angle=20, density=c(500,20), 
        col=c(strep_col,strep_col,cef_col,cef_col,clinda_col,clinda_col), cex.names=0.9)
abline(h=seq(7.5,168.5,7), lty=5, col='black')
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2))
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.5)

dev.off()



















text(x=12.5, y=5.5, 'Infection', cex=1.1)
legend('bottomright', legend=c(expression(italic('C. difficile')),'Mock'), pt.bg=c('black','white'),
       pch=22, pt.cex=1.5, cex=0.9)







# Final name changes to add pathways
colnames(strep_annotated) <- c('(A) ATP synthase subunit B',
                               '(A) F0F1 ATP synthase subunit A',
                               'Cell division initiation protein',
                               '(B,D) L-lactate dehydrogenase',
                               '(C) Transcriptional regulator',
                               '(G) Phosphoenolpyruvate phosphotransferase',
                               '(B) Phosphoglycerate kinase',
                               '(I) Alkyl hydroperoxide reductase',
                               '(D) Acetyltransferase',
                               '(B,F,K) Phosphoglycerate mutase')
colnames(cef_annotated) <- c('(D,J) 2-isopropylmalate synthase',
                             '(C,E) Superfamily II DNA and RNA helicases',
                             '(C) Translation elongation factor P',
                             'Na+/glucose cotransporter',
                             'Outer membrane receptor protein',
                             '(J) Methylmalonyl-CoA mutase',
                             'DnaK suppressor protein',
                             'Two-component histidine sensor regulator',
                             '(H) SusD family protein',
                             '(B,D) Pyruvate kinase')
colnames(clinda_annotated) <- c('(A) ATP synthase subunit B',
                                'Thioredoxin',
                                '(D) Acetyltransferase',
                                '(G) Phosphoenolpyruvate phosphotransferase',
                                '(D) Formate Acetyltransferase',
                                '(B,F,K) Phosphoglycerate mutase',
                                'Cell division protein GpsB',
                                '(D) Pyruvate formate-lyase activating enzyme',
                                '(B,F) Phosphopyruvate hydratase',
                                '(C) Ribosome-associated factor Y')


#--------------------------------------------------------------------------------------------------#

# Generate figure
pdf(file=plot_file, width=5, height=11)
layout(matrix(c(1,1,
                1,1,
                2,2,
                2,2,
                3,3,
                3,3,
                4,4),
              nrow=7, ncol=2, byrow = TRUE))

#------------------#


# Streptomycin
par(mar=c(3,15,2,1), mgp=c(2.5, 0.75, 0), las=1, xaxs='i')
barplot(strep_annotated, xaxt='n', xlim=c(0,14), beside=TRUE, horiz=TRUE,
        xlab='', ylab='', col=c('white','black'), cex.names=0.8)
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2))
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.2, cex=0.75)
title('Streptomycin-pretreated', line=0.5, cex.main=1.2, col.main=strep_col, font.main=2)
mtext('A', side=2, padj=-10, adj=16, font=2)
text(x=12.5, y=5.5, 'Infection', cex=1.1)
legend('bottomright', legend=c(expression(italic('C. difficile')),'Mock'), pt.bg=c('black','white'),
       pch=22, pt.cex=1.5, cex=0.9)

#------------------#

# Cefoperazone
par(mar=c(3,15,2,1), mgp=c(2.5, 0.75, 0), las=1, xaxs='i')
barplot(cef_annotated, xaxt='n', xlim=c(0,14), beside=TRUE, horiz=TRUE,
        xlab='', ylab='', col=c('white','black'), cex.names=0.8)
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2))
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.2, cex=0.75)
title('Cefoperazone-pretreated', line=0.5, cex.main=1.2, col.main=cef_col, font.main=2)
mtext('B', side=2, padj=-9, adj=16, font=2)
text(x=12.5, y=5.5, 'Infection', cex=1.1)
legend('bottomright', legend=c(expression(italic('C. difficile')),'Mock'), pt.bg=c('black','white'),
       pch=22, pt.cex=1.5, cex=0.9)

#------------------#

# Clindamycin
par(mar=c(3,15,2,1), mgp=c(2.5, 0.75, 0), las=1, xaxs='i')
barplot(clinda_annotated, xaxt='n', xlim=c(0,14), beside=TRUE, horiz=TRUE,
        xlab='', ylab='', col=c('white','black'), cex.names=0.8)
box()
axis(1, at=seq(0,14,2), label=seq(0,14,2))
minor.ticks.axis(1, 10, mn=0, mx=14)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.2, cex=0.75)
title('Clindamycin-pretreated', line=0.5, cex.main=1.2, col.main=clinda_col, font.main=2)
mtext('C', side=2, padj=-9, adj=16, font=2)
text(x=12.5, y=5.5, 'Infection', cex=1.1)
legend('bottomright', legend=c(expression(italic('C. difficile')),'Mock'), pt.bg=c('black','white'),
       pch=22, pt.cex=1.5, cex=0.9)

#------------------#

# Pathway legend
par(mar=c(5,0,1,0))
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-12,12), ylim=c(-2,2))
legend('center', ncol=3, cex=0.75, pt.cex=0.9,
       pch=c('A','B','C','D','E','F','G','H','I','J','K'),
       c('Oxidative phosphorylation','Glycolysis / Gluconeogenesis','RNA folding, sorting & degradation',
         'Pyruvate metabolism','Homologous recombination',
         'Methane metabolism','Phosphotransferase system','Starch & sucrose metabolism',
         'Glutathione metabolism','Valine, leucine, and isoleucine biosyn.','Glycine, serine, and threonine metab.'))

text(x=-6.5, y=1.8, expression(paste('Most frequent pathways among genes in ', bold('A-C'), ':')), cex=0.9, xpd=TRUE)

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


