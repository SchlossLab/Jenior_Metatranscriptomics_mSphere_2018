
# Set up environment

# Load dependencies
deps <- c('vegan', 'biglm')
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Set seed for RNG
set.seed(6189)

# Load in functions
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define files

# Output plot
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/figures/figure_3.pdf'

# KEGG organism IDs
kegg_org_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/kegg_organisms.tsv'

# Metagenomes
cef_metagenome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Cefoperazone.DNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_metagenome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Clindamycin.DNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_metagenome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Streptomycin.DNA_reads2pangenome.all.norm.remove.annotated.txt'
conv_metagenome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Conventional.DNA_reads2pangenome.all.norm.remove.annotated.txt'

# Metatranscriptomes
cef_630_metatranscriptome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/cefoperazone_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'
cef_mock_metatranscriptome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/cefoperazone_mock.RNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_630_metatranscriptome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/clindamycin_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_mock_metatranscriptome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/clindamycin_mock.RNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_630_metatranscriptome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/streptomycin_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_mock_metatranscriptome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/streptomycin_mock.RNA_reads2pangenome.all.norm.remove.annotated.txt'
conv_metatranscriptome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/conventional.RNA_reads2pangenome.all.norm.remove.annotated.txt'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data

# Metagenomes
cef_metagenome <- read.delim(cef_metagenome_file, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(cef_metagenome) <- c('cef_metaG_reads', 'ko', 'gene', 'pathway')
clinda_metagenome <- read.delim(clinda_metagenome_file, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(clinda_metagenome) <- c('clinda_metaG_reads', 'ko', 'gene', 'pathway')
strep_metagenome <- read.delim(strep_metagenome_file, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(strep_metagenome) <- c('strep_metaG_reads', 'ko', 'gene', 'pathway')
conv_metagenome <- read.delim(conv_metagenome_file, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(conv_metagenome) <- c('conv_metaG_reads', 'ko', 'gene', 'pathway')
rm(cef_metagenome_file, clinda_metagenome_file, strep_metagenome_file, conv_metagenome_file, na.strings=c('','NA'))

# Metatranscriptomes
cef_630_metatranscriptome <- read.delim(cef_630_metatranscriptome_file, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(cef_630_metatranscriptome) <- c('cef_630_metaT_reads', 'ko', 'gene', 'pathway')
cef_mock_metatranscriptome <- read.delim(cef_mock_metatranscriptome_file, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(cef_mock_metatranscriptome) <- c('cef_mock_metaT_reads', 'ko', 'gene', 'pathway')
rm(cef_630_metatranscriptome_file, cef_mock_metatranscriptome_file)
clinda_630_metatranscriptome <- read.delim(clinda_630_metatranscriptome_file, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(clinda_630_metatranscriptome) <- c('clinda_630_metaT_reads', 'ko', 'gene', 'pathway')
clinda_mock_metatranscriptome <- read.delim(clinda_mock_metatranscriptome_file, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(clinda_mock_metatranscriptome) <- c('clinda_mock_metaT_reads', 'ko', 'gene', 'pathway')
rm(clinda_630_metatranscriptome_file, clinda_mock_metatranscriptome_file)
strep_630_metatranscriptome <- read.delim(strep_630_metatranscriptome_file, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(strep_630_metatranscriptome) <- c('strep_630_metaT_reads', 'ko', 'gene', 'pathway')
strep_mock_metatranscriptome <- read.delim(strep_mock_metatranscriptome_file, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(strep_mock_metatranscriptome) <- c('strep_mock_metaT_reads', 'ko', 'gene', 'pathway')
rm(strep_630_metatranscriptome_file, strep_mock_metatranscriptome_file)
conv_metatranscriptome <- read.delim(conv_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(conv_metatranscriptome) <- c('conv_metaT_reads', 'ko', 'gene', 'pathway')
rm(conv_final_reads_file)

# KEGG organism file
kegg_org <- read.delim(kegg_org_file, sep='\t', header=TRUE)
rm(kegg_org_file)

#-------------------------------------------------------------------------------------------------------------------------#

# Format metagenomic data for merging
cef_metagenome$ko <- NULL
cef_metagenome$gene <- NULL
cef_metagenome$pathway <- NULL
clinda_metagenome$ko <- NULL
clinda_metagenome$gene <- NULL
clinda_metagenome$pathway <- NULL
strep_metagenome$ko <- NULL
strep_metagenome$gene <- NULL
strep_metagenome$pathway <- NULL
conv_metagenome$ko <- NULL
conv_metagenome$gene <- NULL
conv_metagenome$pathway <- NULL

# Format metatranscriptomic data for merging
cef_630_metatranscriptome$ko <- NULL
cef_630_metatranscriptome$gene <- NULL
cef_630_metatranscriptome$pathway <- NULL
clinda_630_metatranscriptome$ko <- NULL
clinda_630_metatranscriptome$gene <- NULL
clinda_630_metatranscriptome$pathway <- NULL
strep_630_metatranscriptome$ko <- NULL
strep_630_metatranscriptome$gene <- NULL
strep_630_metatranscriptome$pathway <- NULL

# Merge metagenomic and metatranscriptomic data
cef_raw_reads <- clean_merge(cef_metagenome, cef_630_metatranscriptome)
cef_raw_reads <- clean_merge(cef_raw_reads, cef_mock_metatranscriptome)
clinda_raw_reads <- clean_merge(clinda_metagenome, clinda_630_metatranscriptome)
clinda_raw_reads <- clean_merge(clinda_raw_reads, clinda_mock_metatranscriptome)
strep_raw_reads <- clean_merge(strep_metagenome, strep_630_metatranscriptome)
strep_raw_reads <- clean_merge(strep_raw_reads, strep_mock_metatranscriptome)
conv_raw_reads <- clean_merge(conv_metagenome, conv_metatranscriptome)

rm(cef_metagenome, clinda_metagenome, strep_metagenome, conv_metagenome,
   cef_630_metatranscriptome, cef_mock_metatranscriptome, clinda_630_metatranscriptome, 
   clinda_mock_metatranscriptome, strep_630_metatranscriptome, strep_mock_metatranscriptome, 
   conv_metatranscriptome)

#-------------------------------------------------------------------------------------------------------------------------#

# Remove residual C. difficile mappings (str. 630 & 196)
cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdf:CD630', rownames(cef_raw_reads)),]), ]
cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdc:CD196', rownames(cef_raw_reads)),]), ]
clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdf:CD630', rownames(clinda_raw_reads)),]), ]
clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdc:CD196', rownames(clinda_raw_reads)),]), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdf:CD630', rownames(strep_raw_reads)),]), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdc:CD196', rownames(strep_raw_reads)),]), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdf:CD630', rownames(conv_raw_reads)),]), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdc:CD196', rownames(conv_raw_reads)),]), ]

# Remove genes with no metagenomic coverage
cef_raw_reads <- subset(cef_raw_reads, cef_metaG_reads != 0)
clinda_raw_reads <- subset(clinda_raw_reads, clinda_metaG_reads != 0)
strep_raw_reads <- subset(strep_raw_reads, strep_metaG_reads != 0)
conv_raw_reads <- subset(conv_raw_reads, conv_metaG_reads != 0)


round(median(c(nrow(subset(cef_raw_reads, cef_metaG_reads == 0 & (cef_mock_metaT_reads + cef_630_metaT_reads ) != 0)),
nrow(subset(clinda_raw_reads, clinda_metaG_reads == 0 & (clinda_mock_metaT_reads  + clinda_630_metaT_reads ) != 0)),
nrow(subset(strep_raw_reads, strep_metaG_reads == 0 & (strep_mock_metaT_reads  + strep_630_metaT_reads ) != 0)),
nrow(subset(conv_raw_reads, conv_metaG_reads == 0 & conv_metaT_reads  != 0)))))


# Rarefy read abundances
size <- round(min(colSums(cef_raw_reads[,c(1:3)]))*0.9)
cef_raw_reads$cef_metaG_reads <- t(rrarefy(cef_raw_reads$cef_metaG_reads, sample=size)) + 1
cef_raw_reads$cef_630_metaT_reads <- t(rrarefy(cef_raw_reads$cef_630_metaT_reads, sample=size))
cef_raw_reads$cef_mock_metaT_reads <- t(rrarefy(cef_raw_reads$cef_mock_metaT_reads, sample=size))
size <- round(min(colSums(clinda_raw_reads[,c(1:3)]))*0.9)
clinda_raw_reads$clinda_metaG_reads <- t(rrarefy(clinda_raw_reads$clinda_metaG_reads, sample=size)) + 1
clinda_raw_reads$clinda_630_metaT_reads <- t(rrarefy(clinda_raw_reads$clinda_630_metaT_reads, sample=size))
clinda_raw_reads$clinda_mock_metaT_reads <- t(rrarefy(clinda_raw_reads$clinda_mock_metaT_reads, sample=size))
size <- round(min(colSums(strep_raw_reads[,c(1:3)]))*0.9)
strep_raw_reads$strep_metaG_reads <- t(rrarefy(strep_raw_reads$strep_metaG_reads, sample=size)) + 1
strep_raw_reads$strep_630_metaT_reads <- t(rrarefy(strep_raw_reads$strep_630_metaT_reads, sample=size))
strep_raw_reads$strep_mock_metaT_reads <- t(rrarefy(strep_raw_reads$strep_mock_metaT_reads, sample=size))
size <- round(min(colSums(conv_raw_reads[,c(1:2)]))*0.9)
conv_raw_reads$conv_metaG_reads <- t(rrarefy(conv_raw_reads$conv_metaG_reads, sample=size)) + 1
conv_raw_reads$conv_metaT_reads <- t(rrarefy(conv_raw_reads$conv_metaT_reads, sample=size))
rm(size)

# Normalize metatranscriptomes to metagenomic coverage
cef_raw_reads$cef_630_metaT_reads <- cef_raw_reads$cef_630_metaT_reads / cef_raw_reads$cef_metaG_reads
cef_raw_reads$cef_mock_metaT_reads <- cef_raw_reads$cef_mock_metaT_reads / cef_raw_reads$cef_metaG_reads
cef_raw_reads$cef_metaG_reads <- NULL
clinda_raw_reads$clinda_630_metaT_reads <- clinda_raw_reads$clinda_630_metaT_reads / clinda_raw_reads$clinda_metaG_reads
clinda_raw_reads$clinda_mock_metaT_reads <- clinda_raw_reads$clinda_mock_metaT_reads / clinda_raw_reads$clinda_metaG_reads
clinda_raw_reads$clinda_metaG_reads <- NULL
strep_raw_reads$strep_630_metaT_reads <- strep_raw_reads$strep_630_metaT_reads / strep_raw_reads$strep_metaG_reads
strep_raw_reads$strep_mock_metaT_reads <- strep_raw_reads$strep_mock_metaT_reads / strep_raw_reads$strep_metaG_reads
strep_raw_reads$strep_metaG_reads <- NULL
conv_raw_reads$conv_metaT_reads <- conv_raw_reads$conv_metaT_reads / conv_raw_reads$conv_metaG_reads
conv_raw_reads$conv_metaG_reads <- NULL

# Log2 transform the data
cef_raw_reads[,c(1,2)] <- log2(cef_raw_reads[,c(1,2)] + 1)
clinda_raw_reads[,c(1,2)] <- log2(clinda_raw_reads[,c(1,2)] + 1)
strep_raw_reads[,c(1,2)] <- log2(strep_raw_reads[,c(1,2)] + 1)
conv_raw_reads[,1] <- log2(conv_raw_reads[,1] + 1)

# Screen for active transcription in either condition
cef_raw_reads <- subset(cef_raw_reads, cef_raw_reads$cef_630_metaT_reads > 0 | cef_raw_reads$cef_mock_metaT_reads > 0)
clinda_raw_reads <- subset(clinda_raw_reads, clinda_raw_reads$clinda_630_metaT_reads > 0 | clinda_raw_reads$clinda_mock_metaT_reads > 0)
strep_raw_reads <- subset(strep_raw_reads, strep_raw_reads$strep_630_metaT_reads > 0 | strep_raw_reads$strep_mock_metaT_reads > 0)

# Screen for those genes that were able to be annotated
cef_annotated <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('unknown_\\d', cef_raw_reads$gene),]), ]
clinda_annotated <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('unknown_\\d', clinda_raw_reads$gene),]), ]
strep_annotated <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('unknown_\\d', strep_raw_reads$gene),]), ]
conv_annotated <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('unknown_\\d', strep_raw_reads$gene),]), ]

# Also save those that remain unknown
cef_unknown <- cef_raw_reads[rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('unknown_\\d', cef_raw_reads$gene),]), ]
clinda_unknown <- clinda_raw_reads[rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('unknown_\\d', clinda_raw_reads$gene),]), ]
strep_unknown <- strep_raw_reads[rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('unknown_\\d', strep_raw_reads$gene),]), ]
conv_unknown <- strep_raw_reads[rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('unknown_\\d', strep_raw_reads$gene),]), ]

rm(cef_raw_reads, clinda_raw_reads, strep_raw_reads, conv_raw_reads)

#-------------------------------------------------------------------------------------------------------------------------#

# Order data to find distinct patterns of expression
cef_630_top <- cef_annotated[order(cef_annotated$cef_mock_metaT_reads, -cef_annotated$cef_630_metaT_reads),]
cef_mock_top <- cef_annotated[order(cef_annotated$cef_630_metaT_reads, -cef_annotated$cef_mock_metaT_reads),]
clinda_630_top <- clinda_annotated[order(clinda_annotated$clinda_mock_metaT_reads, -clinda_annotated$clinda_630_metaT_reads),]
clinda_mock_top <- clinda_annotated[order(clinda_annotated$clinda_630_metaT_reads, -clinda_annotated$clinda_mock_metaT_reads),]
strep_630_top <- strep_annotated[order(strep_annotated$strep_mock_metaT_reads, -strep_annotated$strep_630_metaT_reads),]
strep_mock_top <- strep_annotated[order(strep_annotated$strep_630_metaT_reads, -strep_annotated$strep_mock_metaT_reads),]

# Screen for largest distinctions between conditions
size <- max(cef_630_top$cef_630_metaT_reads - cef_630_top$cef_mock_metaT_reads) * 0.75
cef_630_top <- subset(cef_630_top, (cef_630_top$cef_630_metaT_reads - cef_630_top$cef_mock_metaT_reads) > size)
size <- max(cef_mock_top$cef_mock_metaT_reads - cef_mock_top$cef_630_metaT_reads) * 0.75
cef_mock_top <- subset(cef_mock_top, (cef_mock_top$cef_mock_metaT_reads - cef_mock_top$cef_630_metaT_reads) > size)
size <- max(clinda_630_top$clinda_630_metaT_reads - clinda_630_top$clinda_mock_metaT_reads) * 0.75
clinda_630_top <- subset(clinda_630_top, (clinda_630_top$clinda_630_metaT_reads - clinda_630_top$clinda_mock_metaT_reads) > size)
size <- max(clinda_mock_top$clinda_mock_metaT_reads - clinda_mock_top$clinda_630_metaT_reads) * 0.75
clinda_mock_top <- subset(clinda_mock_top, (clinda_mock_top$clinda_mock_metaT_reads - clinda_mock_top$clinda_630_metaT_reads) > size)
size <- max(strep_630_top$strep_630_metaT_reads - strep_630_top$strep_mock_metaT_reads) * 0.75
strep_630_top <- subset(strep_630_top, (strep_630_top$strep_630_metaT_reads - strep_630_top$strep_mock_metaT_reads) > size)
size <- max(strep_mock_top$strep_mock_metaT_reads - strep_mock_top$strep_630_metaT_reads) * 0.75
strep_mock_top <- subset(strep_mock_top, (strep_mock_top$strep_mock_metaT_reads - strep_mock_top$strep_630_metaT_reads) > size)
rm(size)

#-------------------------------------------------------------------------------------------------------------------------#





# Check for linear correlation
biglm(strep_annotated$strep_630_metaT_reads ~ strep_annotated$strep_mock_metaT_reads)
biglm(cef_annotated$cef_630_metaT_reads ~ cef_annotated$cef_mock_metaT_reads)
biglm(clinda_annotated$clinda_630_metaT_reads ~ clinda_annotated$clinda_mock_metaT_reads)
biglm(strep_unknown$strep_630_metaT_reads ~ strep_annotated$strep_mock_metaT_reads)
biglm(cef_unknown$cef_630_metaT_reads ~ cef_annotated$cef_mock_metaT_reads)
biglm(clinda_unknown$clinda_630_metaT_reads ~ clinda_annotated$clinda_mock_metaT_reads)




# Calculate the distance of all points from x=y


# Will reveal which genes were most effected by c. diff colonization




#good.dist <- sqrt((good.ord - typ.ord)^2 / 2)



#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_file, width=8, height=12)
layout(matrix(c(1,2,
                3,4,
                5,5), 
              nrow=3, ncol=2, byrow = TRUE))

#-------------------#

# Streptomycin
par(mar=c(4.5, 5, 1, 1), mgp=c(3,0.7,0))
plot(x=strep_annotated$strep_mock_metaT_reads, y=strep_annotated$strep_630_metaT_reads, 
     xlim=c(0,12), ylim=c(0,12), pch=20, cex=1.3, col='gray40', xaxt='n', yaxt='n', xlab='', ylab='')
segments(-2, -2, 14, 14, lty=2)
minor.ticks.axis(1, 12, mn=0, mx=12)
minor.ticks.axis(2, 12, mn=0, mx=12)
mtext('Fold Normalized cDNA Abundance', side=1, padj=2.2, cex=0.7)
mtext('Mock-Infected', side=1, padj=3.5, font=2, cex=0.9)
mtext('Fold Normalized cDNA Abundance', side=2, padj=-2.2, cex=0.7)
mtext(expression(bolditalic('C. difficile')~bold('630-Infected')), side=2, padj=-3.5, font=2, cex=0.9)
legend('topleft', 'Streptomycin-pretreated', bty='n', cex=1.2) 

# Points for pathways of interest
points(x=strep_pathways$strep_mock_metaT_reads, y=strep_pathways$strep_630_metaT_reads, 
       cex=1.9, pch=21, bg=adjustcolor(strep_pathways$colors, alpha=0.6), col='black')

mtext('A', side=2, line=2, las=2, adj=6, padj=-2, cex=1.3)

#-------------------#

# Cefoperazone
par(mar=c(4.5, 5, 1, 1), mgp=c(3,0.7,0))
plot(x=cef_annotated$cef_mock_metaT_reads, y=cef_annotated$cef_630_metaT_reads, 
     xlim=c(0,12), ylim=c(0,12), pch=20, cex=1.3, col='gray40', xaxt='n', yaxt='n', xlab='', ylab='')
segments(-2, -2, 14, 14, lty=2)
minor.ticks.axis(1, 12, mn=0, mx=12)
minor.ticks.axis(2, 12, mn=0, mx=12)
mtext('Fold Normalized cDNA Abundance', side=1, padj=2.2, cex=0.7)
mtext('Mock-Infected', side=1, padj=3.5, font=2, cex=0.9)
mtext('Fold Normalized cDNA Abundance', side=2, padj=-2.2, cex=0.7)
mtext(expression(bolditalic('C. difficile')~bold('630-Infected')), side=2, padj=-3.5, font=2, cex=0.9)
legend('topleft', 'Cefoperazone-pretreated', bty='n', cex=1.2) 

# Points for pathways of interest
points(x=cef_pathways$cef_mock_metaT_reads, y=cef_pathways$cef_630_metaT_reads, 
       cex=1.9, pch=21, bg=adjustcolor(strep_pathways$colors, alpha=0.6), col='black')

#-------------------#

# Clindamycin
par(mar=c(4.5, 5, 1, 1), mgp=c(3,0.7,0))
plot(x=clinda_annotated$clinda_mock_metaT_reads, y=clinda_annotated$clinda_630_metaT_reads, 
     xlim=c(0,12), ylim=c(0,12), pch=20, cex=1.3, col='gray40', xaxt='n', yaxt='n', xlab='', ylab='')
segments(-2, -2, 14, 14, lty=2)
minor.ticks.axis(1, 12, mn=0, mx=12)
minor.ticks.axis(2, 12, mn=0, mx=12)
mtext('Fold Normalized cDNA Abundance', side=1, padj=2.2, cex=0.7)
mtext('Mock-Infected', side=1, padj=3.5, font=2, cex=0.9)
mtext('Fold Normalized cDNA Abundance', side=2, padj=-2.2, cex=0.7)
mtext(expression(bolditalic('C. difficile')~bold('630-Infected')), side=2, padj=-3.5, font=2, cex=0.9)
legend('topleft', 'Clindamycin-pretreated', bty='n', cex=1.2) 

# Points for pathways of interest
points(x=clinda_pathways$clinda_mock_metaT_reads, y=clinda_pathways$clinda_630_metaT_reads, 
       cex=1.9, pch=21, bg=adjustcolor(strep_pathways$colors, alpha=0.6), col='black')

#-------------------#

par(mar=c(1,1,1,1))
plot(0, type='n', axes=FALSE, xlab='', ylab='')
legend('center', legend=c('Amino sugar metabolism', 'Fructose/Mannose metabolism', 'Proline/Glycine metabolism', 'Galactose metabolism'), 
       pt.bg=c('chartreuse3', 'firebrick3', 'darkgoldenrod1', 'darkorchid3'), 
       pch=21, cex=1.6, pt.cex=3, ncol=1)

#-------------------#

# Pathway plot like in Cody's paper - overrepresentation in infected vs mock

par(las=1, mar=c(5,4,1,1), mgp=c(2.5,0.7,0), yaxs='i')
plot(0, type='n', xlab='', ylab='Fold Transcript Abundance Difference', 
     xlim=c(1,15), ylim=c(-5,5), xaxt='n', yaxt='n')
box()
axis(side=2, at=c(-5:5), labels=c(-5:5))
abline(h=0, lty=2, lwd=2, col='gray35')
abline(h=c(-2,2), lty=3, lwd=2, col='firebrick3')
legend('topright', legend=c("630 infected", "Mock infected"), cex=1.2,
       pch=c(21, 21), pt.bg=c("mediumorchid3","mediumseagreen"), bg='white', pt.cex=2)


text(x=seq(2,14,2), y=par()$usr[3]-0.035*(par()$usr[4]-par()$usr[3]), srt=45, adj=1, xpd=TRUE, cex=1.2,
     labels=c('Gene 1', 'Gene 2', 'Gene 3', 'Gene 4', 'Gene 5', 'Gene 6', 'Gene 7'))



mtext('B', side=2, line=2, las=2, adj=6, padj=-2, cex=1.3)





dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up

#for (dep in deps){
#  pkg <- paste('package:', dep, sep='')
#  detach(pkg, character.only = TRUE)}
#rm(list=ls())
#gc()
