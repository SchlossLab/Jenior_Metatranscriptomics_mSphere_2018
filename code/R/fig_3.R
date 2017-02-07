# Set up environment

# Load dependencies
deps <- c('vegan', 'plotrix')
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Set seed for RNG
set.seed(6189)

#----------------#

# Load in functions
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

#----------------#

# Define file names

plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/figures/figure_3.pdf'

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
conv_final_reads_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/conventional.RNA_reads2pangenome.all.norm.remove.annotated.txt'

# Metabolomes
metabolome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolome/metabolomics.tsv'

# Experimental design metadata
metadata_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metadata.tsv'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data

# Metagenomes
cef_metagenome <- read.delim(cef_metagenome_file, sep='\t', header=FALSE, row.names=1)
colnames(cef_metagenome) <- c('cef_metaG_reads', 'ko', 'gene', 'pathway')
clinda_metagenome <- read.delim(clinda_metagenome_file, sep='\t', header=FALSE, row.names=1)
colnames(clinda_metagenome) <- c('clinda_metaG_reads', 'ko', 'gene', 'pathway')
strep_metagenome <- read.delim(strep_metagenome_file, sep='\t', header=FALSE, row.names=1)
colnames(strep_metagenome) <- c('strep_metaG_reads', 'ko', 'gene', 'pathway')
conv_metagenome <- read.delim(conv_metagenome_file, sep='\t', header=FALSE, row.names=1)
colnames(conv_metagenome) <- c('conv_metaG_reads', 'ko', 'gene', 'pathway')
rm(cef_metagenome_file, clinda_metagenome_file, strep_metagenome_file, conv_metagenome_file)

# Metatranscriptomes
cef_630_metatranscriptome <- read.delim(cef_630_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(cef_630_metatranscriptome) <- c('cef_630_metaT_reads', 'ko', 'gene', 'pathway')
cef_mock_metatranscriptome <- read.delim(cef_mock_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(cef_mock_metatranscriptome) <- c('cef_mock_metaT_reads', 'ko', 'gene', 'pathway')
rm(cef_630_metatranscriptome_file, cef_mock_metatranscriptome_file)
clinda_630_metatranscriptome <- read.delim(clinda_630_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(clinda_630_metatranscriptome) <- c('clinda_630_metaT_reads', 'ko', 'gene', 'pathway')
clinda_mock_metatranscriptome <- read.delim(clinda_mock_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(clinda_mock_metatranscriptome) <- c('clinda_mock_metaT_reads', 'ko', 'gene', 'pathway')
rm(clinda_630_metatranscriptome_file, clinda_mock_metatranscriptome_file)
strep_630_metatranscriptome <- read.delim(strep_630_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(strep_630_metatranscriptome) <- c('strep_630_metaT_reads', 'ko', 'gene', 'pathway')
strep_mock_metatranscriptome <- read.delim(strep_mock_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(strep_mock_metatranscriptome) <- c('strep_mock_metaT_reads', 'ko', 'gene', 'pathway')
rm(strep_630_metatranscriptome_file, strep_mock_metatranscriptome_file)
conv_final_reads <- read.delim(conv_final_reads_file, sep='\t', header=FALSE, row.names=1)
colnames(conv_final_reads) <- c('conv_metaT_reads', 'ko', 'gene', 'pathway')
rm(conv_final_reads_file)

# Metabolomes
metabolome <- read.delim(metabolome_file, sep='\t', header=TRUE)
rm(metabolome_file)

# Metadata
metadata <- read.delim(metadata_file, sep='\t', header=T, row.names=1)
rm(metadata_file)

#-------------------------------------------------------------------------------------------------------------------------#

# Format metabolomics


metadata <- metadata[!rownames(metadata) %in% c('CefC5M2'), ] # Remove contaminated sample
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL

rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome <- t(metabolome)
metabolome <- metabolome[!rownames(metabolome) %in% c('CefC5M2'), ]

metabolome_nmds <- metaMDS(metabolome, k=2, trymax=100)$points
metabolome_nmds[,1] <- metabolome_nmds[,1] + 0.15
metabolome_nmds <- clean_merge(metadata, metabolome_nmds)


susceptible <- subset(metabolome_nmds, susceptibility == 'susceptible')
resistant <- subset(metabolome_nmds, susceptibility == 'resistant')
cefoperazone <- subset(metabolome_nmds, abx == 'cefoperazone')
clindamycin <- subset(metabolome_nmds, abx == 'clindamycin')
streptomycin <- subset(metabolome_nmds, abx == 'streptomycin')
germfree <- subset(metabolome_nmds, abx == 'germfree')
untreated <- subset(metabolome_nmds, abx == 'none')


pdf(file='~/Desktop/metabolome_nmds.pdf', width=6, height=6)

dev.off()



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
conv_raw_reads <- clean_merge(conv_metagenome, conv_final_reads)

rm(cef_metagenome, clinda_metagenome, strep_metagenome, conv_metagenome,
   cef_630_metatranscriptome, cef_mock_metatranscriptome, clinda_630_metatranscriptome,
   clinda_mock_metatranscriptome, strep_630_metatranscriptome, strep_mock_metatranscriptome,
   conv_final_reads)

#-------------------------------------------------------------------------------------------------------------------------#

# Remove residual C. difficile 630 mappings
cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% grep('cdf:CD630', rownames(cef_raw_reads)), ]
clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% grep('cdf:CD630', rownames(clinda_raw_reads)), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% grep('cdf:CD630', rownames(strep_raw_reads)), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% grep('cdf:CD630', rownames(conv_raw_reads)), ]

# Remove genes with no metagenomic coverage
cef_raw_reads <- subset(cef_raw_reads, cef_metaG_reads != 0)
clinda_raw_reads <- subset(clinda_raw_reads, clinda_metaG_reads != 0)
strep_raw_reads <- subset(strep_raw_reads, strep_metaG_reads != 0)
conv_raw_reads <- subset(conv_raw_reads, conv_metaG_reads != 0)

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
cef_final_reads <- subset(cef_raw_reads, cef_raw_reads$cef_630_metaT_reads > 0 | cef_raw_reads$cef_mock_metaT_reads > 0)
clinda_final_reads <- subset(clinda_raw_reads, clinda_raw_reads$clinda_630_metaT_reads > 0 | clinda_raw_reads$clinda_mock_metaT_reads > 0)
strep_final_reads <- subset(strep_raw_reads, strep_raw_reads$strep_630_metaT_reads > 0 | strep_raw_reads$strep_mock_metaT_reads > 0)
conv_final_reads <- conv_raw_reads
rm(cef_raw_reads, clinda_raw_reads, strep_raw_reads, conv_raw_reads)

#-------------------------------------------------------------------------------------------------------------------------#

# Calculate the distance of all points from x=y
# Will reveal which genes were most effected by c. diff colonization




#good.dist <- sqrt((good.ord - typ.ord)^2 / 2)



#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_file, width=8.5, height=11)
layout(matrix(c(1,2,
                3,3,
                3,3),
              nrow=3, ncol=2, byrow=TRUE))

#-------------------#

# Metabolomics alone

par(mar=c(3,4,1,1), las=1, mgp=c(2,0.75,0), xaxs='i', yaxs='i')
plot(x=metabolome_nmds$MDS1, y=metabolome_nmds$MDS2, xlim=c(-0.4,0.4), ylim=c(-0.3,0.3),
     xlab='NMDS axis 1', ylab='NMDS axis 2', xaxt='n', yaxt='n', pch=19, cex=0.2)
axis(side=1, at=seq(-0.4,0.4,0.1), labels=seq(-0.4,0.4,0.1))
axis(side=2, at=seq(-0.3,0.3,0.1), labels=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3))
points(x=cefoperazone$MDS1, y=cefoperazone$MDS2, bg='chartreuse3', pch=21, cex=1.7, lwd=1.2)
points(x=clindamycin$MDS1, y=clindamycin$MDS2, bg='blue2', pch=21, cex=1.7, lwd=1.2)
points(x=streptomycin$MDS1, y=streptomycin$MDS2, bg='firebrick1', pch=21, cex=1.7, lwd=1.2)
points(x=germfree$MDS1, y=germfree$MDS2, bg='gold1', pch=21, cex=1.7, lwd=1.2)
points(x=untreated$MDS1, y=untreated$MDS2, bg='azure2', pch=21, cex=1.7, lwd=1.2)
draw.ellipse(x=0.19, y=-0.01, a=0.28, b=0.17, angle=-60, lty=2, lwd=2) # susceptible
text(x=0.25, y=0.25, labels='Susceptible', font=2)
draw.ellipse(x=-0.27, y=-0.13, a=0.15, b=0.1, angle=-100, lty=2, lwd=2) # resistant
text(x=-0.27, y=0.04, labels='Resistant', font=2)
legend('topleft', legend=c('Untreated (SPF)','Streptomycin-treated (SPF)','Cefoperzone-treated (SPF)','Clindamycin-treated (SPF)','Germfree'), 
       pt.bg=c('azure2','firebrick1','chartreuse3','blue2','gold1'), 
       pch=21, pt.cex=1.7, bty='n')

mtext('A', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.3)

#-------------------#

plot(0, type='n', axes=FALSE, xlab='', ylab='')
mtext('B', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.3)

#-------------------#

plot(0, type='n', axes=FALSE, xlab='', ylab='')
# Heatmap or correlation somehow...
mtext('C', side=2, line=2, las=2, adj=1.5, padj=-10, cex=1.3)


dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up

#for (dep in deps){
#  pkg <- paste('package:', dep, sep='')
#  detach(pkg, character.only = TRUE)
#}
#rm(list=ls())
#gc()

