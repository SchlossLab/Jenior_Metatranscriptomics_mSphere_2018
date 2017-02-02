# Set up environment

# Load dependencies
deps <- c('vegan', 'shape')
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Set seed for RNG
set.seed(6189)

#----------------#

# Define functions

# Neatly merge 2 matices with shared row names
clean_merge <- function(data_1, data_2){
  
  clean_merged <- merge(data_1, data_2, by = 'row.names')
  rownames(clean_merged) <- clean_merged$Row.names
  clean_merged$Row.names <- NULL
  
  return(clean_merged)
}

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

# Extract specific pathway annotations


# Find more specific pathways or genes...



cef1 <- cbind(cef_final_reads[grep('Amino_sugar', cef_final_reads$pathway), ][,c(1,2)], rep('amino sugars',length(grep('Amino_sugar', cef_final_reads$pathway))), rep('chartreuse3',length(grep('Amino_sugar', cef_final_reads$pathway))))
colnames(cef1) <- c('cef_mock_metaT_reads', 'cef_630_metaT_reads', 'pathways', 'colors')
cef2 <- cbind(cef_final_reads[grep('Fructose', cef_final_reads$pathway), ][,c(1,2)], rep('amino sugars',length(grep('Fructose', cef_final_reads$pathway))), rep('firebrick3',length(grep('Fructose', cef_final_reads$pathway))))
colnames(cef2) <- c('cef_mock_metaT_reads', 'cef_630_metaT_reads', 'pathways', 'colors')
cef3 <- cbind(cef_final_reads[grep('proline', cef_final_reads$pathway), ][,c(1,2)], rep('proline',length(grep('proline', cef_final_reads$pathway))), rep('darkgoldenrod1',length(grep('proline', cef_final_reads$pathway))))
colnames(cef3) <- c('cef_mock_metaT_reads', 'cef_630_metaT_reads', 'pathways', 'colors')
cef4 <- cbind(cef_final_reads[grep('Glycine', cef_final_reads$pathway), ][,c(1,2)], rep('glycine',length(grep('Glycine', cef_final_reads$pathway))), rep('darkgoldenrod1',length(grep('Glycine', cef_final_reads$pathway))))
colnames(cef4) <- c('cef_mock_metaT_reads', 'cef_630_metaT_reads', 'pathways', 'colors')
cef5 <- cbind(cef_final_reads[grep('Galactose', cef_final_reads$pathway), ][,c(1,2)], rep('galactose',length(grep('Galactose', cef_final_reads$pathway))), rep('chartreuse3',length(grep('Galactose', cef_final_reads$pathway))))
colnames(cef5) <- c('cef_mock_metaT_reads', 'cef_630_metaT_reads', 'pathways', 'colors')
cef_pathways <- rbind(cef1, cef2, cef3, cef4, cef5)
rm(cef1, cef2, cef3, cef4, cef5)

clinda1 <- cbind(clinda_final_reads[grep('Amino_sugar', clinda_final_reads$pathway), ][,c(1,2)], rep('amino sugars',length(grep('Amino_sugar', clinda_final_reads$pathway))), rep('chartreuse3',length(grep('Amino_sugar', clinda_final_reads$pathway))))
colnames(clinda1) <- c('clinda_mock_metaT_reads', 'clinda_630_metaT_reads', 'pathways', 'colors')
clinda2 <- cbind(clinda_final_reads[grep('Fructose', clinda_final_reads$pathway), ][,c(1,2)], rep('amino sugars',length(grep('Fructose', clinda_final_reads$pathway))), rep('firebrick3',length(grep('Fructose', clinda_final_reads$pathway))))
colnames(clinda2) <- c('clinda_mock_metaT_reads', 'clinda_630_metaT_reads', 'pathways', 'colors')
clinda3 <- cbind(clinda_final_reads[grep('proline', clinda_final_reads$pathway), ][,c(1,2)], rep('proline',length(grep('proline', clinda_final_reads$pathway))), rep('darkgoldenrod1',length(grep('proline', clinda_final_reads$pathway))))
colnames(clinda3) <- c('clinda_mock_metaT_reads', 'clinda_630_metaT_reads', 'pathways', 'colors')
clinda4 <- cbind(clinda_final_reads[grep('Glycine', clinda_final_reads$pathway), ][,c(1,2)], rep('glycine',length(grep('Glycine', clinda_final_reads$pathway))), rep('darkgoldenrod1',length(grep('Glycine', clinda_final_reads$pathway))))
colnames(clinda4) <- c('clinda_mock_metaT_reads', 'clinda_630_metaT_reads', 'pathways', 'colors')
clinda5 <- cbind(clinda_final_reads[grep('Galactose', clinda_final_reads$pathway), ][,c(1,2)], rep('galactose',length(grep('Galactose', clinda_final_reads$pathway))), rep('chartreuse3',length(grep('Galactose', clinda_final_reads$pathway))))
colnames(clinda5) <- c('clinda_mock_metaT_reads', 'clinda_630_metaT_reads', 'pathways', 'colors')
clinda_pathways <- rbind(clinda1, clinda2, clinda3, clinda4, clinda5)
rm(clinda1, clinda2, clinda3, clinda4, clinda5)

strep1 <- cbind(strep_final_reads[grep('Amino_sugar', strep_final_reads$pathway), ][,c(1,2)], rep('amino sugars',length(grep('Amino_sugar', strep_final_reads$pathway))), rep('chartreuse3',length(grep('Amino_sugar', strep_final_reads$pathway))))
colnames(strep1) <- c('strep_mock_metaT_reads', 'strep_630_metaT_reads', 'pathways', 'colors')
strep2 <- cbind(strep_final_reads[grep('Fructose', strep_final_reads$pathway), ][,c(1,2)], rep('amino sugars',length(grep('Fructose', strep_final_reads$pathway))), rep('firebrick3',length(grep('Fructose', strep_final_reads$pathway))))
colnames(strep2) <- c('strep_mock_metaT_reads', 'strep_630_metaT_reads', 'pathways', 'colors')
strep3 <- cbind(strep_final_reads[grep('proline', strep_final_reads$pathway), ][,c(1,2)], rep('proline',length(grep('proline', strep_final_reads$pathway))), rep('darkgoldenrod1',length(grep('proline', strep_final_reads$pathway))))
colnames(strep3) <- c('strep_mock_metaT_reads', 'strep_630_metaT_reads', 'pathways', 'colors')
strep4 <- cbind(strep_final_reads[grep('Glycine', strep_final_reads$pathway), ][,c(1,2)], rep('glycine',length(grep('Glycine', strep_final_reads$pathway))), rep('darkgoldenrod1',length(grep('Glycine', strep_final_reads$pathway))))
colnames(strep4) <- c('strep_mock_metaT_reads', 'strep_630_metaT_reads', 'pathways', 'colors')
strep5 <- cbind(strep_final_reads[grep('Galactose', strep_final_reads$pathway), ][,c(1,2)], rep('galactose',length(grep('Galactose', strep_final_reads$pathway))), rep('chartreuse3',length(grep('Galactose', strep_final_reads$pathway))))
colnames(strep5) <- c('strep_mock_metaT_reads', 'strep_630_metaT_reads', 'pathways', 'colors')
strep_pathways <- rbind(strep1, strep2, strep3, strep4, strep5)
rm(strep1, strep2, strep3, strep4, strep5)

conv1 <- cbind(conv_final_reads[grep('Amino_sugar', conv_final_reads$pathway), ][,1], rep('amino sugars',length(grep('Amino_sugar', conv_final_reads$pathway))), rep('chartreuse3',length(grep('Amino_sugar', conv_final_reads$pathway))))
colnames(conv1) <- c('conv_metaT_reads', 'pathways', 'colors')
conv2 <- cbind(conv_final_reads[grep('Fructose', conv_final_reads$pathway), ][,1], rep('amino sugars',length(grep('Fructose', conv_final_reads$pathway))), rep('firebrick3',length(grep('Fructose', conv_final_reads$pathway))))
colnames(conv2) <- c('conv_metaT_reads', 'pathways', 'colors')
conv3 <- cbind(conv_final_reads[grep('proline', conv_final_reads$pathway), ][,1], rep('proline',length(grep('proline', conv_final_reads$pathway))), rep('darkgoldenrod1',length(grep('proline', conv_final_reads$pathway))))
colnames(conv3) <- c('conv_metaT_reads', 'pathways', 'colors')
conv4 <- cbind(conv_final_reads[grep('Glycine', conv_final_reads$pathway), ][,1], rep('glycine',length(grep('Glycine', conv_final_reads$pathway))), rep('darkgoldenrod1',length(grep('Glycine', conv_final_reads$pathway))))
colnames(conv4) <- c('conv_metaT_reads', 'pathways', 'colors')
conv5 <- cbind(conv_final_reads[grep('Galactose', conv_final_reads$pathway), ][,1], rep('galactose',length(grep('Galactose', conv_final_reads$pathway))), rep('chartreuse3',length(grep('Galactose', conv_final_reads$pathway))))
colnames(conv5) <- c('conv_metaT_reads', 'pathways', 'colors')
conv_pathways <- as.data.frame(rbind(conv1, conv2, conv3, conv4, conv5))
rm(conv1, conv2, conv3, conv4, conv5)

# Remove annotated points from general points
cef_final_reads <- cef_final_reads[!rownames(cef_final_reads) %in% rownames(cef_pathways), ]
clinda_final_reads <- clinda_final_reads[!rownames(clinda_final_reads) %in% rownames(clinda_pathways), ]
strep_final_reads <- strep_final_reads[!rownames(strep_final_reads) %in% rownames(strep_pathways), ]
conv_final_reads <- conv_final_reads[!rownames(conv_final_reads) %in% rownames(conv_pathways), ]

#-------------------------------------------------------------------------------------------------------------------------#

# Calculate the distance of all points from x=y
# Will reveal which genes were most effected by c. diff colonization




#good.dist <- sqrt((good.ord - typ.ord)^2 / 2)



#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_file, width=5, height=12)
layout(matrix(c(1,
                2,
                3),
              nrow=3, ncol=1, byrow = TRUE))

#-------------------#

# Heatmap or correlation somehow...


dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up

#for (dep in deps){
#  pkg <- paste('package:', dep, sep='')
#  detach(pkg, character.only = TRUE)
#}
#rm(list=ls())
#gc()

