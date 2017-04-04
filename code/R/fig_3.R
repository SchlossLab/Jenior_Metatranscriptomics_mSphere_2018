
# Set up environment
#rm(list=ls())
#gc()

# Load dependencies
deps <- c('wesanderson', 'vegan', 'shape')
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Set seed for RNG
set.seed(8619)

# Load in functions
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define files

# Metagenomes
cef_metagenome <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Cefoperazone.DNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_metagenome <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Clindamycin.DNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_metagenome <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Streptomycin.DNA_reads2pangenome.all.norm.remove.annotated.txt'
conv_metagenome <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Conventional.DNA_reads2pangenome.all.norm.remove.annotated.txt'

# Metatranscriptomes
cef_630_metatranscriptome <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/cefoperazone_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'
cef_mock_metatranscriptome <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/cefoperazone_mock.RNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_630_metatranscriptome <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/clindamycin_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_mock_metatranscriptome <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/clindamycin_mock.RNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_630_metatranscriptome <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/streptomycin_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_mock_metatranscriptome <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/streptomycin_mock.RNA_reads2pangenome.all.norm.remove.annotated.txt'
conv_metatranscriptome <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/conventional.RNA_reads2pangenome.all.norm.remove.annotated.txt'

# KEGG taxonomy IDs
kegg_tax <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/kegg_taxonomy.tsv'

# Output plot
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/figures/figure_3.pdf'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data

# Metagenomes
cef_metagenome <- read.delim(cef_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(cef_metagenome) <- c('cef_metaG_reads', 'ko', 'gene', 'pathway')
clinda_metagenome <- read.delim(clinda_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(clinda_metagenome) <- c('clinda_metaG_reads', 'ko', 'gene', 'pathway')
strep_metagenome <- read.delim(strep_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(strep_metagenome) <- c('strep_metaG_reads', 'ko', 'gene', 'pathway')
conv_metagenome <- read.delim(conv_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(conv_metagenome) <- c('conv_metaG_reads', 'ko', 'gene', 'pathway')

# Metatranscriptomes
cef_630_metatranscriptome <- read.delim(cef_630_metatranscriptome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(cef_630_metatranscriptome) <- c('cef_630_metaT_reads', 'ko', 'gene', 'pathway')
cef_mock_metatranscriptome <- read.delim(cef_mock_metatranscriptome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(cef_mock_metatranscriptome) <- c('cef_mock_metaT_reads', 'ko', 'gene', 'pathway')
clinda_630_metatranscriptome <- read.delim(clinda_630_metatranscriptome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(clinda_630_metatranscriptome) <- c('clinda_630_metaT_reads', 'ko', 'gene', 'pathway')
clinda_mock_metatranscriptome <- read.delim(clinda_mock_metatranscriptome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(clinda_mock_metatranscriptome) <- c('clinda_mock_metaT_reads', 'ko', 'gene', 'pathway')
strep_630_metatranscriptome <- read.delim(strep_630_metatranscriptome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(strep_630_metatranscriptome) <- c('strep_630_metaT_reads', 'ko', 'gene', 'pathway')
strep_mock_metatranscriptome <- read.delim(strep_mock_metatranscriptome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(strep_mock_metatranscriptome) <- c('strep_mock_metaT_reads', 'ko', 'gene', 'pathway')
conv_metatranscriptome <- read.delim(conv_metatranscriptome, sep='\t', header=FALSE, row.names=1)
colnames(conv_metatranscriptome) <- c('conv_metaT_reads', 'ko', 'gene', 'pathway')

# KEGG organism file
kegg_tax <- read.delim(kegg_tax, sep='\t', header=TRUE)
kegg_tax[] <- lapply(kegg_tax, as.character)

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

# Remove residual C. difficile mappings
cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdf:CD630', rownames(cef_raw_reads)),]), ]
cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdc:CD196', rownames(cef_raw_reads)),]), ]
cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdf:CDP07', rownames(cef_raw_reads)),]), ]
cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdf:CDP02', rownames(cef_raw_reads)),]), ]
cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdf:CDP06', rownames(cef_raw_reads)),]), ]
cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdf:CDP03', rownames(cef_raw_reads)),]), ]
cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdg:CDBI1_08660', rownames(cef_raw_reads)),]), ]
cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdg:CDBI1_06165', rownames(cef_raw_reads)),]), ]
cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdg:CDBI1_09680', rownames(cef_raw_reads)),]), ]
cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdg:CDBI1_02800', rownames(cef_raw_reads)),]), ]
cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdg:CDBI1_13505', rownames(cef_raw_reads)),]), ]
cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdl:CDR20291_1755', rownames(cef_raw_reads)),]), ]

clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdf:CD630', rownames(clinda_raw_reads)),]), ]
clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdc:CD196', rownames(clinda_raw_reads)),]), ]
clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdf:CDP07', rownames(clinda_raw_reads)),]), ]
clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdf:CDP02', rownames(clinda_raw_reads)),]), ]
clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdf:CDP06', rownames(clinda_raw_reads)),]), ]
clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdf:CDP03', rownames(clinda_raw_reads)),]), ]
clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdg:CDBI1_08660', rownames(clinda_raw_reads)),]), ]
clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdg:CDBI1_06165', rownames(clinda_raw_reads)),]), ]
clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdg:CDBI1_09680', rownames(clinda_raw_reads)),]), ]
clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdg:CDBI1_02800', rownames(clinda_raw_reads)),]), ]
clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdg:CDBI1_13505', rownames(clinda_raw_reads)),]), ]
clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdl:CDR20291_1755', rownames(clinda_raw_reads)),]), ]

strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdf:CD630', rownames(strep_raw_reads)),]), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdc:CD196', rownames(strep_raw_reads)),]), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdf:CD630', rownames(strep_raw_reads)),]), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdc:CD196', rownames(strep_raw_reads)),]), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdf:CDP07', rownames(strep_raw_reads)),]), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdf:CDP02', rownames(strep_raw_reads)),]), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdf:CDP06', rownames(strep_raw_reads)),]), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdf:CDP03', rownames(strep_raw_reads)),]), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdg:CDBI1_08660', rownames(strep_raw_reads)),]), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdg:CDBI1_06165', rownames(strep_raw_reads)),]), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdg:CDBI1_09680', rownames(strep_raw_reads)),]), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdg:CDBI1_02800', rownames(strep_raw_reads)),]), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdg:CDBI1_13505', rownames(strep_raw_reads)),]), ]
strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdl:CDR20291_1755', rownames(strep_raw_reads)),]), ]

conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdf:CD630', rownames(conv_raw_reads)),]), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdc:CD196', rownames(conv_raw_reads)),]), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdf:CD630', rownames(conv_raw_reads)),]), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdc:CD196', rownames(conv_raw_reads)),]), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdf:CDP07', rownames(conv_raw_reads)),]), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdf:CDP02', rownames(conv_raw_reads)),]), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdf:CDP06', rownames(conv_raw_reads)),]), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdf:CDP03', rownames(conv_raw_reads)),]), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdg:CDBI1_08660', rownames(conv_raw_reads)),]), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdg:CDBI1_06165', rownames(conv_raw_reads)),]), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdg:CDBI1_09680', rownames(conv_raw_reads)),]), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdg:CDBI1_02800', rownames(conv_raw_reads)),]), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdg:CDBI1_13505', rownames(conv_raw_reads)),]), ]
conv_raw_reads <- conv_raw_reads[!rownames(conv_raw_reads) %in% rownames(conv_raw_reads[grep('cdl:CDR20291_1755', rownames(conv_raw_reads)),]), ]

# Remove genes with no metagenomic coverage
cef_raw_reads <- subset(cef_raw_reads, cef_metaG_reads != 0)
clinda_raw_reads <- subset(clinda_raw_reads, clinda_metaG_reads != 0)
strep_raw_reads <- subset(strep_raw_reads, strep_metaG_reads != 0)
conv_raw_reads <- subset(conv_raw_reads, conv_metaG_reads != 0)

# Rarefy read abundances
size <- round(min(colSums(cef_raw_reads[,c(1:3)]))*0.9) # Determine subsample level
cef_raw_reads$cef_metaG_reads <- t(rrarefy(cef_raw_reads$cef_metaG_reads, sample=size)) + 1
cef_raw_reads$cef_630_metaT_reads <- t(rrarefy(cef_raw_reads$cef_630_metaT_reads, sample=size))
cef_raw_reads$cef_mock_metaT_reads <- t(rrarefy(cef_raw_reads$cef_mock_metaT_reads, sample=size))
cef_normalized_reads <- cef_raw_reads
rm(cef_raw_reads)
size <- round(min(colSums(clinda_raw_reads[,c(1:3)]))*0.9) # Determine subsample level
clinda_raw_reads$clinda_metaG_reads <- t(rrarefy(clinda_raw_reads$clinda_metaG_reads, sample=size)) + 1
clinda_raw_reads$clinda_630_metaT_reads <- t(rrarefy(clinda_raw_reads$clinda_630_metaT_reads, sample=size))
clinda_raw_reads$clinda_mock_metaT_reads <- t(rrarefy(clinda_raw_reads$clinda_mock_metaT_reads, sample=size))
clinda_normalized_reads <- clinda_raw_reads
rm(clinda_raw_reads)
size <- round(min(colSums(strep_raw_reads[,c(1:3)]))*0.9) # Determine subsample level
strep_raw_reads$strep_metaG_reads <- t(rrarefy(strep_raw_reads$strep_metaG_reads, sample=size)) + 1
strep_raw_reads$strep_630_metaT_reads <- t(rrarefy(strep_raw_reads$strep_630_metaT_reads, sample=size))
strep_raw_reads$strep_mock_metaT_reads <- t(rrarefy(strep_raw_reads$strep_mock_metaT_reads, sample=size))
strep_normalized_reads <- strep_raw_reads
rm(strep_raw_reads)
size <- round(min(colSums(conv_raw_reads[,c(1:2)]))*0.9) # Determine subsample level
conv_raw_reads$conv_metaG_reads <- t(rrarefy(conv_raw_reads$conv_metaG_reads, sample=size)) + 1
conv_raw_reads$conv_metaT_reads <- t(rrarefy(conv_raw_reads$conv_metaT_reads, sample=size))
conv_normalized_reads <- conv_raw_reads
rm(conv_raw_reads)
rm(size)

# Normalize metatranscriptomes to metagenomic coverage
cef_normalized_reads$cef_630_metaT_reads <- cef_normalized_reads$cef_630_metaT_reads / cef_normalized_reads$cef_metaG_reads
cef_normalized_reads$cef_mock_metaT_reads <- cef_normalized_reads$cef_mock_metaT_reads / cef_normalized_reads$cef_metaG_reads
cef_normalized_reads$cef_metaG_reads <- NULL
clinda_normalized_reads$clinda_630_metaT_reads <- clinda_normalized_reads$clinda_630_metaT_reads / clinda_normalized_reads$clinda_metaG_reads
clinda_normalized_reads$clinda_mock_metaT_reads <- clinda_normalized_reads$clinda_mock_metaT_reads / clinda_normalized_reads$clinda_metaG_reads
clinda_normalized_reads$clinda_metaG_reads <- NULL
strep_normalized_reads$strep_630_metaT_reads <- strep_normalized_reads$strep_630_metaT_reads / strep_normalized_reads$strep_metaG_reads
strep_normalized_reads$strep_mock_metaT_reads <- strep_normalized_reads$strep_mock_metaT_reads / strep_normalized_reads$strep_metaG_reads
strep_normalized_reads$strep_metaG_reads <- NULL
conv_normalized_reads$conv_metaT_reads <- conv_normalized_reads$conv_metaT_reads / conv_normalized_reads$conv_metaG_reads
conv_normalized_reads$conv_metaG_reads <- NULL

# Log2 transform the data
cef_normalized_reads[,c(1,2)] <- log2(cef_normalized_reads[,c(1,2)] + 1)
clinda_normalized_reads[,c(1,2)] <- log2(clinda_normalized_reads[,c(1,2)] + 1)
strep_normalized_reads[,c(1,2)] <- log2(strep_normalized_reads[,c(1,2)] + 1)
conv_normalized_reads[,1] <- log2(conv_normalized_reads[,1] + 1)

# Screen for active transcription in either condition
cef_normalized_reads <- subset(cef_normalized_reads, cef_normalized_reads$cef_630_metaT_reads > 0 | cef_normalized_reads$cef_mock_metaT_reads > 0)
clinda_normalized_reads <- subset(clinda_normalized_reads, clinda_normalized_reads$clinda_630_metaT_reads > 0 | clinda_normalized_reads$clinda_mock_metaT_reads > 0)
strep_normalized_reads <- subset(strep_normalized_reads, strep_normalized_reads$strep_630_metaT_reads > 0 | strep_normalized_reads$strep_mock_metaT_reads > 0)

# Screen for those genes that were able to be annotated
cef_annotated <- cef_normalized_reads[!rownames(cef_normalized_reads) %in% rownames(cef_normalized_reads[grep('unknown_\\d', cef_normalized_reads$gene),]), ]
clinda_annotated <- clinda_normalized_reads[!rownames(clinda_normalized_reads) %in% rownames(clinda_normalized_reads[grep('unknown_\\d', clinda_normalized_reads$gene),]), ]
strep_annotated <- strep_normalized_reads[!rownames(strep_normalized_reads) %in% rownames(strep_normalized_reads[grep('unknown_\\d', strep_normalized_reads$gene),]), ]
#conv_annotated <- strep_normalized_reads[!rownames(strep_normalized_reads) %in% rownames(strep_normalized_reads[grep('unknown_\\d', strep_normalized_reads$gene),]), ]

# Also save those that remain unknown
#cef_unknown <- cef_normalized_reads[rownames(cef_normalized_reads) %in% rownames(cef_normalized_reads[grep('unknown_\\d', cef_normalized_reads$gene),]), ]
#clinda_unknown <- clinda_normalized_reads[rownames(clinda_normalized_reads) %in% rownames(clinda_normalized_reads[grep('unknown_\\d', clinda_normalized_reads$gene),]), ]
#strep_unknown <- strep_normalized_reads[rownames(strep_normalized_reads) %in% rownames(strep_normalized_reads[grep('unknown_\\d', strep_normalized_reads$gene),]), ]
#conv_unknown <- strep_normalized_reads[rownames(strep_normalized_reads) %in% rownames(strep_normalized_reads[grep('unknown_\\d', strep_normalized_reads$gene),]), ]
rm(cef_normalized_reads, clinda_normalized_reads, strep_normalized_reads, conv_normalized_reads)

#-------------------------------------------------------------------------------------------------------------------------#

# Calculate correlation coefficients
strep_corr <- as.character(round(cor.test(strep_annotated[,1], strep_annotated[,2], method='spearman', exact=FALSE)$estimate, digits=3))
cef_corr <- as.character(round(cor.test(cef_annotated[,1], cef_annotated[,2], method='spearman', exact=FALSE)$estimate, digits=3))
clinda_corr <- as.character(round(cor.test(clinda_annotated[,1], clinda_annotated[,2], method='spearman', exact=FALSE)$estimate, digits=3))

# Using previously defined lines, find outliers to y = x
strep_630_outliers <- subset(strep_annotated, strep_annotated$strep_mock_metaT_reads > strep_annotated$strep_630_metaT_reads + 2)
strep_mock_outliers <- subset(strep_annotated, strep_annotated$strep_mock_metaT_reads < strep_annotated$strep_630_metaT_reads - 2)
cef_630_outliers <- subset(cef_annotated, cef_annotated$cef_mock_metaT_reads > cef_annotated$cef_630_metaT_reads + 2)
cef_mock_outliers <- subset(cef_annotated, cef_annotated$cef_mock_metaT_reads < cef_annotated$cef_630_metaT_reads - 2)
clinda_630_outliers <- subset(clinda_annotated, clinda_annotated$clinda_mock_metaT_reads > clinda_annotated$clinda_630_metaT_reads + 2)
clinda_mock_outliers <- subset(clinda_annotated, clinda_annotated$clinda_mock_metaT_reads < clinda_annotated$clinda_630_metaT_reads - 2)

#-------------------------------------------------------------------------------------------------------------------------#

# Break down outliers into taxonomic groups

# Split out origin organism code
strep_630_outliers <- get_kegg_org(strep_630_outliers)
strep_mock_outliers <- get_kegg_org(strep_mock_outliers)
cef_630_outliers <- get_kegg_org(cef_630_outliers)
cef_mock_outliers <- get_kegg_org(cef_mock_outliers)
clinda_630_outliers <- get_kegg_org(clinda_630_outliers)
clinda_mock_outliers <- get_kegg_org(clinda_mock_outliers)

# Drop levels
strep_630_outliers[] <- lapply(strep_630_outliers, as.character)
strep_mock_outliers[] <- lapply(strep_mock_outliers, as.character)
cef_630_outliers[] <- lapply(cef_630_outliers, as.character)
cef_mock_outliers[] <- lapply(cef_mock_outliers, as.character)
clinda_630_outliers[] <- lapply(clinda_630_outliers, as.character)
clinda_mock_outliers[] <- lapply(clinda_mock_outliers, as.character)

# Save KEGG ID names
strep_630_outliers$kegg_id <- rownames(strep_630_outliers)
strep_mock_outliers$kegg_id <- rownames(strep_mock_outliers)
cef_630_outliers$kegg_id <- rownames(cef_630_outliers)
cef_mock_outliers$kegg_id <- rownames(cef_mock_outliers)
clinda_630_outliers$kegg_id <- rownames(clinda_630_outliers)
clinda_mock_outliers$kegg_id <- rownames(clinda_mock_outliers)

# Merge with KEGG taxonomy
strep_630_outliers <- merge(x=strep_630_outliers, y=kegg_tax, by.x='org_code', by.y='org_code')
strep_630_outliers$org_code <- NULL
strep_630_outliers$strep_630_metaT_reads <- as.numeric(strep_630_outliers$strep_630_metaT_reads)
strep_630_outliers$strep_mock_metaT_reads <- as.numeric(strep_630_outliers$strep_mock_metaT_reads)
strep_mock_outliers <- merge(x=strep_mock_outliers, y=kegg_tax, by.x='org_code', by.y='org_code')
strep_mock_outliers$org_code <- NULL
strep_mock_outliers$strep_630_metaT_reads <- as.numeric(strep_mock_outliers$strep_630_metaT_reads)
strep_mock_outliers$strep_mock_metaT_reads <- as.numeric(strep_mock_outliers$strep_mock_metaT_reads)
cef_630_outliers <- merge(x=cef_630_outliers, y=kegg_tax, by.x='org_code', by.y='org_code')
cef_630_outliers$org_code <- NULL
cef_630_outliers$cef_630_metaT_reads <- as.numeric(cef_630_outliers$cef_630_metaT_reads)
cef_630_outliers$cef_mock_metaT_reads <- as.numeric(cef_630_outliers$cef_mock_metaT_reads)
cef_mock_outliers <- merge(x=cef_mock_outliers, y=kegg_tax, by.x='org_code', by.y='org_code')
cef_mock_outliers$org_code <- NULL
cef_mock_outliers$cef_630_metaT_reads <- as.numeric(cef_mock_outliers$cef_630_metaT_reads)
cef_mock_outliers$cef_mock_metaT_reads <- as.numeric(cef_mock_outliers$cef_mock_metaT_reads)
clinda_630_outliers <- merge(x=clinda_630_outliers, y=kegg_tax, by.x='org_code', by.y='org_code')
clinda_630_outliers$org_code <- NULL
clinda_630_outliers$clinda_630_metaT_reads <- as.numeric(clinda_630_outliers$clinda_630_metaT_reads)
clinda_630_outliers$clinda_mock_metaT_reads <- as.numeric(clinda_630_outliers$clinda_mock_metaT_reads)
clinda_mock_outliers <- merge(x=clinda_mock_outliers, y=kegg_tax, by.x='org_code', by.y='org_code')
clinda_mock_outliers$org_code <- NULL
clinda_mock_outliers$clinda_630_metaT_reads <- as.numeric(clinda_mock_outliers$clinda_630_metaT_reads)
clinda_mock_outliers$clinda_mock_metaT_reads <- as.numeric(clinda_mock_outliers$clinda_mock_metaT_reads)
rm(kegg_tax)

# Define colors based on genus







#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_file, width=10, height=10)
layout(matrix(c(1,2,
                3,4), 
              nrow=2, ncol=2, byrow = TRUE))
par(mar=c(3.5, 3.5, 1, 1), mgp=c(3,0.7,0))

#-------------------#

# Streptomycin
plot(0, type='n', xlim=c(0,12), ylim=c(0,12), pch=20, xaxt='n', yaxt='n', xlab='', ylab='')
filledrectangle(wx=20, wy=2.8, col='gray80', mid=c(6,6), angle=45)
box()
points(x=strep_annotated$strep_mock_metaT_reads, y=strep_annotated$strep_630_metaT_reads, pch=20, cex=1.3, col='gray40')
segments(-2, -2, 14, 14, lty=2)
minor.ticks.axis(1, 12, mn=0, mx=12)
minor.ticks.axis(2, 12, mn=0, mx=12)
mtext('Fold Normalized cDNA Abundance', side=1, padj=2.8, cex=0.7)
mtext('Mock-Infected', side=1, padj=3.7, font=2, cex=0.9)
mtext('Fold Normalized cDNA Abundance', side=2, padj=-2.7, cex=0.7)
mtext(expression(bolditalic('C. difficile')~bold('630-Infected')), side=2, padj=-3.5, font=2, cex=0.9)
legend('topleft', c('Streptomycin-pretreated', as.expression(bquote(paste(italic('rho'),' = ',.(strep_corr))))), bty='n', cex=1.2, text.col=c(wes_palette("FantasticFox")[1],'black'))


points(x=strep_630_outliers$strep_mock_metaT_reads, y=strep_630_outliers$strep_630_metaT_reads, cex=1.9, pch=21, bg='red', col='black')
points(x=strep_mock_outliers$strep_mock_metaT_reads, y=strep_mock_outliers$strep_630_metaT_reads, cex=1.9, pch=21, bg='red', col='black')



mtext('a', side=2, line=2, las=2, adj=1, padj=-12, cex=1.4, font=2)

#-------------------#

# Cefoperazone
plot(0, type='n', xlim=c(0,12), ylim=c(0,12), pch=20, xaxt='n', yaxt='n', xlab='', ylab='')
filledrectangle(wx=20, wy=2.8, col='gray80', mid=c(6,6), angle=45)
box()
points(x=cef_annotated$cef_mock_metaT_reads, y=cef_annotated$cef_630_metaT_reads, pch=20, cex=1.3, col='gray40')
segments(-2, -2, 14, 14, lty=2)
minor.ticks.axis(1, 12, mn=0, mx=12)
minor.ticks.axis(2, 12, mn=0, mx=12)
mtext('Fold Normalized cDNA Abundance', side=1, padj=2.8, cex=0.7)
mtext('Mock-Infected', side=1, padj=3.7, font=2, cex=0.9)
mtext('Fold Normalized cDNA Abundance', side=2, padj=-2.7, cex=0.7)
mtext(expression(bolditalic('C. difficile')~bold('630-Infected')), side=2, padj=-3.5, font=2, cex=0.9)
legend('topleft', c('Cefoperazone-pretreated', as.expression(bquote(paste(italic('rho'),' = ',.(cef_corr))))), bty='n', cex=1.2, text.col=c(wes_palette("FantasticFox")[3],'black'))


points(x=cef_630_outliers$cef_mock_metaT_reads, y=cef_630_outliers$cef_630_metaT_reads, cex=1.9, pch=21, bg='red', col='black')
points(x=cef_mock_outliers$cef_mock_metaT_reads, y=cef_mock_outliers$cef_630_metaT_reads, cex=1.9, pch=21, bg='red', col='black')

mtext('b', side=2, line=2, las=2, adj=1, padj=-12, cex=1.4, font=2)

#-------------------#

# Clindamycin
plot(0, type='n', xlim=c(0,12), ylim=c(0,12), pch=20, xaxt='n', yaxt='n', xlab='', ylab='')
filledrectangle(wx=20, wy=2.8, col='gray80', mid=c(6,6), angle=45)
box()
points(x=clinda_annotated$clinda_mock_metaT_reads, y=clinda_annotated$clinda_630_metaT_reads, pch=20, cex=1.3, col='gray40')
segments(-2, -2, 14, 14, lty=2)
minor.ticks.axis(1, 12, mn=0, mx=12)
minor.ticks.axis(2, 12, mn=0, mx=12)
mtext('Fold Normalized cDNA Abundance', side=1, padj=2.8, cex=0.7)
mtext('Mock-Infected', side=1, padj=3.7, font=2, cex=0.9)
mtext('Fold Normalized cDNA Abundance', side=2, padj=-2.7, cex=0.7)
mtext(expression(bolditalic('C. difficile')~bold('630-Infected')), side=2, padj=-3.5, font=2, cex=0.9)
legend('topleft', c('Clindamycin-pretreated', as.expression(bquote(paste(italic('rho'),' = ',.(clinda_corr))))), bty='n', cex=1.2, text.col=c(wes_palette("FantasticFox")[5],'black'))


points(x=clinda_630_outliers$clinda_mock_metaT_reads, y=clinda_630_outliers$clinda_630_metaT_reads, cex=1.9, pch=21, bg='red', col='black')
points(x=clinda_mock_outliers$clinda_mock_metaT_reads, y=clinda_mock_outliers$clinda_630_metaT_reads, cex=1.9, pch=21, bg='red', col='black')


mtext('c', side=2, line=2, las=2, adj=1, padj=-12, cex=1.4, font=2)

#-------------------#

# Taxonomic group legend
par(mar=c(0,0,0,0))
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-5,5), ylim=c(-10,10))

rect(xleft=-4.5, ybottom=10, xright=3.5, ytop=-5.5, border='black')


text(x=c(-3,-3,1,1,1,1), y=c(9,2.4,9,4.5,1.4,-1), labels=c('Bacteroidetes', 'Firmicutes','Actinobacteria','Proteobacteria','Verrucomicrobia','Other'), cex=1.2) # Phyla
text(x=-2.5, y=c(8.3,7.6,6.9,6.2,5.5,4.8,4.1), labels=c('Allistipes','Bacteroides','Odoribacter','Parabacteroides','Prevotella','Porphymonas','Other'), cex=0.9) # Bacteroietes
points(x=rep(-1.2,7), y=c(8.3,7.6,6.9,6.2,5.5,4.8,4.1), pch=22, cex=2.1, col='black', bg=c('#000066','#000099','#0000CC','#0000FF','#3333FF','#6666FF','#9999FF')) # blue

text(x=-2.5, y=c(1.7,1,0.3,-0.4,-1.1,-1.8,-2.5,-3.2,-3.9,-4.6), labels=c('Clostridium','Enterococcus','Eubacterium','Lactobacillus','Lactococcus','Roseburia','Ruminococcus','Staphylococcus','Streptococcus','Other'), cex=0.9) # Firmicutes
points(x=rep(-1.2,10), y=c(1.7,1,0.3,-0.4,-1.1,-1.8,-2.5,-3.2,-3.9,-4.6), pch=22, cex=2.1, col='black', bg=c('#330000','#660000','#990000','#CC0000','#FF0000','#FF3333','#FF6666','#FF9999','#FF99CC','#FFCCCC')) # red

text(x=1.5, y=c(8.3,7.6,6.9,6.2), labels=c('Bifidobacterium','Corynebacterium','Olsenella','Other'), cex=0.9) # Actinobacteria
points(x=rep(2.8,4), y=c(8.3,7.6,6.9,6.2), pch=22, cex=2.1, col='black', bg=c('#006600','#009900','#00CC00','#33FF33'))

text(x=1.5, y=c(3.8,3.1), labels=c('Escherichia','Other'), cex=0.9) # Proteobacteria
points(x=rep(2.8,2), y=c(3.8,3.1), pch=22, cex=2.1, col='black', bg=c('#CCCC00','#FFFF33'))

text(x=1.5, y=0.7, labels=c('Akkermansia'), cex=0.9) # Verrucomicrobia
points(x=2.8, y=0.7, pch=22, cex=2.1, col='black', bg='#990099')


points(x=2.8, y=-1, pch=22, cex=2.1, col='black', bg='#FF8000') # Other

dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
#for (dep in deps){
#  pkg <- paste('package:', dep, sep='')
#  detach(pkg, character.only = TRUE)}
#rm(list=ls())
#gc()
