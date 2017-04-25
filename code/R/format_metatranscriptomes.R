
# Set up environment
rm(list=ls())
gc()

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

rm(cef_metagenome, clinda_metagenome, strep_metagenome, 
   cef_630_metatranscriptome, cef_mock_metatranscriptome, clinda_630_metatranscriptome, 
   clinda_mock_metatranscriptome, strep_630_metatranscriptome, strep_mock_metatranscriptome)

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


# Add rarefaction curves with rarecurve()


# Rarefy read abundances
size <- round(min(colSums(cef_raw_reads[,c(1:3)]))*0.9) # Determine subsample level
cef_raw_reads$cef_metaG_reads <- t(rrarefy(cef_raw_reads$cef_metaG_reads, sample=size)) + 1
cef_raw_reads$cef_630_metaT_reads <- t(rrarefy(cef_raw_reads$cef_630_metaT_reads, sample=size)) + 1
cef_raw_reads$cef_mock_metaT_reads <- t(rrarefy(cef_raw_reads$cef_mock_metaT_reads, sample=size)) + 1
cef_normalized_reads <- cef_raw_reads
rm(cef_raw_reads)
size <- round(min(colSums(clinda_raw_reads[,c(1:3)]))*0.9) # Determine subsample level
clinda_raw_reads$clinda_metaG_reads <- t(rrarefy(clinda_raw_reads$clinda_metaG_reads, sample=size)) + 1
clinda_raw_reads$clinda_630_metaT_reads <- t(rrarefy(clinda_raw_reads$clinda_630_metaT_reads, sample=size)) + 1
clinda_raw_reads$clinda_mock_metaT_reads <- t(rrarefy(clinda_raw_reads$clinda_mock_metaT_reads, sample=size)) + 1
clinda_normalized_reads <- clinda_raw_reads
rm(clinda_raw_reads)
size <- round(min(colSums(strep_raw_reads[,c(1:3)]))*0.9) # Determine subsample level
strep_raw_reads$strep_metaG_reads <- t(rrarefy(strep_raw_reads$strep_metaG_reads, sample=size)) + 1
strep_raw_reads$strep_630_metaT_reads <- t(rrarefy(strep_raw_reads$strep_630_metaT_reads, sample=size)) + 1
strep_raw_reads$strep_mock_metaT_reads <- t(rrarefy(strep_raw_reads$strep_mock_metaT_reads, sample=size)) + 1
strep_normalized_reads <- strep_raw_reads
rm(strep_raw_reads)
size <- round(min(colSums(conv_raw_reads[,c(1:2)]))*0.9) # Determine subsample level
conv_raw_reads$conv_metaG_reads <- t(rrarefy(conv_raw_reads$conv_metaG_reads, sample=size)) + 1
conv_raw_reads$conv_metaT_reads <- t(rrarefy(conv_raw_reads$conv_metaT_reads, sample=size)) + 1
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
conv_normalized_reads <- subset(conv_normalized_reads, conv_normalized_reads$conv_metaT_reads > 0)

# Move row names to column
cef_normalized_reads$kegg_id <- rownames(cef_normalized_reads)
clinda_normalized_reads$kegg_id <- rownames(clinda_normalized_reads)
strep_normalized_reads$kegg_id <- rownames(strep_normalized_reads)
conv_normalized_reads$kegg_id <- rownames(conv_normalized_reads)

# Write normalized reads to files
write.table(cef_normalized_reads, file='~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/cef_normalized.tsv', sep='\t', row.names=FALSE, quote=FALSE)
write.table(clinda_normalized_reads, file='~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/clinda_normalized.tsv', sep='\t', row.names=FALSE, quote=FALSE)
write.table(strep_normalized_reads, file='~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/strep_normalized.tsv', sep='\t', row.names=FALSE, quote=FALSE)
write.table(conv_normalized_reads, file='~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/conv_normalized.tsv', sep='\t', row.names=FALSE, quote=FALSE)

# Clean up
rm(list=ls())
gc()

