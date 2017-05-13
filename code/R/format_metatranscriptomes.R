
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define files
# Metagenomes
cef_630_metagenome <- 'data/read_mapping/metagenome/cefoperazone_630.Cefoperazone.metaG.final.pool.norm.txt'
cef_mock_metagenome <- 'data/read_mapping/metagenome/cefoperazone_mock.Cefoperazone.metaG.final.pool.norm.txt'
clinda_630_metagenome <- 'data/read_mapping/metagenome/clindamycin_630.Clindamycin.metaG.final.pool.norm.txt'
clinda_mock_metagenome <- 'data/read_mapping/metagenome/clindamycin_mock.Clindamycin.metaG.final.pool.norm.txt'
strep_630_metagenome <- 'data/read_mapping/metagenome/streptomycin_630.Streptomycin.metaG.final.pool.norm.txt'
strep_mock_metagenome <- 'data/read_mapping/metagenome/streptomycin_mock.Streptomycin.metaG.final.pool.norm.txt'
noabx_mock_metagenome <- 'data/read_mapping/metagenome/conventional.Conventional.metaG.final.pool.norm.txt'
# Metatranscriptomes
cef_630_metatranscriptome <- 'data/read_mapping/metatranscriptome/cefoperazone_630.Cefoperazone.metaT.final.pool.norm.txt'
cef_mock_metatranscriptome <- 'data/read_mapping/metatranscriptome/cefoperazone_mock.Cefoperazone.metaT.final.pool.norm.txt'
clinda_630_metatranscriptome <- 'data/read_mapping/metatranscriptome/clindamycin_630.Clindamycin.metaT.final.pool.norm.txt'
clinda_mock_metatranscriptome <- 'data/read_mapping/metatranscriptome/clindamycin_mock.Clindamycin.metaT.final.pool.norm.txt'
strep_630_metatranscriptome <- 'data/read_mapping/metatranscriptome/streptomycin_630.Streptomycin.metaT.final.pool.norm.txt'
strep_mock_metatranscriptome <- 'data/read_mapping/metatranscriptome/streptomycin_mock.Streptomycin.metaT.final.pool.norm.txt'
noabx_mock_metatranscriptome <- 'data/read_mapping/metatranscriptome/conventional.Conventional.metaT.final.pool.norm.txt'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data

# Metagenomes
cef_630_metagenome <- read.delim(cef_630_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(cef_630_metagenome) <- c('cef_630_metaG_reads', 'ko', 'gene', 'pathway')
cef_mock_metagenome <- read.delim(cef_mock_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(cef_mock_metagenome) <- c('cef_mock_metaG_reads', 'ko', 'gene', 'pathway')
clinda_630_metagenome <- read.delim(clinda_630_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(clinda_630_metagenome) <- c('clinda_630_metaG_reads', 'ko', 'gene', 'pathway')
clinda_mock_metagenome <- read.delim(clinda_mock_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(clinda_mock_metagenome) <- c('clinda_mock_metaG_reads', 'ko', 'gene', 'pathway')
strep_630_metagenome <- read.delim(strep_630_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(strep_630_metagenome) <- c('strep_630_metaG_reads', 'ko', 'gene', 'pathway')
strep_mock_metagenome <- read.delim(strep_mock_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(strep_mock_metagenome) <- c('strep_mock_metaG_reads', 'ko', 'gene', 'pathway')
noabx_mock_metagenome <- read.delim(noabx_mock_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(noabx_mock_metagenome) <- c('noabx_mock_metaG_reads', 'ko', 'gene', 'pathway')

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
noabx_metatranscriptome <- read.delim(noabx_metatranscriptome, sep='\t', header=FALSE, row.names=1)
colnames(noabx_metatranscriptome) <- c('noabx_metaT_reads', 'ko', 'gene', 'pathway')

#-------------------------------------------------------------------------------------------------------------------------#

# Format metagenomic data for merging
cef_630_metagenome$ko <- NULL
cef_630_metagenome$gene <- NULL
cef_630_metagenome$pathway <- NULL
cef_mock_metagenome$ko <- NULL
cef_mock_metagenome$gene <- NULL
cef_mock_metagenome$pathway <- NULL
clinda_630_metagenome$ko <- NULL
clinda_630_metagenome$gene <- NULL
clinda_630_metagenome$pathway <- NULL
clinda_mock_metagenome$ko <- NULL
clinda_mock_metagenome$gene <- NULL
clinda_mock_metagenome$pathway <- NULL
strep_630_metagenome$ko <- NULL
strep_630_metagenome$gene <- NULL
strep_630_metagenome$pathway <- NULL
strep_mock_metagenome$ko <- NULL
strep_mock_metagenome$gene <- NULL
strep_mock_metagenome$pathway <- NULL
noabx_mock_metagenome$ko <- NULL
noabx_mock_metagenome$gene <- NULL
noabx_mock_metagenome$pathway <- NULL

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
cef_raw_reads <- clean_merge(cef_630_metagenome, cef_mock_metagenome)
cef_raw_reads <- clean_merge(cef_raw_reads, cef_630_metatranscriptome)
cef_raw_reads <- clean_merge(cef_raw_reads, cef_mock_metatranscriptome)
clinda_raw_reads <- clean_merge(clinda_630_metagenome, clinda_mock_metagenome)
clinda_raw_reads <- clean_merge(clinda_raw_reads, clinda_630_metatranscriptome)
clinda_raw_reads <- clean_merge(clinda_raw_reads, clinda_mock_metatranscriptome)
strep_raw_reads <- clean_merge(strep_630_metagenome, strep_mock_metagenome)
strep_raw_reads <- clean_merge(strep_raw_reads, strep_630_metatranscriptome)
strep_raw_reads <- clean_merge(strep_raw_reads, strep_mock_metatranscriptome)
noabx_raw_reads <- clean_merge(noabx_mock_metagenome, noabx_mock_metatranscriptome)

rm(cef_630_metagenome, clinda_630_metagenome, strep_630_metagenome, 
   cef_mock_metagenome, clinda_mock_metagenome, strep_mock_metagenome,
   cef_630_metatranscriptome, cef_mock_metatranscriptome, clinda_630_metatranscriptome, 
   clinda_mock_metatranscriptome, strep_630_metatranscriptome, strep_mock_metatranscriptome)

#-------------------------------------------------------------------------------------------------------------------------#

# Remove residual C. difficile mappings
#cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdf:CD630', rownames(cef_raw_reads)),]), ]
#cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdc:CD196', rownames(cef_raw_reads)),]), ]
#cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdf:CDP07', rownames(cef_raw_reads)),]), ]
#cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdf:CDP02', rownames(cef_raw_reads)),]), ]
#cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdf:CDP06', rownames(cef_raw_reads)),]), ]
#cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdf:CDP03', rownames(cef_raw_reads)),]), ]
#cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdg:CDBI1_08660', rownames(cef_raw_reads)),]), ]
#cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdg:CDBI1_06165', rownames(cef_raw_reads)),]), ]
#cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdg:CDBI1_09680', rownames(cef_raw_reads)),]), ]
#cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdg:CDBI1_02800', rownames(cef_raw_reads)),]), ]
#cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdg:CDBI1_13505', rownames(cef_raw_reads)),]), ]
#cef_raw_reads <- cef_raw_reads[!rownames(cef_raw_reads) %in% rownames(cef_raw_reads[grep('cdl:CDR20291_1755', rownames(cef_raw_reads)),]), ]
#clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdf:CD630', rownames(clinda_raw_reads)),]), ]
#clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdc:CD196', rownames(clinda_raw_reads)),]), ]
#clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdf:CDP07', rownames(clinda_raw_reads)),]), ]
#clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdf:CDP02', rownames(clinda_raw_reads)),]), ]
#clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdf:CDP06', rownames(clinda_raw_reads)),]), ]
#clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdf:CDP03', rownames(clinda_raw_reads)),]), ]
#clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdg:CDBI1_08660', rownames(clinda_raw_reads)),]), ]
#clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdg:CDBI1_06165', rownames(clinda_raw_reads)),]), ]
#clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdg:CDBI1_09680', rownames(clinda_raw_reads)),]), ]
#clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdg:CDBI1_02800', rownames(clinda_raw_reads)),]), ]
#clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdg:CDBI1_13505', rownames(clinda_raw_reads)),]), ]
#clinda_raw_reads <- clinda_raw_reads[!rownames(clinda_raw_reads) %in% rownames(clinda_raw_reads[grep('cdl:CDR20291_1755', rownames(clinda_raw_reads)),]), ]
#strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdf:CD630', rownames(strep_raw_reads)),]), ]
#strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdc:CD196', rownames(strep_raw_reads)),]), ]
#strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdf:CD630', rownames(strep_raw_reads)),]), ]
#strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdc:CD196', rownames(strep_raw_reads)),]), ]
#strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdf:CDP07', rownames(strep_raw_reads)),]), ]
#strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdf:CDP02', rownames(strep_raw_reads)),]), ]
#strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdf:CDP06', rownames(strep_raw_reads)),]), ]
#strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdf:CDP03', rownames(strep_raw_reads)),]), ]
#strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdg:CDBI1_08660', rownames(strep_raw_reads)),]), ]
#strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdg:CDBI1_06165', rownames(strep_raw_reads)),]), ]
#strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdg:CDBI1_09680', rownames(strep_raw_reads)),]), ]
#strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdg:CDBI1_02800', rownames(strep_raw_reads)),]), ]
#strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdg:CDBI1_13505', rownames(strep_raw_reads)),]), ]
#strep_raw_reads <- strep_raw_reads[!rownames(strep_raw_reads) %in% rownames(strep_raw_reads[grep('cdl:CDR20291_1755', rownames(strep_raw_reads)),]), ]
#noabx_raw_reads <- noabx_raw_reads[!rownames(noabx_raw_reads) %in% rownames(noabx_raw_reads[grep('cdf:CD630', rownames(noabx_raw_reads)),]), ]
#noabx_raw_reads <- noabx_raw_reads[!rownames(noabx_raw_reads) %in% rownames(noabx_raw_reads[grep('cdc:CD196', rownames(noabx_raw_reads)),]), ]
#noabx_raw_reads <- noabx_raw_reads[!rownames(noabx_raw_reads) %in% rownames(noabx_raw_reads[grep('cdf:CD630', rownames(noabx_raw_reads)),]), ]
#noabx_raw_reads <- noabx_raw_reads[!rownames(noabx_raw_reads) %in% rownames(noabx_raw_reads[grep('cdc:CD196', rownames(noabx_raw_reads)),]), ]
#noabx_raw_reads <- noabx_raw_reads[!rownames(noabx_raw_reads) %in% rownames(noabx_raw_reads[grep('cdf:CDP07', rownames(noabx_raw_reads)),]), ]
#noabx_raw_reads <- noabx_raw_reads[!rownames(noabx_raw_reads) %in% rownames(noabx_raw_reads[grep('cdf:CDP02', rownames(noabx_raw_reads)),]), ]
#noabx_raw_reads <- noabx_raw_reads[!rownames(noabx_raw_reads) %in% rownames(noabx_raw_reads[grep('cdf:CDP06', rownames(noabx_raw_reads)),]), ]
#noabx_raw_reads <- noabx_raw_reads[!rownames(noabx_raw_reads) %in% rownames(noabx_raw_reads[grep('cdf:CDP03', rownames(noabx_raw_reads)),]), ]
#noabx_raw_reads <- noabx_raw_reads[!rownames(noabx_raw_reads) %in% rownames(noabx_raw_reads[grep('cdg:CDBI1_08660', rownames(noabx_raw_reads)),]), ]
#noabx_raw_reads <- noabx_raw_reads[!rownames(noabx_raw_reads) %in% rownames(noabx_raw_reads[grep('cdg:CDBI1_06165', rownames(noabx_raw_reads)),]), ]
#noabx_raw_reads <- noabx_raw_reads[!rownames(noabx_raw_reads) %in% rownames(noabx_raw_reads[grep('cdg:CDBI1_09680', rownames(noabx_raw_reads)),]), ]
#noabx_raw_reads <- noabx_raw_reads[!rownames(noabx_raw_reads) %in% rownames(noabx_raw_reads[grep('cdg:CDBI1_02800', rownames(noabx_raw_reads)),]), ]
#noabx_raw_reads <- noabx_raw_reads[!rownames(noabx_raw_reads) %in% rownames(noabx_raw_reads[grep('cdg:CDBI1_13505', rownames(noabx_raw_reads)),]), ]
#noabx_raw_reads <- noabx_raw_reads[!rownames(noabx_raw_reads) %in% rownames(noabx_raw_reads[grep('cdl:CDR20291_1755', rownames(noabx_raw_reads)),]), ]

#-------------------------------------------------------------------------------------------------------------------------#

# Remove genes with no metagenomic coverage
cef_raw_reads <- subset(cef_raw_reads, cef_630_metaG_reads != 0 & cef_mock_metaG_reads != 0)
clinda_raw_reads <- subset(clinda_raw_reads, clinda_630_metaG_reads != 0 & clinda_mock_metaG_reads != 0)
strep_raw_reads <- subset(strep_raw_reads, strep_630_metaG_reads != 0 & strep_mock_metaG_reads != 0)
noabx_raw_reads <- subset(noabx_raw_reads, noabx_metaG_reads != 0)

# Rarefy read abundances for metagenomes + metatranscriptomes
cef_size <- round(min(colSums(cef_raw_reads[,c(1:3)]))*0.9) # Determine subsample level
cef_raw_reads$cef_630_metaG_reads <- t(rrarefy(cef_raw_reads$cef_630_metaG_reads, sample=cef_size)) + 1
cef_raw_reads$cef_mock_metaG_reads <- t(rrarefy(cef_raw_reads$cef_mock_metaG_reads, sample=cef_size)) + 1
cef_raw_reads$cef_630_metaT_reads <- t(rrarefy(cef_raw_reads$cef_630_metaT_reads, sample=cef_size)) + 1
cef_raw_reads$cef_mock_metaT_reads <- t(rrarefy(cef_raw_reads$cef_mock_metaT_reads, sample=cef_size)) + 1
cef_normalized_reads <- cef_raw_reads
rm(cef_raw_reads)
clinda_size <- round(min(colSums(clinda_raw_reads[,c(1:3)]))*0.9) # Determine subsample level
clinda_raw_reads$clinda_630_metaG_reads <- t(rrarefy(clinda_raw_reads$clinda_630_metaG_reads, sample=clinda_size)) + 1
clinda_raw_reads$clinda_mock_metaG_reads <- t(rrarefy(clinda_raw_reads$clinda_mock_metaG_reads, sample=clinda_size)) + 1
clinda_raw_reads$clinda_630_metaT_reads <- t(rrarefy(clinda_raw_reads$clinda_630_metaT_reads, sample=clinda_size)) + 1
clinda_raw_reads$clinda_mock_metaT_reads <- t(rrarefy(clinda_raw_reads$clinda_mock_metaT_reads, sample=clinda_size)) + 1
clinda_normalized_reads <- clinda_raw_reads
rm(clinda_raw_reads)
strep_size <- round(min(colSums(strep_raw_reads[,c(1:3)]))*0.9) # Determine subsample level
strep_raw_reads$strep_630_metaG_reads <- t(rrarefy(strep_raw_reads$strep_630_metaG_reads, sample=strep_size)) + 1
strep_raw_reads$strep_mock_metaG_reads <- t(rrarefy(strep_raw_reads$strep_mock_metaG_reads, sample=strep_size)) + 1
strep_raw_reads$strep_630_metaT_reads <- t(rrarefy(strep_raw_reads$strep_630_metaT_reads, sample=strep_size)) + 1
strep_raw_reads$strep_mock_metaT_reads <- t(rrarefy(strep_raw_reads$strep_mock_metaT_reads, sample=strep_size)) + 1
strep_normalized_reads <- strep_raw_reads
rm(strep_raw_reads)
noabx_size <- round(min(colSums(noabx_raw_reads[,c(1:2)]))*0.9) # Determine subsample level
noabx_raw_reads$noabx_metaG_reads <- t(rrarefy(noabx_raw_reads$noabx_metaG_reads, sample=noabx_size)) + 1
noabx_raw_reads$noabx_metaT_reads <- t(rrarefy(noabx_raw_reads$noabx_metaT_reads, sample=noabx_size)) + 1
noabx_normalized_reads <- noabx_raw_reads
rm(noabx_raw_reads)

# Normalize metatranscriptomes to metagenomic coverage
cef_normalized_reads$cef_630_metaT_reads <- cef_normalized_reads$cef_630_metaT_reads / cef_normalized_reads$cef_630_metaG_reads
cef_normalized_reads$cef_mock_metaT_reads <- cef_normalized_reads$cef_mock_metaT_reads / cef_normalized_reads$cef_mock_metaG_reads
cef_normalized_reads$cef_630_metaG_reads <- NULL
cef_normalized_reads$cef_mock_metaG_reads <- NULL
clinda_normalized_reads$clinda_630_metaT_reads <- clinda_normalized_reads$clinda_630_metaT_reads / clinda_normalized_reads$clinda_630_metaG_reads
clinda_normalized_reads$clinda_mock_metaT_reads <- clinda_normalized_reads$clinda_mock_metaT_reads / clinda_normalized_reads$clinda_mock_metaG_reads
clinda_normalized_reads$clinda_630_metaG_reads <- NULL
clinda_normalized_reads$clinda_mock_metaG_reads <- NULL
strep_normalized_reads$strep_630_metaT_reads <- strep_normalized_reads$strep_630_metaT_reads / strep_normalized_reads$strep_630_metaG_reads
strep_normalized_reads$strep_mock_metaT_reads <- strep_normalized_reads$strep_mock_metaT_reads / strep_normalized_reads$strep_mock_metaG_reads
strep_normalized_reads$strep_630_metaG_reads <- NULL
strep_normalized_reads$strep_mock_metaG_reads <- NULL
noabx_normalized_reads$noabx_metaT_reads <- noabx_normalized_reads$noabx_metaT_reads / noabx_normalized_reads$noabx_metaG_reads
noabx_normalized_reads$noabx_metaG_reads <- NULL

# Log2 transform the data
cef_normalized_reads[,c(1,2)] <- log2(cef_normalized_reads[,c(1,2)])
clinda_normalized_reads[,c(1,2)] <- log2(clinda_normalized_reads[,c(1,2)])
strep_normalized_reads[,c(1,2)] <- log2(strep_normalized_reads[,c(1,2)])
noabx_normalized_reads[,1] <- log2(noabx_normalized_reads[,1])

# Screen for active transcription in either condition
cef_normalized_reads <- subset(cef_normalized_reads, cef_normalized_reads$cef_630_metaT_reads > 0 | cef_normalized_reads$cef_mock_metaT_reads > 0)
clinda_normalized_reads <- subset(clinda_normalized_reads, clinda_normalized_reads$clinda_630_metaT_reads > 0 | clinda_normalized_reads$clinda_mock_metaT_reads > 0)
strep_normalized_reads <- subset(strep_normalized_reads, strep_normalized_reads$strep_630_metaT_reads > 0 | strep_normalized_reads$strep_mock_metaT_reads > 0)
noabx_normalized_reads <- subset(noabx_normalized_reads, noabx_normalized_reads$noabx_metaT_reads > 0)

# Move row names to column
cef_normalized_reads$kegg_id <- rownames(cef_normalized_reads)
clinda_normalized_reads$kegg_id <- rownames(clinda_normalized_reads)
strep_normalized_reads$kegg_id <- rownames(strep_normalized_reads)
noabx_normalized_reads$kegg_id <- rownames(noabx_normalized_reads)

# Write normalized reads to files
write.table(cef_normalized_reads, file='data/read_mapping/cef_normalized.tsv', sep='\t', row.names=FALSE, quote=FALSE)
write.table(clinda_normalized_reads, file='data/read_mapping/clinda_normalized.tsv', sep='\t', row.names=FALSE, quote=FALSE)
write.table(strep_normalized_reads, file='data/read_mapping/strep_normalized.tsv', sep='\t', row.names=FALSE, quote=FALSE)
write.table(noabx_normalized_reads, file='data/read_mapping/noabx_normalized.tsv', sep='\t', row.names=FALSE, quote=FALSE)

# Clean up
setwd(starting_dir)
rm(list=ls())
gc()
