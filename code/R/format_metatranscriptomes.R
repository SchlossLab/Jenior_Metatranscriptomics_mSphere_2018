
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')

# Define files
# Metagenomes
cef_630_metagenome <- 'data/read_mapping/metagenome/cefoperazone_630.Cefoperazone.metaG.final.pool.norm.txt'
cef_mock_metagenome <- 'data/read_mapping/metagenome/cefoperazone_mock.Cefoperazone.metaG.final.pool.norm.txt'
clinda_630_metagenome <- 'data/read_mapping/metagenome/clindamycin_630.Clindamycin.metaG.final.pool.norm.txt'
clinda_mock_metagenome <- 'data/read_mapping/metagenome/clindamycin_mock.Clindamycin.metaG.final.pool.norm.txt'
strep_630_metagenome <- 'data/read_mapping/metagenome/streptomycin_630.Streptomycin.metaG.final.pool.norm.txt'
strep_mock_metagenome <- 'data/read_mapping/metagenome/streptomycin_mock.Streptomycin.metaG.final.pool.norm.txt'
noabx_mock_metagenome <- 'data/read_mapping/metagenome/conventional_mock.Conventional.metaG.final.pool.norm.txt'
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
colnames(cef_630_metagenome) <- c('cef_630_metaG_reads')
cef_mock_metagenome <- read.delim(cef_mock_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(cef_mock_metagenome) <- c('cef_mock_metaG_reads')
clinda_630_metagenome <- read.delim(clinda_630_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(clinda_630_metagenome) <- c('clinda_630_metaG_reads')
clinda_mock_metagenome <- read.delim(clinda_mock_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(clinda_mock_metagenome) <- c('clinda_mock_metaG_reads')
strep_630_metagenome <- read.delim(strep_630_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(strep_630_metagenome) <- c('strep_630_metaG_reads')
strep_mock_metagenome <- read.delim(strep_mock_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(strep_mock_metagenome) <- c('strep_mock_metaG_reads')
noabx_mock_metagenome <- read.delim(noabx_mock_metagenome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(noabx_mock_metagenome) <- c('noabx_mock_metaG_reads')

# Metatranscriptomes
cef_630_metatranscriptome <- read.delim(cef_630_metatranscriptome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(cef_630_metatranscriptome) <- c('cef_630_metaT_reads')
cef_mock_metatranscriptome <- read.delim(cef_mock_metatranscriptome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(cef_mock_metatranscriptome) <- c('cef_mock_metaT_reads')
clinda_630_metatranscriptome <- read.delim(clinda_630_metatranscriptome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(clinda_630_metatranscriptome) <- c('clinda_630_metaT_reads')
clinda_mock_metatranscriptome <- read.delim(clinda_mock_metatranscriptome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(clinda_mock_metatranscriptome) <- c('clinda_mock_metaT_reads')
strep_630_metatranscriptome <- read.delim(strep_630_metatranscriptome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(strep_630_metatranscriptome) <- c('strep_630_metaT_reads')
strep_mock_metatranscriptome <- read.delim(strep_mock_metatranscriptome, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
colnames(strep_mock_metatranscriptome) <- c('strep_mock_metaT_reads')
noabx_mock_metatranscriptome <- read.delim(noabx_mock_metatranscriptome, sep='\t', header=FALSE, row.names=1)
colnames(noabx_mock_metatranscriptome) <- c('noabx_mock_metaT_reads')

#-------------------------------------------------------------------------------------------------------------------------#

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
   clinda_mock_metatranscriptome, strep_630_metatranscriptome, strep_mock_metatranscriptome,
   noabx_mock_metagenome, noabx_mock_metatranscriptome)

#-------------------------------------------------------------------------------------------------------------------------#

# Remove genes with no metagenomic coverage
cef_raw_reads <- subset(cef_raw_reads, cef_630_metaG_reads != 0 & cef_mock_metaG_reads != 0)
clinda_raw_reads <- subset(clinda_raw_reads, clinda_630_metaG_reads != 0 & clinda_mock_metaG_reads != 0)
strep_raw_reads <- subset(strep_raw_reads, strep_630_metaG_reads != 0 & strep_mock_metaG_reads != 0)
noabx_raw_reads <- subset(noabx_raw_reads, noabx_mock_metaG_reads != 0)

# Rarefy read abundances for metagenomes + metatranscriptomes
cef_size <- round(min(colSums(cef_raw_reads[,c(1:4)]))*0.95) # Determine subsample level
cef_raw_reads$cef_630_metaG_reads <- as.vector(rrarefy(cef_raw_reads$cef_630_metaG_reads, sample=cef_size)) + 1
cef_raw_reads$cef_mock_metaG_reads <- as.vector(rrarefy(cef_raw_reads$cef_mock_metaG_reads, sample=cef_size)) + 1
cef_raw_reads$cef_630_metaT_reads <- as.vector(rrarefy(cef_raw_reads$cef_630_metaT_reads, sample=cef_size)) + 1
cef_raw_reads$cef_mock_metaT_reads <- as.vector(rrarefy(cef_raw_reads$cef_mock_metaT_reads, sample=cef_size)) + 1
clinda_size <- round(min(colSums(clinda_raw_reads[,c(1:4)]))*0.95) # Determine subsample level
clinda_raw_reads$clinda_630_metaG_reads <- as.vector(rrarefy(clinda_raw_reads$clinda_630_metaG_reads, sample=clinda_size)) + 1
clinda_raw_reads$clinda_mock_metaG_reads <- as.vector(rrarefy(clinda_raw_reads$clinda_mock_metaG_reads, sample=clinda_size)) + 1
clinda_raw_reads$clinda_630_metaT_reads <- as.vector(rrarefy(clinda_raw_reads$clinda_630_metaT_reads, sample=clinda_size)) + 1
clinda_raw_reads$clinda_mock_metaT_reads <- as.vector(rrarefy(clinda_raw_reads$clinda_mock_metaT_reads, sample=clinda_size)) + 1
strep_size <- round(min(colSums(strep_raw_reads[,c(1:4)]))*0.95) # Determine subsample level
strep_raw_reads$strep_630_metaG_reads <- as.vector(rrarefy(strep_raw_reads$strep_630_metaG_reads, sample=strep_size)) + 1
strep_raw_reads$strep_mock_metaG_reads <- as.vector(rrarefy(strep_raw_reads$strep_mock_metaG_reads, sample=strep_size)) + 1
strep_raw_reads$strep_630_metaT_reads <- as.vector(rrarefy(strep_raw_reads$strep_630_metaT_reads, sample=strep_size)) + 1
strep_raw_reads$strep_mock_metaT_reads <- as.vector(rrarefy(strep_raw_reads$strep_mock_metaT_reads, sample=strep_size)) + 1
noabx_size <- round(min(colSums(noabx_raw_reads[,c(1:2)]))*0.95) # Determine subsample level
noabx_raw_reads$noabx_mock_metaG_reads <- as.vector(rrarefy(noabx_raw_reads$noabx_mock_metaG_reads, sample=noabx_size)) + 1
noabx_raw_reads$noabx_mock_metaT_reads <- as.vector(rrarefy(noabx_raw_reads$noabx_mock_metaT_reads, sample=noabx_size)) + 1

# Normalize metatranscriptomes to metagenomic coverage
cef_raw_reads$cef_630_metaT_reads <- cef_raw_reads$cef_630_metaT_reads / cef_raw_reads$cef_630_metaG_reads
cef_raw_reads$cef_mock_metaT_reads <- cef_raw_reads$cef_mock_metaT_reads / cef_raw_reads$cef_mock_metaG_reads
cef_raw_reads$cef_630_metaG_reads <- NULL
cef_raw_reads$cef_mock_metaG_reads <- NULL
clinda_raw_reads$clinda_630_metaT_reads <- clinda_raw_reads$clinda_630_metaT_reads / clinda_raw_reads$clinda_630_metaG_reads
clinda_raw_reads$clinda_mock_metaT_reads <- clinda_raw_reads$clinda_mock_metaT_reads / clinda_raw_reads$clinda_mock_metaG_reads
clinda_raw_reads$clinda_630_metaG_reads <- NULL
clinda_raw_reads$clinda_mock_metaG_reads <- NULL
strep_raw_reads$strep_630_metaT_reads <- strep_raw_reads$strep_630_metaT_reads / strep_raw_reads$strep_630_metaG_reads
strep_raw_reads$strep_mock_metaT_reads <- strep_raw_reads$strep_mock_metaT_reads / strep_raw_reads$strep_mock_metaG_reads
strep_raw_reads$strep_630_metaG_reads <- NULL
strep_raw_reads$strep_mock_metaG_reads <- NULL
noabx_raw_reads$noabx_mock_metaT_reads <- noabx_raw_reads$noabx_mock_metaT_reads / noabx_raw_reads$noabx_mock_metaG_reads
noabx_raw_reads$noabx_mock_metaG_reads <- NULL

# Add KEGG annotations
cef_kegg <- read.delim('~/Desktop/rows/cef_formatted.txt', sep='\t', header=TRUE, row.names=1)
clinda_kegg <- read.delim('~/Desktop/rows/clinda_formatted.txt', sep='\t', header=TRUE, row.names=1)
strep_kegg <- read.delim('~/Desktop/rows/strep_formatted.txt', sep='\t', header=TRUE, row.names=1)
noabx_kegg <- read.delim('~/Desktop/rows/noabx_formatted.txt', sep='\t', header=TRUE, row.names=1)
cef_raw_reads <- clean_merge(cef_raw_reads, cef_kegg)
clinda_raw_reads <- clean_merge(clinda_raw_reads, clinda_kegg)
strep_raw_reads <- clean_merge(strep_raw_reads, strep_kegg)
noabx_raw_reads <- clean_merge(noabx_raw_reads, noabx_kegg)
rm(cef_kegg,clinda_kegg,strep_kegg,noabx_kegg)

# Move row names to column
cef_raw_reads$kegg_hit <- rownames(cef_raw_reads)
clinda_raw_reads$kegg_hit <- rownames(clinda_raw_reads)
strep_raw_reads$kegg_hit <- rownames(strep_raw_reads)
noabx_raw_reads$kegg_hit <- rownames(noabx_raw_reads)

# Add KEGG pathways
pathways <- read.delim('~/Desktop/kegg_path/gene_paths.tsv', sep='\t', header=TRUE)
cef_raw_reads <- merge(cef_raw_reads, pathways, by='gene', all.x=TRUE, all.y=FALSE)
clinda_raw_reads <- merge(clinda_raw_reads, pathways, by='gene', all.x=TRUE, all.y=FALSE)
strep_raw_reads <- merge(strep_raw_reads, pathways, by='gene', all.x=TRUE, all.y=FALSE)
noabx_raw_reads <- merge(noabx_raw_reads, pathways, by='gene', all.x=TRUE, all.y=FALSE)
rm(pathways)

# Remove introduced duplicates
cef_raw_reads <- cef_raw_reads[!duplicated(cef_raw_reads$kegg_hit),]
clinda_raw_reads <- clinda_raw_reads[!duplicated(clinda_raw_reads$kegg_hit),]
strep_raw_reads <- strep_raw_reads[!duplicated(strep_raw_reads$kegg_hit),]
noabx_raw_reads <- noabx_raw_reads[!duplicated(noabx_raw_reads$kegg_hit),]

# Write raw reads to files
write.table(cef_raw_reads, file='data/read_mapping/cef_normalized_metaT.tsv', sep='\t', row.names=FALSE, quote=FALSE)
write.table(clinda_raw_reads, file='data/read_mapping/clinda_normalized_metaT.tsv', sep='\t', row.names=FALSE, quote=FALSE)
write.table(strep_raw_reads, file='data/read_mapping/strep_normalized_metaT.tsv', sep='\t', row.names=FALSE, quote=FALSE)
write.table(noabx_raw_reads, file='data/read_mapping/noabx_normalized_metaT.tsv', sep='\t', row.names=FALSE, quote=FALSE)

# Clean up
setwd(starting_dir)
rm(list=ls())
