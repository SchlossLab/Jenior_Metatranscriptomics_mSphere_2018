
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

# Process and reformat data

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

# Add KEGG annotations
cef_kegg <- read.delim('data/kegg/cef_formatted.txt', sep='\t', header=TRUE, row.names=1)
clinda_kegg <- read.delim('data/kegg/clinda_formatted.txt', sep='\t', header=TRUE, row.names=1)
strep_kegg <- read.delim('data/kegg/strep_formatted.txt', sep='\t', header=TRUE, row.names=1)
noabx_kegg <- read.delim('data/kegg/noabx_formatted.txt', sep='\t', header=TRUE, row.names=1)
cef_raw_reads <- clean_merge(cef_raw_reads, cef_kegg)
clinda_raw_reads <- clean_merge(clinda_raw_reads, clinda_kegg)
strep_raw_reads <- clean_merge(strep_raw_reads, strep_kegg)
noabx_raw_reads <- clean_merge(noabx_raw_reads, noabx_kegg)
rm(cef_kegg,clinda_kegg,strep_kegg,noabx_kegg)

# Remove introduced duplicates
cef_raw_reads$kegg_hit <- rownames(cef_raw_reads)
clinda_raw_reads$kegg_hit <- rownames(clinda_raw_reads)
strep_raw_reads$kegg_hit <- rownames(strep_raw_reads)
noabx_raw_reads$kegg_hit <- rownames(noabx_raw_reads)
cef_raw_reads <- cef_raw_reads[!duplicated(cef_raw_reads$kegg_hit),]
clinda_raw_reads <- clinda_raw_reads[!duplicated(clinda_raw_reads$kegg_hit),]
strep_raw_reads <- strep_raw_reads[!duplicated(strep_raw_reads$kegg_hit),]
noabx_raw_reads <- noabx_raw_reads[!duplicated(noabx_raw_reads$kegg_hit),]
cef_raw_reads$kegg_hit <- NULL
clinda_raw_reads$kegg_hit <- NULL
strep_raw_reads$kegg_hit <- NULL
noabx_raw_reads$kegg_hit <- NULL

# Add KEGG pathways
pathways <- read.delim('data/kegg/ko_paths.tsv', sep='\t', header=TRUE)
cef_raw_reads <- merge(cef_raw_reads, pathways, by='ko', all.x=TRUE)
clinda_raw_reads <- merge(clinda_raw_reads, pathways, by='ko', all.x=TRUE)
strep_raw_reads <- merge(strep_raw_reads, pathways, by='ko', all.x=TRUE)
noabx_raw_reads <- merge(noabx_raw_reads, pathways, by='ko', all.x=TRUE)
rm(pathways)

# Separate in groups and aggregate (Remove genes with no metagenomic coverage)
simp <- read.delim('data/kegg/simp_pathways.tsv', sep='\t', header=TRUE)
cef_mock_metaG <- as.data.frame(cbind(as.character(cef_raw_reads$cef_mock_metaG_reads), as.character(cef_raw_reads$pathways)), stringsAsFactors=FALSE)
colnames(cef_mock_metaG) <- c('abundance', 'pathway')
cef_mock_metaG$pathway[is.na(cef_mock_metaG$pathway)] <- 'Unannotated'
cef_mock_metaG$abundance <- as.numeric(cef_mock_metaG$abundance)
cef_mock_metaG <- merge(cef_mock_metaG, simp, by='pathway')
cef_630_metaG <- as.data.frame(cbind(as.character(cef_raw_reads$cef_630_metaG_reads), as.character(cef_raw_reads$pathways)), stringsAsFactors=FALSE)
colnames(cef_630_metaG) <- c('abundance', 'pathway')
cef_630_metaG$pathway[is.na(cef_630_metaG$pathway)] <- 'Unannotated'
cef_630_metaG$abundance <- as.numeric(cef_630_metaG$abundance)
cef_630_metaG <- merge(cef_630_metaG, simp, by='pathway')
clinda_mock_metaG <- as.data.frame(cbind(as.character(clinda_raw_reads$clinda_mock_metaG_reads), as.character(clinda_raw_reads$pathways)), stringsAsFactors=FALSE)
colnames(clinda_mock_metaG) <- c('abundance', 'pathway')
clinda_mock_metaG$pathway[is.na(clinda_mock_metaG$pathway)] <- 'Unannotated'
clinda_mock_metaG$abundance <- as.numeric(clinda_mock_metaG$abundance)
clinda_mock_metaG <- merge(clinda_mock_metaG, simp, by='pathway')
clinda_630_metaG <- as.data.frame(cbind(as.character(clinda_raw_reads$clinda_630_metaG_reads), as.character(clinda_raw_reads$pathways)), stringsAsFactors=FALSE)
colnames(clinda_630_metaG) <- c('abundance', 'pathway')
clinda_630_metaG$pathway[is.na(clinda_630_metaG$pathway)] <- 'Unannotated'
clinda_630_metaG$abundance <- as.numeric(clinda_630_metaG$abundance)
clinda_630_metaG <- merge(clinda_630_metaG, simp, by='pathway')
strep_mock_metaG <- as.data.frame(cbind(as.character(strep_raw_reads$strep_mock_metaG_reads), as.character(strep_raw_reads$pathways)), stringsAsFactors=FALSE)
colnames(strep_mock_metaG) <- c('abundance', 'pathway')
strep_mock_metaG$pathway[is.na(strep_mock_metaG$pathway)] <- 'Unannotated'
strep_mock_metaG$abundance <- as.numeric(strep_mock_metaG$abundance)
strep_mock_metaG <- merge(strep_mock_metaG, simp, by='pathway')
strep_630_metaG <- as.data.frame(cbind(as.character(strep_raw_reads$strep_630_metaG_reads), as.character(strep_raw_reads$pathways)), stringsAsFactors=FALSE)
colnames(strep_630_metaG) <- c('abundance', 'pathway')
strep_630_metaG$pathway[is.na(strep_630_metaG$pathway)] <- 'Unannotated'
strep_630_metaG$abundance <- as.numeric(strep_630_metaG$abundance)
strep_630_metaG <- merge(strep_630_metaG, simp, by='pathway')
noabx_mock_metaG <- as.data.frame(cbind(as.character(noabx_raw_reads$noabx_mock_metaG_reads), as.character(noabx_raw_reads$pathways)), stringsAsFactors=FALSE)
colnames(noabx_mock_metaG) <- c('abundance', 'pathway')
noabx_mock_metaG$pathway[is.na(noabx_mock_metaG$pathway)] <- 'Unannotated'
noabx_mock_metaG$abundance <- as.numeric(noabx_mock_metaG$abundance)
noabx_mock_metaG <- merge(noabx_mock_metaG, simp, by='pathway')
cef_mock_metaT <- as.data.frame(cbind(as.character(cef_raw_reads$cef_mock_metaT_reads), as.character(cef_raw_reads$pathways)), stringsAsFactors=FALSE)
colnames(cef_mock_metaT) <- c('abundance', 'pathway')
cef_mock_metaT$pathway[is.na(cef_mock_metaT$pathway)] <- 'Unannotated'
cef_mock_metaT$abundance <- as.numeric(cef_mock_metaT$abundance)
cef_mock_metaT <- merge(cef_mock_metaT, simp, by='pathway')
cef_630_metaT <- as.data.frame(cbind(as.character(cef_raw_reads$cef_630_metaT_reads), as.character(cef_raw_reads$pathways)), stringsAsFactors=FALSE)
colnames(cef_630_metaT) <- c('abundance', 'pathway')
cef_630_metaT$pathway[is.na(cef_630_metaT$pathway)] <- 'Unannotated'
cef_630_metaT$abundance <- as.numeric(cef_630_metaT$abundance)
cef_630_metaT <- merge(cef_630_metaT, simp, by='pathway')
clinda_mock_metaT <- as.data.frame(cbind(as.character(clinda_raw_reads$clinda_mock_metaT_reads), as.character(clinda_raw_reads$pathways)), stringsAsFactors=FALSE)
colnames(clinda_mock_metaT) <- c('abundance', 'pathway')
clinda_mock_metaT$pathway[is.na(clinda_mock_metaT$pathway)] <- 'Unannotated'
clinda_mock_metaT$abundance <- as.numeric(clinda_mock_metaT$abundance)
clinda_mock_metaT <- merge(clinda_mock_metaT, simp, by='pathway')
clinda_630_metaT <- as.data.frame(cbind(as.character(clinda_raw_reads$clinda_630_metaT_reads), as.character(clinda_raw_reads$pathways)), stringsAsFactors=FALSE)
colnames(clinda_630_metaT) <- c('abundance', 'pathway')
clinda_630_metaT$pathway[is.na(clinda_630_metaT$pathway)] <- 'Unannotated'
clinda_630_metaT$abundance <- as.numeric(clinda_630_metaT$abundance)
clinda_630_metaT <- merge(clinda_630_metaT, simp, by='pathway')
strep_mock_metaT <- as.data.frame(cbind(as.character(strep_raw_reads$strep_mock_metaT_reads), as.character(strep_raw_reads$pathways)), stringsAsFactors=FALSE)
colnames(strep_mock_metaT) <- c('abundance', 'pathway')
strep_mock_metaT$pathway[is.na(strep_mock_metaT$pathway)] <- 'Unannotated'
strep_mock_metaT$abundance <- as.numeric(strep_mock_metaT$abundance)
strep_mock_metaT <- merge(strep_mock_metaT, simp, by='pathway')
strep_630_metaT <- as.data.frame(cbind(as.character(strep_raw_reads$strep_630_metaT_reads), as.character(strep_raw_reads$pathways)), stringsAsFactors=FALSE)
colnames(strep_630_metaT) <- c('abundance', 'pathway')
strep_630_metaT$pathway[is.na(strep_630_metaT$pathway)] <- 'Unannotated'
strep_630_metaT$abundance <- as.numeric(strep_630_metaT$abundance)
strep_630_metaT <- merge(strep_630_metaT, simp, by='pathway')
noabx_mock_metaT <- as.data.frame(cbind(as.character(noabx_raw_reads$noabx_mock_metaT_reads), as.character(noabx_raw_reads$pathways)), stringsAsFactors=FALSE)
colnames(noabx_mock_metaT) <- c('abundance', 'pathway')
noabx_mock_metaT$pathway[is.na(noabx_mock_metaT$pathway)] <- 'Unannotated'
noabx_mock_metaT$abundance <- as.numeric(noabx_mock_metaT$abundance)
noabx_mock_metaT <- merge(noabx_mock_metaT, simp, by='pathway')
rm(simp, cef_raw_reads, clinda_raw_reads, strep_raw_reads, noabx_raw_reads)

# Remove groups with no reads
cef_630_metaG <- subset(cef_630_metaG, abundance > 0)
cef_mock_metaG <- subset(cef_mock_metaG, abundance > 0)
cef_630_metaT <- subset(cef_630_metaT, abundance > 0)
cef_mock_metaT <- subset(cef_mock_metaT, abundance > 0)
strep_630_metaG <- subset(strep_630_metaG, abundance > 0)
strep_mock_metaG <- subset(strep_mock_metaG, abundance > 0)
strep_630_metaT <- subset(strep_630_metaT, abundance > 0)
strep_mock_metaT <- subset(strep_mock_metaT, abundance > 0)
clinda_630_metaG <- subset(clinda_630_metaG, abundance > 0)
clinda_mock_metaG <- subset(clinda_mock_metaG, abundance > 0)
clinda_630_metaT <- subset(clinda_630_metaT, abundance > 0)
clinda_mock_metaT <- subset(clinda_mock_metaT, abundance > 0)
noabx_mock_metaG <- subset(noabx_mock_metaG, abundance > 0)
noabx_mock_metaT <- subset(noabx_mock_metaT, abundance > 0)

# Tabulate results
cef_mock_metaG <- as.data.frame(table(cef_mock_metaG$simplified))
cef_mock_metaG <- subset(cef_mock_metaG, Freq != 0)
cef_630_metaG <- as.data.frame(table(cef_630_metaG$simplified))
cef_630_metaG <- subset(cef_630_metaG, Freq != 0)
clinda_mock_metaG <- as.data.frame(table(clinda_mock_metaG$simplified))
clinda_mock_metaG <- subset(clinda_mock_metaG, Freq != 0)
clinda_630_metaG <- as.data.frame(table(clinda_630_metaG$simplified))
clinda_630_metaG <- subset(clinda_630_metaG, Freq != 0)
strep_mock_metaG <- as.data.frame(table(strep_mock_metaG$simplified))
strep_mock_metaG <- subset(strep_mock_metaG, Freq != 0)
strep_630_metaG <- as.data.frame(table(strep_630_metaG$simplified))
strep_630_metaG <- subset(strep_630_metaG, Freq != 0)
noabx_mock_metaG <- as.data.frame(table(noabx_mock_metaG$simplified))
noabx_mock_metaG <- subset(noabx_mock_metaG, Freq != 0)
cef_mock_metaT <- as.data.frame(table(cef_mock_metaT$simplified))
cef_mock_metaT <- subset(cef_mock_metaT, Freq != 0)
cef_630_metaT <- as.data.frame(table(cef_630_metaT$simplified))
cef_630_metaT <- subset(cef_630_metaT, Freq != 0)
clinda_mock_metaT <- as.data.frame(table(clinda_mock_metaT$simplified))
clinda_mock_metaT <- subset(clinda_mock_metaT, Freq != 0)
clinda_630_metaT <- as.data.frame(table(clinda_630_metaT$simplified))
clinda_630_metaT <- subset(clinda_630_metaT, Freq != 0)
strep_mock_metaT <- as.data.frame(table(strep_mock_metaT$simplified))
strep_mock_metaT <- subset(strep_mock_metaT, Freq != 0)
strep_630_metaT <- as.data.frame(table(strep_630_metaT$simplified))
strep_630_metaT <- subset(strep_630_metaT, Freq != 0)
noabx_mock_metaT <- as.data.frame(table(noabx_mock_metaT$simplified))
noabx_mock_metaT <- subset(noabx_mock_metaT, Freq != 0)

# Combine sample types
metaG <- merge(strep_mock_metaG, cef_mock_metaG, by='Var1', all=TRUE)
colnames(metaG) <- c('Var1','strep_mock','cef_mock')
metaG <- merge(metaG, clinda_mock_metaG, by='Var1', all=TRUE)
colnames(metaG) <- c('Var1','strep_mock','cef_mock','clinda_mock')
metaG <- merge(metaG, noabx_mock_metaG, by='Var1', all=TRUE)
colnames(metaG) <- c('pathway','strep_mock','cef_mock','clinda_mock','noabx_mock')
metaG[is.na(metaG)] <- 0
metaG$pathway <- gsub('_', ' ', metaG$pathway)
rownames(metaG) <- metaG$pathway
metaG$pathway <- NULL
metaT <- merge(strep_630_metaT, strep_mock_metaT, by='Var1', all=TRUE)
colnames(metaT) <- c('Var1','strep_630','strep_mock')
metaT <- merge(metaT, cef_630_metaT, by='Var1', all=TRUE)
colnames(metaT) <- c('Var1','strep_630','strep_mock','cef_630')
metaT <- merge(metaT, cef_mock_metaT, by='Var1', all=TRUE)
colnames(metaT) <- c('Var1','strep_630','strep_mock','cef_630','cef_mock')
metaT <- merge(metaT, clinda_630_metaT, by='Var1', all=TRUE)
colnames(metaT) <- c('Var1','strep_630','strep_mock','cef_630','cef_mock','clinda_630')
metaT <- merge(metaT, clinda_mock_metaT, by='Var1', all=TRUE)
colnames(metaT) <- c('Var1','strep_630','strep_mock','cef_630','cef_mock','clinda_630','clinda_mock')
metaT <- merge(metaT, noabx_mock_metaT, by='Var1', all=TRUE)
colnames(metaT) <- c('pathway','strep_630','strep_mock','cef_630','cef_mock','clinda_630','clinda_mock','noabx_mock')
metaT[is.na(metaT)] <- 0
metaT$pathway <- gsub('_', ' ', metaT$pathway)
rownames(metaT) <- metaT$pathway
metaT$pathway <- NULL
rm(strep_630_metaG, strep_mock_metaG, cef_630_metaG, cef_mock_metaG, clinda_630_metaG, clinda_mock_metaG, noabx_mock_metaG,
   strep_630_metaT, strep_mock_metaT, cef_630_metaT, cef_mock_metaT, clinda_630_metaT, clinda_mock_metaT, noabx_mock_metaT)

# Remove Unannotated catagory and reassign for just those groups with annotations
metaG <- subset(metaG, rownames(metaG) != 'Unannotated')
metaG_other <- subset(metaG, rowSums(metaG) <= 5)
metaG <- subset(metaG, rowSums(metaG) > 5)
metaG_other <- colSums(metaG) - colSums(metaG_other)
metaG_rows <- c(rownames(metaG), 'Other') 
metaG <- as.data.frame(rbind(metaG, metaG_other))
rownames(metaG) <- metaG_rows
rm(metaG_rows, metaG_other)
metaT <- subset(metaT, rownames(metaT) != 'Unannotated')
metaT_other <- subset(metaT, rowSums(metaT) <= 5)
metaT <- subset(metaT, rowSums(metaT) > 5)
metaT_other <- colSums(metaT) - colSums(metaT_other)
metaT_rows <- c(rownames(metaT), 'Other') 
metaT <- as.data.frame(rbind(metaT, metaT_other))
rownames(metaT) <- metaT_rows
rm(metaT_rows, metaT_other)

# Break into treatment groups
cef <- as.data.frame((cbind(metaG$cef_mock, metaT$cef_630, metaT$cef_mock)))
rownames(cef) <- 



metaG <- as.matrix(metaG)

metaT <- as.matrix(metaT)

#-------------------------------------------------------------------------------------------------------------------------#

# Generate plot
bar_palette <- viridis(n=max(c(nrow(metaG), nrow(metaT))))
#layout(matrix(c(1,3,
#                2,3), 
#              nrow=2, ncol=2, byrow=FALSE))
#par(mar=c(5, 5, 1, 1), mgp=c(3,0.7,0))



# Stacked Bar Plot with Colors and Legend
barplot(metaG, main='', xlab='', col=bar_palette) 



barplot(metaT, main='', xlab='', col=bar_palette) 









# Clean up
#rm(list=ls())
#gc()
