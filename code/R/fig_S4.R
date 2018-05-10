
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_mSphere_2018/code/R/functions.R')

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

# Output file
plot_file <- 'results/supplement/figures/figure_S4.pdf'

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
rm(cef_kegg, clinda_kegg, strep_kegg, noabx_kegg)

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

# Screen for annotated genes
cef_annotated <- subset(cef_raw_reads, !is.na(organism))
cef_annotated <- subset(cef_annotated, !description %in% c('unknown_function', 'unlikely', 'hypothetical_protein'))
cef_annotated <- subset(cef_annotated, !grepl("uncharacterized_", cef_annotated$description))
clinda_annotated <- subset(clinda_raw_reads, !is.na(organism))
clinda_annotated <- subset(clinda_annotated, !description %in% c('unknown_function', 'unlikely', 'hypothetical_protein'))
clinda_annotated <- subset(clinda_annotated, !grepl("uncharacterized_", clinda_annotated$description))
strep_annotated <- subset(strep_raw_reads, !is.na(organism))
strep_annotated <- subset(strep_annotated, !description %in% c('unknown_function', 'unlikely', 'hypothetical_protein'))
strep_annotated <- subset(strep_annotated, !grepl("uncharacterized_", strep_annotated$description))
noabx_annotated <- subset(noabx_raw_reads, !is.na(organism))
noabx_annotated <- subset(noabx_annotated, !description %in% c('unknown_function', 'unlikely', 'hypothetical_protein'))
noabx_annotated <- subset(noabx_annotated, !grepl("uncharacterized_", noabx_annotated$description))
cef_raw_reads <- cef_annotated
clinda_raw_reads <- clinda_annotated
strep_raw_reads <- strep_annotated
noabx_raw_reads <- noabx_annotated
rm(cef_annotated, clinda_annotated, strep_annotated, noabx_annotated)

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

# Combine treatment types
cef <- merge(cef_mock_metaG, cef_mock_metaT, by='Var1', all=TRUE)
colnames(cef) <- c('Var1','cef_mock_metaG','cef_mock_metaT')
cef <- merge(cef, cef_630_metaT, by='Var1', all=TRUE)
colnames(cef) <- c('pathway','Metagenome','Mock','CDI')
cef[is.na(cef)] <- 0
cef$pathway <- gsub('_', ' ', cef$pathway)
rownames(cef) <- cef$pathway
cef$pathway <- NULL
strep <- merge(strep_mock_metaG, strep_mock_metaT, by='Var1', all=TRUE)
colnames(strep) <- c('Var1','strep_mock_metaG','strep_mock_metaT')
strep <- merge(strep, strep_630_metaT, by='Var1', all=TRUE)
colnames(strep) <- c('pathway','Metagenome','Mock','CDI')
strep[is.na(strep)] <- 0
strep$pathway <- gsub('_', ' ', strep$pathway)
rownames(strep) <- strep$pathway
strep$pathway <- NULL
clinda <- merge(clinda_mock_metaG, clinda_mock_metaT, by='Var1', all=TRUE)
colnames(clinda) <- c('Var1','clinda_mock_metaG','clinda_mock_metaT')
clinda <- merge(clinda, clinda_630_metaT, by='Var1', all=TRUE)
colnames(clinda) <- c('pathway','Metagenome','Mock','CDI')
clinda[is.na(clinda)] <- 0
clinda$pathway <- gsub('_', ' ', clinda$pathway)
rownames(clinda) <- clinda$pathway
clinda$pathway <- NULL
noabx <- merge(noabx_mock_metaG, noabx_mock_metaT, by='Var1', all=TRUE)
colnames(noabx) <- c('pathway','Metagenome','Metatranscriptome')
noabx[is.na(noabx)] <- 0
noabx$pathway <- gsub('_', ' ', noabx$pathway)
rownames(noabx) <- noabx$pathway
noabx$pathway <- NULL
rm(strep_630_metaG, strep_mock_metaG, cef_630_metaG, cef_mock_metaG, clinda_630_metaG, clinda_mock_metaG, noabx_mock_metaG,
   strep_630_metaT, strep_mock_metaT, cef_630_metaT, cef_mock_metaT, clinda_630_metaT, clinda_mock_metaT, noabx_mock_metaT)

# Bin low abundance in Other
cef_sums <- colSums(cef[,1:3])
cef <- cef[rowSums(cef) > 15,]
cef_sums <- cef_sums - colSums(cef) 
cef <- as.data.frame(t(cef))
cef$Other <- as.numeric(cef_sums)
cef <- t(cef)
strep_sums <- colSums(strep[,1:3])
strep <- strep[rowSums(strep) > 15,]
strep_sums <- strep_sums - colSums(strep) 
strep <- as.data.frame(t(strep))
strep$Other <- as.numeric(strep_sums)
strep <- t(strep)
clinda_sums <- colSums(clinda[,1:3])
clinda <- clinda[rowSums(clinda) > 15,]
clinda_sums <- clinda_sums - colSums(clinda) 
clinda <- as.data.frame(t(clinda))
clinda$Other <- as.numeric(clinda_sums)
clinda <- t(clinda)
noabx_sums <- colSums(noabx[,1:2])
noabx <- noabx[rowSums(noabx) > 15,]
noabx_sums <- noabx_sums - colSums(noabx) 
noabx <- as.data.frame(t(noabx))
noabx$Other <- as.numeric(noabx_sums)
noabx <- t(noabx)

# Assign colors 
bar_palette <- viridis(n=nrow(noabx))
bar_palette <- as.data.frame(cbind(rownames(noabx), bar_palette))
colnames(bar_palette) <- c('pathway', 'color')
cef <- merge(cef, bar_palette, by.x='row.names', by.y='pathway')
rownames(cef) <- cef$Row.names
cef$Row.names <- NULL
strep <- merge(strep, bar_palette, by.x='row.names', by.y='pathway')
rownames(strep) <- strep$Row.names
strep$Row.names <- NULL
clinda <- merge(clinda, bar_palette, by.x='row.names', by.y='pathway')
rownames(clinda) <- clinda$Row.names
clinda$Row.names <- NULL
noabx <- merge(noabx, bar_palette, by.x='row.names', by.y='pathway')
rownames(noabx) <- noabx$Row.names
noabx$Row.names <- NULL

# Reorder for largest groups on top
cef <- cef[order(cef$Metagenome),] 
strep <- strep[order(strep$Metagenome),] 
clinda <- clinda[order(clinda$Metagenome),] 
noabx <- noabx[order(noabx$Metagenome),] 


#-------------------------------------------------------------------------------------------------------------------------#

# Generate plot
pdf(file=plot_file, width=6, height=8)
layout(matrix(c(1,2,
                3,4,
                5,5), 
              nrow=3, ncol=2, byrow=TRUE))

par(mar=c(3, 4, 2, 1), mgp=c(3,0.7,0), las=1)
barplot(as.matrix(noabx[,1:2]), main='No Antibiotics', xlab='', ylab='Unique Genes', col=as.character(noabx$color), ylim=c(0,10000))
box(lwd=1.5)
abline(v=1.3, lty=2, lwd=1.5)
mtext('A', side=2, line=2, las=2, adj=1.5, padj=-8, cex=1.2, font=2)

barplot(as.matrix(strep[,1:3]), main='Streptomycin-pretreated', xlab='', ylab='Unique Genes', col=as.character(strep$color), ylim=c(0,3000))
box(lwd=1.5)
mtext('Metatranscriptome', side=1, line=2, las=2, adj=0.75, padj=-0.5, cex=0.7, las=1)
abline(v=1.3, lty=2, lwd=1.5)
mtext('B', side=2, line=2, las=2, adj=1.5, padj=-8, cex=1.2, font=2)

barplot(as.matrix(cef[,1:3]), main='Cefoperazone-pretreated', xlab='', ylab='Unique Genes', col=as.character(cef$color), ylim=c(0,2000))
box(lwd=1.5)
mtext('Metatranscriptome', side=1, line=2, las=2, adj=0.75, padj=-0.5, cex=0.7, las=1)
abline(v=1.3, lty=2, lwd=1.5)
mtext('C', side=2, line=2, las=2, adj=1.5, padj=-8, cex=1.2, font=2)

barplot(as.matrix(clinda[,1:3]), main='Clindamycin-pretreated', xlab='', ylab='Unique Genes', col=as.character(clinda$color), ylim=c(0,2000))
box(lwd=1.5)
mtext('Metatranscriptome', side=1, line=2, las=2, adj=0.75, padj=-0.5, cex=0.7, las=1)
abline(v=1.3, lty=2, lwd=1.5)
mtext('D', side=2, line=2, las=2, adj=1.5, padj=-8, cex=1.2, font=2)

par(mar=c(0,0,1,0))
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-6,6))
legend('top', legend=as.character(bar_palette$pathway), bty='n', 
       pt.bg=rev(as.character(bar_palette$color)), pch=22, cex=1.1, pt.cex=1.9, ncol=2)

dev.off()

#----------------#

# Clean up
rm(list=ls())
gc()

#-------------------------------------------------------------------------------------------------------------------------#

# Panels E and F

# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')

# Define files
# MetaGs
cef_630_metaG <- 'data/read_mapping/metagenome/cefoperazone_630.Cefoperazone.metaG.final.pool.norm.txt'
cef_mock_metaG <- 'data/read_mapping/metagenome/cefoperazone_mock.Cefoperazone.metaG.final.pool.norm.txt'
clinda_630_metaG <- 'data/read_mapping/metagenome/clindamycin_630.Clindamycin.metaG.final.pool.norm.txt'
clinda_mock_metaG <- 'data/read_mapping/metagenome/clindamycin_mock.Clindamycin.metaG.final.pool.norm.txt'
strep_630_metaG <- 'data/read_mapping/metagenome/streptomycin_630.Streptomycin.metaG.final.pool.norm.txt'
strep_mock_metaG <- 'data/read_mapping/metagenome/streptomycin_mock.Streptomycin.metaG.final.pool.norm.txt'
noabx_mock_metaG <- 'data/read_mapping/metagenome/conventional_mock.Conventional.metaG.final.pool.norm.txt'
# MetaTs
cef_630_metaT <- 'data/read_mapping/metatranscriptome/cefoperazone_630.Cefoperazone.metaT.final.pool.norm.txt'
cef_mock_metaT <- 'data/read_mapping/metatranscriptome/cefoperazone_mock.Cefoperazone.metaT.final.pool.norm.txt'
clinda_630_metaT <- 'data/read_mapping/metatranscriptome/clindamycin_630.Clindamycin.metaT.final.pool.norm.txt'
clinda_mock_metaT <- 'data/read_mapping/metatranscriptome/clindamycin_mock.Clindamycin.metaT.final.pool.norm.txt'
strep_630_metaT <- 'data/read_mapping/metatranscriptome/streptomycin_630.Streptomycin.metaT.final.pool.norm.txt'
strep_mock_metaT <- 'data/read_mapping/metatranscriptome/streptomycin_mock.Streptomycin.metaT.final.pool.norm.txt'
noabx_mock_metaT <- 'data/read_mapping/metatranscriptome/conventional.Conventional.metaT.final.pool.norm.txt'

plot_file <- 'results/supplement/figures/figure_S5EF.pdf'


#--------------------------------------------------------------------------------------------------#

# Read in data
# MetaGs
cef_630_metaG <- read.delim(cef_630_metaG, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
cef_630_metaG <- cef_630_metaG[cef_630_metaG != 0]
cef_mock_metaG <- read.delim(cef_mock_metaG, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
cef_mock_metaG <- cef_mock_metaG[cef_mock_metaG != 0]
clinda_630_metaG <- read.delim(clinda_630_metaG, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
clinda_630_metaG <- clinda_630_metaG[clinda_630_metaG != 0]
clinda_mock_metaG <- read.delim(clinda_mock_metaG, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
clinda_mock_metaG <- clinda_mock_metaG[clinda_mock_metaG != 0]
strep_630_metaG <- read.delim(strep_630_metaG, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
strep_630_metaG <- strep_630_metaG[strep_630_metaG != 0]
strep_mock_metaG <- read.delim(strep_mock_metaG, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
strep_mock_metaG <- strep_mock_metaG[strep_mock_metaG != 0]
noabx_mock_metaG <- read.delim(noabx_mock_metaG, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
noabx_mock_metaG <- noabx_mock_metaG[noabx_mock_metaG != 0]

# MetaTs
cef_630_metaT <- read.delim(cef_630_metaT, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
cef_630_metaT <- cef_630_metaT[cef_630_metaT != 0]
cef_mock_metaT <- read.delim(cef_mock_metaT, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
cef_mock_metaT <- cef_mock_metaT[cef_mock_metaT != 0]
clinda_630_metaT <- read.delim(clinda_630_metaT, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
clinda_630_metaT <- clinda_630_metaT[clinda_630_metaT != 0]
clinda_mock_metaT <- read.delim(clinda_mock_metaT, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
clinda_mock_metaT <- clinda_mock_metaT[clinda_mock_metaT != 0]
strep_630_metaT <- read.delim(strep_630_metaT, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
strep_630_metaT <- strep_630_metaT[strep_630_metaT != 0]
strep_mock_metaT <- read.delim(strep_mock_metaT, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
strep_mock_metaT <- strep_mock_metaT[strep_mock_metaT != 0]
noabx_mock_metaT <- read.delim(noabx_mock_metaT, sep='\t', header=FALSE, row.names=1)[,1]
noabx_mock_metaT <- noabx_mock_metaT[noabx_mock_metaT != 0]

#--------------------------------------------------------------------------------------------------#

# Perform rarefaction
cef_630_metaG <- stepRarefy(cef_630_metaG)
cef_mock_metaG <- stepRarefy(cef_mock_metaG)
strep_630_metaG <- stepRarefy(strep_630_metaG)
strep_mock_metaG <- stepRarefy(strep_mock_metaG)
clinda_630_metaG <- stepRarefy(clinda_630_metaG)
clinda_mock_metaG <- stepRarefy(clinda_mock_metaG)
noabx_mock_metaG <- stepRarefy(noabx_mock_metaG)
cef_630_metaT <- stepRarefy(cef_630_metaT)
cef_mock_metaT <- stepRarefy(cef_mock_metaT)
strep_630_metaT <- stepRarefy(strep_630_metaT)
strep_mock_metaT <- stepRarefy(strep_mock_metaT)
clinda_630_metaT <- stepRarefy(clinda_630_metaT)
clinda_mock_metaT <- stepRarefy(clinda_mock_metaT)
noabx_mock_metaT <- stepRarefy(noabx_mock_metaT)

# Log transform
cef_630_metaG <- log10(cef_630_metaG)
cef_mock_metaG <- log10(cef_mock_metaG)
strep_630_metaG <- log10(strep_630_metaG)
strep_mock_metaG <- log10(strep_mock_metaG)
clinda_630_metaG <- log10(clinda_630_metaG)
clinda_mock_metaG <- log10(clinda_mock_metaG)
noabx_mock_metaG <- log10(noabx_mock_metaG)
cef_630_metaT <- log10(cef_630_metaT)
cef_mock_metaT <- log10(cef_mock_metaT)
strep_630_metaT <- log10(strep_630_metaT)
strep_mock_metaT <- log10(strep_mock_metaT)
clinda_630_metaT <- log10(clinda_630_metaT)
clinda_mock_metaT <- log10(clinda_mock_metaT)
noabx_mock_metaT <- log10(noabx_mock_metaT)

#--------------------------------------------------------------------------------------------------#

# Plot collector's curves
pdf(file=plot_file, width=10, height=5)
layout(matrix(c(1,2),
              nrow=1, ncol=2, byrow=TRUE))

# Metagenomes
par(mar=c(4,4,1,2), mgp=c(2.5, 0.75, 0), las=1, xaxs='i', yaxs='i')
plot(1, type='n', xlim=c(1,10), ylim=c(0,6), xaxt='n', yaxt='n',
     xlab='Sampling Depth', ylab='Detected Genes')
axis(1, at=seq(1,12,2), label=c(1,100,10000,1000000,100000000,1000000000), cex.axis=0.9)
minor.ticks.axis(1, 10, mn=0, mx=10)
axis(2, at=seq(0,6,2), label=c(1,100,10000,1000000), cex.axis=0.9, las=1)
minor.ticks.axis(2, 10, mn=0, mx=6)
legend('topleft', 'Metagenomes', bty='n', cex=1.2)
mtext('E', side=2, line=2, las=2, adj=2, padj=-10, cex=1.5, font=2)
legend('right', legend=c(expression(italic('C. difficile')),'Mock'), pch=c(17,19), pt.cex=1.5)
legend('bottomright', legend=c('Streptomycin-pretreated','Cefoperazone-pretreated','Clindamycin-pretreated'), 
       col=c(strep_col,cef_col,clinda_col), pch=22, pt.cex=0, lwd=3)
box()
# Add lines
lines(cef_630_metaG, col=cef_col, type='p', pch=17, cex=1.5)
lines(cef_630_metaG, col=cef_col, lwd=2.5)
lines(cef_mock_metaG, col=cef_col, type='p', pch=19, cex=1.5)
lines(cef_mock_metaG, col=cef_col, lwd=2.5)
lines(strep_630_metaG, col=strep_col, type='p', pch=17, cex=1.5)
lines(strep_630_metaG, col=strep_col, lwd=2.5)
lines(strep_mock_metaG, col=strep_col, type='p', pch=19, cex=1.5)
lines(strep_mock_metaG, col=strep_col, lwd=2.5)
lines(clinda_630_metaG, col=clinda_col, type='p', pch=17, cex=1.5)
lines(clinda_630_metaG, col=clinda_col, lwd=2.5)
lines(clinda_mock_metaG, col=clinda_col, type='p', pch=19, cex=1.5)
lines(clinda_mock_metaG, col=clinda_col, lwd=2.5)

# Metatranscriptomes
par(mar=c(4,4,1,2), mgp=c(2.5, 0.75, 0), las=1, xaxs='i', yaxs='i')
plot(1, type='n', xlim=c(1,10), ylim=c(0,6), xaxt='n', yaxt='n',
     xlab='Sampling Depth', ylab='Detected Genes')
axis(1, at=seq(1,12,2), label=c(1,100,10000,1000000,100000000,1000000000), cex.axis=0.9)
minor.ticks.axis(1, 10, mn=0, mx=10)
axis(2, at=seq(0,6,2), label=c(1,100,10000,1000000), cex.axis=0.9, las=1)
minor.ticks.axis(2, 10, mn=0, mx=6)
legend('topleft', 'Metatranscriptomes', bty='n', cex=1.2)
mtext('F', side=2, line=2, las=2, adj=2, padj=-10, cex=1.5, font=2)
legend('right', legend=c(expression(italic('C. difficile')),'Mock'), pch=c(17,19), pt.cex=1.5)
legend('bottomright', legend=c('Streptomycin-pretreated','Cefoperazone-pretreated','Clindamycin-pretreated'), 
       col=c(strep_col,cef_col,clinda_col), pch=22, pt.cex=0, lwd=3)
box()
# Add lines
lines(cef_630_metaT, col=cef_col, type='p', pch=17, cex=1.5)
lines(cef_630_metaT, col=cef_col, lwd=2.5)
lines(cef_mock_metaT, col=cef_col, type='p', pch=19, cex=1.5)
lines(cef_mock_metaT, col=cef_col, lwd=2.5)
lines(strep_630_metaT, col=strep_col, type='p', pch=17, cex=1.5)
lines(strep_630_metaT, col=strep_col, lwd=2.5)
lines(strep_mock_metaT, col=strep_col, type='p', pch=19, cex=1.5)
lines(strep_mock_metaT, col=strep_col, lwd=2.5)
lines(clinda_630_metaT, col=clinda_col, type='p', pch=17, cex=1.5)
lines(clinda_630_metaT, col=clinda_col, lwd=2.5)
lines(clinda_mock_metaT, col=clinda_col, type='p', pch=19, cex=1.5)
lines(clinda_mock_metaT, col=clinda_col, lwd=2.5)

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


