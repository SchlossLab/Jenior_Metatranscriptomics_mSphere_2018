
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define files
# Sequencing
cef_metatranscriptome <- 'data/read_mapping/cef_normalized_metaT.tsv'
strep_metatranscriptome <- 'data/read_mapping/strep_normalized_metaT.tsv'
clinda_metatranscriptome <- 'data/read_mapping/clinda_normalized_metaT.tsv'

# Metabolomes
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'

# Metadata
metadata <- 'data/metadata.tsv'

# Output plot
plot_file <- 'results/figures/figure_5.pdf'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data
# Normalized Metatranscriptomes
cef_normalized_reads <- read.delim(cef_normalized_reads, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=7)
clinda_normalized_reads <- read.delim(clinda_normalized_reads, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=7)
strep_normalized_reads <- read.delim(strep_normalized_reads, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=7)

# Metabolomes
metabolome <- read.delim(metabolome, sep='\t', header=TRUE)

# Metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data
# Metadata
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL

# Metabolomes
metabolome <- metabolome[,!colnames(metabolome) %in% c('GfC1M1','GfC1M2','GfC1M3', 
                                                       'GfC2M1','GfC2M2','GfC2M3', 
                                                       'GfC3M1','GfC3M2','GfC3M3', 
                                                       'GfC4M1','GfC4M2','GfC4M3', 
                                                       'GfC5M1','GfC5M2','GfC5M3', 
                                                       'GfC6M1','GfC6M2','GfC6M3')] # Germfree
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.data.frame(t(metabolome))
colnames(metabolome) <- gsub('_', ' ', colnames(metabolome))
substr(colnames(metabolome), 1, 1) <- toupper(substr(colnames(metabolome), 1, 1))

# Metatranscriptomes
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

# Screen out ribosomal genes
cef_annotated <- subset(cef_annotated, !grepl('Ribosomal_RNA*', cef_annotated$description))
cef_annotated <- subset(cef_annotated, !grepl('ribosomal_RNA*', cef_annotated$description))
cef_annotated <- subset(cef_annotated, !grepl('*ribosomal_RNA*', cef_annotated$description))
clinda_annotated <- subset(clinda_annotated, !grepl('Ribosomal_RNA*', clinda_annotated$description))
clinda_annotated <- subset(clinda_annotated, !grepl('ribosomal_RNA*', clinda_annotated$description))
clinda_annotated <- subset(clinda_annotated, !grepl('*ribosomal_RNA*', clinda_annotated$description))
strep_annotated <- subset(strep_annotated, !grepl('Ribosomal_RNA*', strep_annotated$description))
strep_annotated <- subset(strep_annotated, !grepl('ribosomal_RNA*', strep_annotated$description))
strep_annotated <- subset(strep_annotated, !grepl('*ribosomal_RNA*', strep_annotated$description))

# Remove hypothetical and uncharacterized annotations
cef_annotated <- subset(cef_annotated, description != 'hypothetical_protein')
cef_annotated <- subset(cef_annotated, description != 'uncharacterized_*')
clinda_annotated <- subset(clinda_annotated, description != 'hypothetical_protein')
clinda_annotated <- subset(clinda_annotated, description != 'uncharacterized_*')
strep_annotated <- subset(strep_annotated, description != 'hypothetical_protein')
strep_annotated <- subset(strep_annotated, description != 'uncharacterized_*')





#-------------------------------------------------------------------------------------------------------------------------#

# Generate figure

heat_palette <- colorRampPalette(c('red', 'white', 'blue'))(n = 299)


#heatmap.2(mat_data...


dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)}
setwd(starting_dir)
rm(list=ls())
gc()
