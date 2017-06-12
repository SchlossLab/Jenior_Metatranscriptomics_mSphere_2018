
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define files
# Metadata
metadata <- 'data/metadata.tsv'

# Sequencing
strep_metatranscriptome <- 'data/read_mapping/strep_normalized_metaT.tsv'
cef_metatranscriptome <- 'data/read_mapping/cef_normalized_metaT.tsv'
clinda_metatranscriptome <- 'data/read_mapping/clinda_normalized_metaT.tsv'

# Metabolomes
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'

# Output plot
plot_file <- 'results/figures/figure_5.pdf'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data
# Metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

# Normalized Metatranscriptomes
strep_metatranscriptome <- read.delim(strep_metatranscriptome, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=7)
cef_metatranscriptome <- read.delim(cef_metatranscriptome, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=7)
clinda_metatranscriptome <- read.delim(clinda_metatranscriptome, sep='\t', header=TRUE, stringsAsFactors=FALSE, row.names=7)

# Metabolomes
metabolome <- read.delim(metabolome, sep='\t', header=TRUE)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data
# Metadata
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$susceptibility <- NULL

# Metatranscriptomes
# Remove C. difficile genes
cdiff_omit <- c('cdf','pdc','cdc','cdl','pdf')
cef_metatranscriptome <- cef_metatranscriptome[!cef_metatranscriptome$organism %in% cdiff_omit,]
clinda_metatranscriptome <- clinda_metatranscriptome[!clinda_metatranscriptome$organism %in% cdiff_omit,]
strep_metatranscriptome <- strep_metatranscriptome[!strep_metatranscriptome$organism %in% cdiff_omit,]
rm(cdiff_omit)

# Screen for those genes that have a gene annotation
cef_annotated <- cef_metatranscriptome[!rownames(cef_metatranscriptome) %in% rownames(cef_metatranscriptome[grep('unknown_\\d', rownames(cef_metatranscriptome)),]), ]
clinda_annotated <- clinda_metatranscriptome[!rownames(clinda_metatranscriptome) %in% rownames(clinda_metatranscriptome[grep('unknown_\\d', rownames(clinda_metatranscriptome)),]), ]
strep_annotated <- strep_metatranscriptome[!rownames(strep_metatranscriptome) %in% rownames(strep_metatranscriptome[grep('unknown_\\d', rownames(strep_metatranscriptome)),]), ]
rm(cef_metatranscriptome, clinda_metatranscriptome, strep_metatranscriptome)

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

# Metabolomes
metabolome$PUBCHEM <- NULL
metabolome_annotation <- metabolome[,c(1:4)]
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome <- metabolome[metabolome$SUPER_PATHWAY %in% c('Carbohydrate','Amino_Acid','Lipid'),] # Subset to amino acids and carbohydrates

metabolome <- metabolome[order(metabolome$SUPER_PATHWAY),]

metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.data.frame(t(metabolome))
colnames(metabolome) <- gsub('_', ' ', colnames(metabolome))
substr(colnames(metabolome), 1, 1) <- toupper(substr(colnames(metabolome), 1, 1))
metabolome <- clean_merge(metadata, metabolome)
metabolome <- subset(metabolome, abx != 'germfree')

# Subset antibiotics
noabx_metabolome <- subset(metabolome, abx == 'none')
noabx_metabolome$abx <- NULL
cef_metabolome <- subset(metabolome, abx == 'cefoperazone')
cef_metabolome$abx <- NULL
clinda_metabolome <- subset(metabolome, abx == 'clindamycin')
clinda_metabolome$abx <- NULL
strep_metabolome <- subset(metabolome, abx == 'streptomycin')
strep_metabolome$abx <- NULL
metabolome$abx <- NULL
metabolome$infection <- NULL

# Find medians within antibiotic groups
noabx_metabolome <- aggregate(noabx_metabolome[,2:ncol(noabx_metabolome)], by=list(noabx_metabolome$infection), FUN=median)
rownames(noabx_metabolome) <- noabx_metabolome$Group.1
noabx_metabolome$Group.1 <- NULL
noabx_metabolome <- as.matrix(noabx_metabolome)
cef_metabolome <- aggregate(cef_metabolome[,2:ncol(cef_metabolome)], by=list(cef_metabolome$infection), FUN=median)
rownames(cef_metabolome) <- cef_metabolome$Group.1
cef_metabolome$Group.1 <- NULL
cef_metabolome <- as.matrix(cef_metabolome)
clinda_metabolome <- aggregate(clinda_metabolome[,2:ncol(clinda_metabolome)], by=list(clinda_metabolome$infection), FUN=median)
rownames(clinda_metabolome) <- clinda_metabolome$Group.1
clinda_metabolome$Group.1 <- NULL
clinda_metabolome <- as.matrix(clinda_metabolome)
strep_metabolome <- aggregate(strep_metabolome[,2:ncol(strep_metabolome)], by=list(strep_metabolome$infection), FUN=median)
rownames(strep_metabolome) <- strep_metabolome$Group.1
strep_metabolome$Group.1 <- NULL
strep_metabolome <- as.matrix(strep_metabolome)
metabolome <- as.matrix(metabolome)

#-------------------------------------------------------------------------------------------------------------------------#

# Set plotting up parameters
heat_palette <- viridis(n=200)
#heat_palette <- colorRampPalette(c('red', 'white', 'blue'))(n=100)


test <- metabolome[,complete.cases(t(metabolome))] 
test <- scale(test) # scale and center columns




# Generate figure
pdf(file='~/Desktop/test.pdf', width=20, height=20)
heatmap.2( test,
           col=heat_palette,
           margins=c(20, 20),
           trace='none',
           density.info='none',
           xlab='Comparison',
           scale='none',
           #key=FALSE,
           symbreaks=min(metabolome, na.rm=TRUE),
           na.color=heat_palette[100],
           cexRow=0.9, cexCol=0.7,
           dendrogram='row', 
           Colv=FALSE )

dev.off()



#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only=TRUE)}
setwd(starting_dir)
rm(list=ls())
gc()
