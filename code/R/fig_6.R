
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
metabolome$BIOCHEMICAL <- gsub('_', ' ', metabolome$BIOCHEMICAL)
substr(metabolome$BIOCHEMICAL, 1, 1) <- toupper(substr(metabolome$BIOCHEMICAL, 1, 1))
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$KEGG <- NULL
metabolome_annotation <- metabolome[,c(1:2)]
metabolome$SUB_PATHWAY <- NULL
metabolome <- metabolome[metabolome$SUPER_PATHWAY %in% c('Carbohydrate','Amino_Acid','Lipid'),] # Subset to amino acids and carbohydrates
metabolome <- metabolome[order(metabolome$SUPER_PATHWAY),]
metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.data.frame(t(metabolome))

# Format names and add metadata
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
cef_metabolome <- aggregate(cef_metabolome[,2:ncol(cef_metabolome)], by=list(cef_metabolome$infection), FUN=median)
rownames(cef_metabolome) <- c('Cefoperzone-Infected','Cefoperzone-Mock')
cef_metabolome$Group.1 <- NULL
clinda_metabolome <- aggregate(clinda_metabolome[,2:ncol(clinda_metabolome)], by=list(clinda_metabolome$infection), FUN=median)
rownames(clinda_metabolome) <- c('Clindamycin-Infected','Clindamycin-Mock')
clinda_metabolome$Group.1 <- NULL
strep_metabolome <- aggregate(strep_metabolome[,2:ncol(strep_metabolome)], by=list(strep_metabolome$infection), FUN=median)
rownames(strep_metabolome) <- c('Streptomycin-Infected','Streptomycin-Mock')
strep_metabolome$Group.1 <- NULL
abx_metabolome <- rbind(strep_metabolome,cef_metabolome,clinda_metabolome)

# Convert to matrices
metabolome <- as.matrix(metabolome)
abx_metabolome <- as.matrix(abx_metabolome)
noabx_metabolome <- as.matrix(noabx_metabolome)
strep_metabolome <- as.matrix(strep_metabolome)
cef_metabolome <- as.matrix(cef_metabolome)
clinda_metabolome <- as.matrix(clinda_metabolome)

# Remove columns with low variance
metabolome <- rm_lowVar(metabolome)
abx_metabolome <- rm_lowVar(abx_metabolome)
noabx_metabolome <- rm_lowVar(noabx_metabolome)
strep_metabolome <- rm_lowVar(strep_metabolome)
cef_metabolome <- rm_lowVar(cef_metabolome)
clinda_metabolome <- rm_lowVar(clinda_metabolome)

# Scale data and center columns
metabolome <- scale(metabolome)
abx_metabolome <- scale(abx_metabolome)
noabx_metabolome <- scale(noabx_metabolome)
strep_metabolome <- scale(strep_metabolome)
cef_metabolome <- scale(cef_metabolome)
clinda_metabolome <- scale(clinda_metabolome)

# Separate, cluster super pathways, and re-join
metabolome_path <- merge(metabolome_annotation, t(metabolome), by='row.names')
rownames(metabolome_path) <- metabolome_path$Row.names
metabolome_path$Row.names <- NULL
metabolome_path$SUB_PATHWAY <- NULL
metabolome_carb <- subset(metabolome_path, SUPER_PATHWAY == 'Carbohydrate')
metabolome_carb$SUPER_PATHWAY <- NULL
carb <- as.vector(hclust(dist(metabolome_carb))$order)
metabolome_carb <- t(metabolome_carb[carb,])
metabolome_amino <- subset(metabolome_path, SUPER_PATHWAY == 'Amino_Acid')
metabolome_amino$SUPER_PATHWAY <- NULL
amino <- as.vector(hclust(dist(metabolome_amino))$order)
metabolome_amino <- t(metabolome_amino[amino,])
metabolome_lipid <- subset(metabolome_path, SUPER_PATHWAY == 'Lipid')
metabolome_lipid$SUPER_PATHWAY <- NULL
lipid <- as.vector(hclust(dist(metabolome_lipid))$order)
metabolome_lipid <- t(metabolome_lipid[lipid,])
metabolome <- as.matrix(cbind(metabolome_carb, metabolome_amino, metabolome_lipid))
lengths <- c(length(carb), length(amino), length(lipid))
rm(metabolome_carb, metabolome_amino, metabolome_lipid,
   carb, amino, lipid, metabolome_path)

abx_metabolome_path <- merge(metabolome_annotation, t(abx_metabolome), by='row.names')
rownames(abx_metabolome_path) <- abx_metabolome_path$Row.names
abx_metabolome_path$Row.names <- NULL
abx_metabolome_path$SUB_PATHWAY <- NULL
abx_metabolome_carb <- subset(abx_metabolome_path, SUPER_PATHWAY == 'Carbohydrate')
abx_metabolome_carb$SUPER_PATHWAY <- NULL
carb <- as.vector(hclust(dist(abx_metabolome_carb))$order)
abx_metabolome_carb <- t(abx_metabolome_carb[carb,])
abx_metabolome_amino <- subset(abx_metabolome_path, SUPER_PATHWAY == 'Amino_Acid')
abx_metabolome_amino$SUPER_PATHWAY <- NULL
amino <- as.vector(hclust(dist(abx_metabolome_amino))$order)
abx_metabolome_amino <- t(abx_metabolome_amino[amino,])
abx_metabolome_lipid <- subset(abx_metabolome_path, SUPER_PATHWAY == 'Lipid')
abx_metabolome_lipid$SUPER_PATHWAY <- NULL
lipid <- as.vector(hclust(dist(abx_metabolome_lipid))$order)
abx_metabolome_lipid <- t(abx_metabolome_lipid[lipid,])
abx_metabolome <- as.matrix(cbind(abx_metabolome_carb, abx_metabolome_amino, abx_metabolome_lipid))
abx_lengths <- c(length(carb), length(amino), length(lipid))
rm(abx_metabolome_carb, abx_metabolome_amino, abx_metabolome_lipid,
   carb, amino, lipid, abx_metabolome_path, metabolome_annotation)

#-------------------------------------------------------------------------------------------------------------------------#

# Set plotting up parameters
heat_palette <- viridis(n=200)

# Generate figure
pdf(file='~/Desktop/test.pdf', width=20, height=20)


heatmap.2( abx_metabolome,
           col=heat_palette,
           trace='none',
           scale='none',
           symm=FALSE,
           symbreaks=FALSE,
           dendrogram='none',
           margins=c(10, 20),
           cexRow=2, 
           Rowv=FALSE,
           Colv=FALSE,
           labCol=FALSE,
           keysize=1,
           density.info='none',
           symkey=FALSE,
           key.xlab='Scaled Intensity',
           key.par=list(cex=1.5))
segments(x0=0.17, y0=0.03, x1=0.22, y1=0.03, lwd=5) # Carbs
text(x=0.19, y=0.01, 'Carbohydrates', cex=1.7)

segments(x0=0.23, y0=0.03, x1=0.45, y1=0.03, lwd=5) # Amino acids
text(x=0.34, y=0.01, 'Amino acids', cex=1.7)

segments(x0=0.46, y0=0.03, x1=0.845, y1=0.03, lwd=5) # Lipids
text(x=0.6525, y=0.01, 'Lipids', cex=1.7)


dev.off()



#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
#for (dep in deps){
#  pkg <- paste('package:', dep, sep='')
#  detach(pkg, character.only=TRUE)}
#setwd(starting_dir)
#rm(list=ls())
#gc()
