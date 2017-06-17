
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define files
# Metadata
metadata <- 'data/metadata.tsv'

# Mapping to CONCOCT clusters
lactobacillus_strep_630 <- 'data/read_mapping/concoct/streptomycin/lactobacillus_streptomycin_630.gene_reads.tsv'
lactobacillus_strep_mock <- 'data/read_mapping/concoct/streptomycin/lactobacillus_streptomycin_mock.gene_reads.tsv'
bacteroides_cef_630 <- 'data/read_mapping/concoct/cefoperazone/bacteroides_cefoperazone_630.gene_reads.tsv'
bacteroides_cef_mock <- 'data/read_mapping/concoct/cefoperazone/bacteroides_cefoperazone_mock.gene_reads.tsv'

# KEGG pathways
gene_paths <- 'data/gene_paths.tsv'

# Metabolomes
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'

# Output plot
plot_file <- 'results/figures/figure_5.pdf'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data
# Metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

# KEGG organism file
gene_paths <- read.delim(gene_paths, sep='\t', header=TRUE)

# cDNA mappings
lactobacillus_strep_630 <- read.delim(lactobacillus_strep_630, sep='\t', header=FALSE)
colnames(lactobacillus_strep_630) <- c('gene', 'reads')
lactobacillus_strep_mock <- read.delim(lactobacillus_strep_mock, sep='\t', header=FALSE)
colnames(lactobacillus_strep_mock) <- c('gene', 'reads')
bacteroides_cef_630 <- read.delim(bacteroides_cef_630, sep='\t', header=FALSE)
colnames(bacteroides_cef_630) <- c('gene', 'reads')
bacteroides_cef_mock <- read.delim(bacteroides_cef_mock, sep='\t', header=FALSE)
colnames(bacteroides_cef_mock) <- c('gene', 'reads')

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

#----------------#

# CONCOCT mappings 
lactobacillus_strep_630 <- merge(lactobacillus_strep_630, gene_paths, by='gene')
lactobacillus_strep_mock <- merge(lactobacillus_strep_mock, gene_paths, by='gene')
bacteroides_cef_630 <- merge(bacteroides_cef_630, gene_paths, by='gene')
bacteroides_cef_mock <- merge(bacteroides_cef_mock, gene_paths, by='gene')
rm(gene_paths)



# Evenly subsample reads
lactobacillus_size <- round(min(colSums(lactobacillus_strep[,c(1:2)]))*0.9) # Determine subsample level
lactobacillus_strep$reads_630 <- as.vector(rrarefy(lactobacillus_strep$reads_630, sample=lactobacillus_size))
lactobacillus_strep$reads_mock <- as.vector(rrarefy(lactobacillus_strep$reads_mock, sample=lactobacillus_size))
rm(lactobacillus_size)
bacteroides_size <- round(min(colSums(bacteroides_cef[,c(1:2)]))*0.9) # Determine subsample level
bacteroides_cef$reads_630 <- as.vector(rrarefy(bacteroides_cef$reads_630, sample=bacteroides_size))
bacteroides_cef$reads_mock <- as.vector(rrarefy(bacteroides_cef$reads_mock, sample=bacteroides_size))
rm(bacteroides_size)




# Aggregate by pathways
lactobacillus_strep_630 <- aggregate(lactobacillus_strep_630$reads, by=list(lactobacillus_strep_630$pathways), FUN=sum)
colnames(lactobacillus_strep_630) <- c('pathway','reads_630')
lactobacillus_strep_mock <- aggregate(lactobacillus_strep_mock$reads, by=list(lactobacillus_strep_mock$pathways), FUN=sum)
colnames(lactobacillus_strep_mock) <- c('pathway','reads_mock')
bacteroides_cef_630 <- aggregate(bacteroides_cef_630$reads, by=list(bacteroides_cef_630$pathways), FUN=sum)
colnames(bacteroides_cef_630) <- c('pathway','reads_630')
bacteroides_cef_mock <- aggregate(bacteroides_cef_mock$reads, by=list(bacteroides_cef_mock$pathways), FUN=sum)
colnames(bacteroides_cef_mock) <- c('pathway','reads_mock')

# Merge mock and infected groups
lactobacillus_strep <- merge(lactobacillus_strep_mock, lactobacillus_strep_630, by='pathway')
rm(lactobacillus_strep_630, lactobacillus_strep_mock)
rownames(lactobacillus_strep) <- lactobacillus_strep$pathway
lactobacillus_strep$pathway <- NULL
bacteroides_cef <- merge(bacteroides_cef_mock, bacteroides_cef_630, by='pathway')
rm(bacteroides_cef_630, bacteroides_cef_mock)
rownames(bacteroides_cef) <- bacteroides_cef$pathway
bacteroides_cef$pathway <- NULL


# Calculate differences between mock and infected metatranscriptomes
lactobacillus_strep$delta <- (lactobacillus_strep$reads_mock - lactobacillus_strep$reads_630) * -1
bacteroides_cef$delta <- (bacteroides_cef$reads_mock - bacteroides_cef$reads_630) * -1




#----------------#

# Metabolomes
metabolome$PUBCHEM <- NULL
metabolome$BIOCHEMICAL <- gsub('_', ' ', metabolome$BIOCHEMICAL)
substr(metabolome$BIOCHEMICAL, 1, 1) <- toupper(substr(metabolome$BIOCHEMICAL, 1, 1))
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$KEGG <- NULL
metabolite_annotation <- metabolome[,c(1:2)]
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.data.frame(t(metabolome))

# Format names and add metadata
colnames(metabolome) <- gsub('_', ' ', colnames(metabolome))
substr(colnames(metabolome), 1, 1) <- toupper(substr(colnames(metabolome), 1, 1))
metabolome <- clean_merge(metadata, metabolome)
rm(metadata)

# Subset antibiotics
cef_metabolome <- subset(metabolome, abx == 'cefoperazone')
cef_metabolome$abx <- NULL
strep_metabolome <- subset(metabolome, abx == 'streptomycin')
strep_metabolome$abx <- NULL
rm(metabolome)

# Remove metabolites with low variance
infection <- strep_metabolome$infection
strep_metabolome$infection <- NULL
strep_metabolome <- rm_lowVar(strep_metabolome)
strep_metabolome <- cbind(infection, strep_metabolome)
infection <- cef_metabolome$infection
cef_metabolome$infection <- NULL
cef_metabolome <- rm_lowVar(cef_metabolome)
cef_metabolome <- cbind(infection, cef_metabolome)

# Find medians within antibiotic groups
cef_metabolome <- aggregate(cef_metabolome[,2:ncol(cef_metabolome)], by=list(cef_metabolome$infection), FUN=median)
rownames(cef_metabolome) <- c('Cefoperzone_Infected','Cefoperzone_Mock')
cef_metabolome$Group.1 <- NULL
cef_metabolome<- t(cef_metabolome)
strep_metabolome <- aggregate(strep_metabolome[,2:ncol(strep_metabolome)], by=list(strep_metabolome$infection), FUN=median)
rownames(strep_metabolome) <- c('Streptomycin_Infected','Streptomycin_Mock')
strep_metabolome$Group.1 <- NULL
strep_metabolome <- t(strep_metabolome)

# Calculate delta between infected and mock
cef_metabolome <- as.data.frame(10 ^ cef_metabolome)
cef_metabolome$delta <- (cef_metabolome$Cefoperzone_Mock - cef_metabolome$Cefoperzone_Infected) * -1
strep_metabolome <- as.data.frame(10 ^ strep_metabolome)
strep_metabolome$delta <- (strep_metabolome$Streptomycin_Mock - strep_metabolome$Streptomycin_Infected) * -1

# Add back pathway annotations
strep_metabolome <- merge(strep_metabolome, metabolite_annotation, by='row.names', all.x=TRUE)
strep_metabolome <- strep_metabolome[order(strep_metabolome$Row.names),]
rownames(strep_metabolome) <- strep_metabolome$Row.names
strep_metabolome$Row.names <- NULL
strep_metabolome$SUB_PATHWAY <- NULL
cef_metabolome <- merge(cef_metabolome, metabolite_annotation, by='row.names', all.x=TRUE)
cef_metabolome <- cef_metabolome[order(cef_metabolome$Row.names),]
rownames(cef_metabolome) <- cef_metabolome$Row.names
cef_metabolome$Row.names <- NULL
cef_metabolome$SUB_PATHWAY <- NULL
rm(metabolite_annotation)

#----------------#

# Create named vectors for each delta dataset
# Metatranscriptomes
strep_metatranscriptome_cdi <- strep_metatranscriptome$strep_630_metaT_reads
names(strep_metatranscriptome_cdi) <- strep_metatranscriptome
strep_metatranscriptome_mock <- strep_metatranscriptome$strep_mock_metaT_reads
names(strep_metatranscriptome_mock) <- strep_metatranscriptome$genus
rm(strep_metatranscriptome)
cef_metatranscriptome_cdi <- cef_metatranscriptome$cef_630_metaT_reads
names(cef_metatranscriptome_cdi) <- cef_metatranscriptome$genus
cef_metatranscriptome_mock <- cef_metatranscriptome$cef_mock_metaT_reads
names(cef_metatranscriptome_mock) <- cef_metatranscriptome$genus
rm(cef_metatranscriptome)
clinda_metatranscriptome_cdi <- clinda_metatranscriptome$clinda_630_metaT_reads
names(clinda_metatranscriptome_cdi) <- clinda_metatranscriptome$genus
clinda_metatranscriptome_mock <- clinda_metatranscriptome$clinda_mock_metaT_reads
names(clinda_metatranscriptome_mock) <- clinda_metatranscriptome$genus
rm(clinda_metatranscriptome)
# Metabolomes
strep_metabolome_cdi <- strep_metabolome$Streptomycin_Infected
names(strep_metabolome_cdi) <- strep_metabolome$SUPER_PATHWAY
strep_metabolome_mock <- strep_metabolome$Streptomycin_Mock
names(strep_metabolome_mock) <- strep_metabolome$SUPER_PATHWAY
rm(strep_metabolome)
cef_metabolome_cdi <- cef_metabolome$Cefoperzone_Infected
names(cef_metabolome_cdi) <- cef_metabolome$SUPER_PATHWAY
cef_metabolome_mock <- cef_metabolome$Cefoperzone_Mock
names(cef_metabolome_mock) <- cef_metabolome$SUPER_PATHWAY
rm(cef_metabolome)
clinda_metabolome_cdi <- clinda_metabolome$Clindamycin_Infected
names(clinda_metabolome_cdi) <- clinda_metabolome$SUPER_PATHWAY
clinda_metabolome_mock <- clinda_metabolome$Clindamycin_Mock
names(clinda_metabolome_mock) <- clinda_metabolome$SUPER_PATHWAY
rm(clinda_metabolome)


# Correlate amounts of change between metatranscriptome and metabolome

test1 <- c()
test2 <- c()
for (i in 1:length(strep_metatranscriptome_cdi)){
  test1[i] <- strep_metatranscriptome_cdi[i]
  test2[i] <- strep_metabolome_cdi[1]
}
test <- cor(x=test1, y=test2, method='spearman')


strep_correlation <- cor(x=strep_metatranscriptome_cdi, y=strep_metabolome_cdi, method='spearman', use='all.obs')
cef_correlation <- cor(x=cef_metatranscriptome_cdi, y=cef_metabolome_cdi, method='spearman')
clinda_correlation <- cor(x=clinda_metatranscriptome_cdi, y=clinda_metabolome_cdi, method='spearman')


d <- dist(as.matrix(mtcars))   # find distance matrix 
hc <- hclust(d)                # apply hirarchical clustering 








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
