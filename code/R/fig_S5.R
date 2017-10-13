
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define files
# Metadata
metadata <- 'data/metadata.tsv'

# Metabolomes
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'

# Output plot
plot_file <- 'results/supplement/figures/figure_S5.pdf'

# Supplementary table
supp_table3 <- 'results/supplement/tables/table_S4.tsv'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data
# Metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

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

# Metabolomes
metabolome$PUBCHEM <- NULL
metabolome$BIOCHEMICAL <- gsub('_', ' ', metabolome$BIOCHEMICAL)
substr(metabolome$BIOCHEMICAL, 1, 1) <- toupper(substr(metabolome$BIOCHEMICAL, 1, 1))
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$KEGG <- NULL
metabolome_annotation <- metabolome[,c(1:2)]
metabolome$SUB_PATHWAY <- NULL
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
rownames(noabx_metabolome) <- 'No_Antibiotics-Mock(median)'
noabx_metabolome$Group.1 <- NULL
cef_metabolome <- aggregate(cef_metabolome[,2:ncol(cef_metabolome)], by=list(cef_metabolome$infection), FUN=median)
rownames(cef_metabolome) <- c('Cefoperzone-Infected(median)','Cefoperzone-Mock(median)')
cef_metabolome$Group.1 <- NULL
clinda_metabolome <- aggregate(clinda_metabolome[,2:ncol(clinda_metabolome)], by=list(clinda_metabolome$infection), FUN=median)
rownames(clinda_metabolome) <- c('Clindamycin-Infected(median)','Clindamycin-Mock(median)')
clinda_metabolome$Group.1 <- NULL
strep_metabolome <- aggregate(strep_metabolome[,2:ncol(strep_metabolome)], by=list(strep_metabolome$infection), FUN=median)
rownames(strep_metabolome) <- c('Streptomycin-Infected(median)','Streptomycin-Mock(median)')
strep_metabolome$Group.1 <- NULL
abx_metabolome <- rbind(strep_metabolome,cef_metabolome,clinda_metabolome)

# Create and write supplementary table
supp_metabolome <- t(rbind(noabx_metabolome, abx_metabolome))
supp_metabolome <- merge(metabolome_annotation, supp_metabolome, by='row.names')
colnames(supp_metabolome)[1] <- 'BIOCHEMICAL'
supp_metabolome$BIOCHEMICAL <- gsub(' ', '_', supp_metabolome$BIOCHEMICAL)
write.table(supp_metabolome, file=supp_table3, sep='\t', row.names=FALSE, quote=FALSE)
rm(supp_metabolome)
rownames(abx_metabolome) <- c('Streptomycin\nInfected','Streptomycin\nMock','Cefoperzone\nInfected',
                              'Cefoperzone\nMock','Clindamycin\nInfected','Clindamycin\nMock')

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
rownames(metabolome) <- c('Cef - Infected 1','Cef - Infected 2','Cef - Infected 3',
                          'Cef - Infected 4','Cef - Infected 5','Cef - Infected 6',
                          'Cef - Infected 7','Cef - Infected 8','Cef - Infected 9',
                          'Cef - Mock 1','Cef - Mock 2','Cef - Mock 3',
                          'Cef - Mock 4','Cef - Mock 5','Cef - Mock 6',
                          'Cef - Mock 7','Cef - Mock 8','Cef - Mock 9',
                          'Clinda - Infected 1','Clinda - Infected 2','Clinda - Infected 3',
                          'Clinda - Infected 4','Clinda - Infected 5','Clinda - Infected 6',
                          'Clinda - Infected 7','Clinda - Infected 8','Clinda - Infected 9',
                          'Clinda - Mock 1','Clinda - Mock 2','Clinda - Mock 3',
                          'Clinda - Mock 4','Clinda - Mock 5','Clinda - Mock 6',
                          'Clinda - Mock 7','Clinda - Mock 8','Clinda - Mock 9',
                          'No Abx - Mock 1','No Abx - Mock 2','No Abx - Mock 3',
                          'No Abx - Mock 4','No Abx - Mock 5','No Abx - Mock 6',
                          'No Abx - Mock 7','No Abx - Mock 8','No Abx - Mock 9',
                          'Strep - Infected 1','Strep - Infected 2','Strep - Infected 3',
                          'Strep - Infected 4','Strep - Infected 5','Strep - Infected 6',
                          'Strep - Infected 7','Strep - Infected 8','Strep - Infected 9',
                          'Strep - Mock 1','Strep - Mock 2','Strep - Mock 3',
                          'Strep - Mock 4','Strep - Mock 5','Strep - Mock 6',
                          'Strep - Mock 7','Strep - Mock 8','Strep - Mock 9')
metabolome <- metabolome[c('No Abx - Mock 1','No Abx - Mock 2','No Abx - Mock 3',
                            'No Abx - Mock 4','No Abx - Mock 5','No Abx - Mock 6',
                            'No Abx - Mock 7','No Abx - Mock 8','No Abx - Mock 9',
                            'Strep - Infected 1','Strep - Infected 2','Strep - Infected 3',
                            'Strep - Infected 4','Strep - Infected 5','Strep - Infected 6',
                            'Strep - Infected 7','Strep - Infected 8','Strep - Infected 9',
                            'Strep - Mock 1','Strep - Mock 2','Strep - Mock 3',
                            'Strep - Mock 4','Strep - Mock 5','Strep - Mock 6',
                            'Strep - Mock 7','Strep - Mock 8','Strep - Mock 9',
                            'Cef - Infected 1','Cef - Infected 2','Cef - Infected 3',
                            'Cef - Infected 4','Cef - Infected 5','Cef - Infected 6',
                            'Cef - Infected 7','Cef - Infected 8','Cef - Infected 9',
                            'Cef - Mock 1','Cef - Mock 2','Cef - Mock 3',
                            'Cef - Mock 4','Cef - Mock 5','Cef - Mock 6',
                            'Cef - Mock 7','Cef - Mock 8','Cef - Mock 9',
                            'Clinda - Infected 1','Clinda - Infected 2','Clinda - Infected 3',
                            'Clinda - Infected 4','Clinda - Infected 5','Clinda - Infected 6',
                            'Clinda - Infected 7','Clinda - Infected 8','Clinda - Infected 9',
                            'Clinda - Mock 1','Clinda - Mock 2','Clinda - Mock 3',
                            'Clinda - Mock 4','Clinda - Mock 5','Clinda - Mock 6',
                            'Clinda - Mock 7','Clinda - Mock 8','Clinda - Mock 9'),]
metabolome <- merge(metabolome_annotation, t(metabolome), by='row.names')
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
metabolome$SUB_PATHWAY <- NULL
pathways <- table(metabolome$SUPER_PATHWAY)
metabolome_amino <- subset(metabolome, SUPER_PATHWAY == 'Amino_Acid')
metabolome_amino$SUPER_PATHWAY <- NULL
amino <- as.vector(hclust(dist(metabolome_amino))$order)
metabolome_amino <- t(metabolome_amino[amino,])
metabolome_carb <- subset(metabolome, SUPER_PATHWAY == 'Carbohydrate')
metabolome_carb$SUPER_PATHWAY <- NULL
carb <- as.vector(hclust(dist(metabolome_carb))$order)
metabolome_carb <- t(metabolome_carb[carb,])
metabolome_vit <- subset(metabolome, SUPER_PATHWAY == 'Cofactors_and_Vitamins')
metabolome_vit$SUPER_PATHWAY <- NULL
vit <- as.vector(hclust(dist(metabolome_vit))$order)
metabolome_vit <- t(metabolome_vit[vit,])
metabolome_energy <- subset(metabolome, SUPER_PATHWAY == 'Energy')
metabolome_energy$SUPER_PATHWAY <- NULL
energy <- as.vector(hclust(dist(metabolome_energy))$order)
metabolome_energy <- t(metabolome_energy[energy,])
metabolome_lipid <- subset(metabolome, SUPER_PATHWAY == 'Lipid')
metabolome_lipid$SUPER_PATHWAY <- NULL
lipid <- as.vector(hclust(dist(metabolome_lipid))$order)
metabolome_lipid <- t(metabolome_lipid[lipid,])
metabolome_nuc <- subset(metabolome, SUPER_PATHWAY == 'Nucleotide')
metabolome_nuc$SUPER_PATHWAY <- NULL
nuc <- as.vector(hclust(dist(metabolome_nuc))$order)
metabolome_nuc <- t(metabolome_nuc[nuc,])
metabolome_pep <- subset(metabolome, SUPER_PATHWAY == 'Peptide')
metabolome_pep$SUPER_PATHWAY <- NULL
pep <- as.vector(hclust(dist(metabolome_pep))$order)
metabolome_pep <- t(metabolome_pep[pep,])
metabolome_xeno <- subset(metabolome, SUPER_PATHWAY == 'Xenobiotics')
metabolome_xeno$SUPER_PATHWAY <- NULL
xeno <- as.vector(hclust(dist(metabolome_xeno))$order)
metabolome_xeno <- t(metabolome_xeno[xeno,])
metabolome <- as.matrix(cbind(metabolome_amino,metabolome_carb,metabolome_vit,
                              metabolome_energy,metabolome_lipid,metabolome_nuc,
                              metabolome_pep,metabolome_xeno))
rm(metabolome_amino,metabolome_carb,metabolome_vit,
   metabolome_energy,metabolome_lipid,metabolome_nuc,
   metabolome_pep,metabolome_xeno)

abx_metabolome <- merge(metabolome_annotation, t(abx_metabolome), by='row.names')
rownames(abx_metabolome) <- abx_metabolome$Row.names
abx_metabolome$Row.names <- NULL
abx_metabolome$SUB_PATHWAY <- NULL
pathways <- table(abx_metabolome$SUPER_PATHWAY)
abx_metabolome_amino <- subset(abx_metabolome, SUPER_PATHWAY == 'Amino_Acid')
abx_metabolome_amino$SUPER_PATHWAY <- NULL
amino <- as.vector(hclust(dist(abx_metabolome_amino))$order)
abx_metabolome_amino <- t(abx_metabolome_amino[amino,])
abx_metabolome_carb <- subset(abx_metabolome, SUPER_PATHWAY == 'Carbohydrate')
abx_metabolome_carb$SUPER_PATHWAY <- NULL
carb <- as.vector(hclust(dist(abx_metabolome_carb))$order)
abx_metabolome_carb <- t(abx_metabolome_carb[carb,])
abx_metabolome_vit <- subset(abx_metabolome, SUPER_PATHWAY == 'Cofactors_and_Vitamins')
abx_metabolome_vit$SUPER_PATHWAY <- NULL
vit <- as.vector(hclust(dist(abx_metabolome_vit))$order)
abx_metabolome_vit <- t(abx_metabolome_vit[vit,])
abx_metabolome_energy <- subset(abx_metabolome, SUPER_PATHWAY == 'Energy')
abx_metabolome_energy$SUPER_PATHWAY <- NULL
energy <- as.vector(hclust(dist(abx_metabolome_energy))$order)
abx_metabolome_energy <- t(abx_metabolome_energy[energy,])
abx_metabolome_lipid <- subset(abx_metabolome, SUPER_PATHWAY == 'Lipid')
abx_metabolome_lipid$SUPER_PATHWAY <- NULL
lipid <- as.vector(hclust(dist(abx_metabolome_lipid))$order)
abx_metabolome_lipid <- t(abx_metabolome_lipid[lipid,])
abx_metabolome_nuc <- subset(abx_metabolome, SUPER_PATHWAY == 'Nucleotide')
abx_metabolome_nuc$SUPER_PATHWAY <- NULL
nuc <- as.vector(hclust(dist(abx_metabolome_nuc))$order)
abx_metabolome_nuc <- t(abx_metabolome_nuc[nuc,])
abx_metabolome_pep <- subset(abx_metabolome, SUPER_PATHWAY == 'Peptide')
abx_metabolome_pep$SUPER_PATHWAY <- NULL
pep <- as.vector(hclust(dist(abx_metabolome_pep))$order)
abx_metabolome_pep <- t(abx_metabolome_pep[pep,])
abx_metabolome_xeno <- subset(abx_metabolome, SUPER_PATHWAY == 'Xenobiotics')
abx_metabolome_xeno$SUPER_PATHWAY <- NULL
xeno <- as.vector(hclust(dist(abx_metabolome_xeno))$order)
abx_metabolome_xeno <- t(abx_metabolome_xeno[xeno,])
abx_metabolome <- as.matrix(cbind(abx_metabolome_amino,abx_metabolome_carb,abx_metabolome_vit,
                              abx_metabolome_energy,abx_metabolome_lipid,abx_metabolome_nuc,
                              abx_metabolome_pep,abx_metabolome_xeno))
rm(abx_metabolome_vit,metabolome_annotation, amino,carb,vit,energy,lipid,nuc,pep,xeno,
   abx_metabolome_energy,abx_metabolome_lipid,abx_metabolome_nuc,
   abx_metabolome_pep,abx_metabolome_xeno)

#-------------------------------------------------------------------------------------------------------------------------#

# Generate figures
pdf(file=plot_file, width=50, height=30)
heatmap.2( metabolome,
           col=heat_palette,
           trace='none',
           scale='none',
           symm=FALSE,
           symbreaks=FALSE,
           dendrogram='none',
           margins=c(10, 20),
           cexRow=2.5, 
           Colv=FALSE,
           Rowv=FALSE,
           labCol=FALSE,
           keysize=0.75,
           density.info='none',
           symkey=FALSE,
           key.xlab='Scaled Intensity',
           key.par=list(cex=2.5))
segments(x0=0.14, y0=0.86, x1=0.14, y1=0.745, lwd=7) # Resistant
text(x=0.11, y=0.8, 'Resistant', cex=3)
segments(x0=0.14, y0=0.738, x1=0.14, y1=0.025, lwd=7) # Susceptible
text(x=0.11, y=0.36, 'Susceptible', cex=3)

segments(x0=0.145, y0=0.017, x1=0.305, y1=0.017, lwd=7) # amino acids
text(x=0.23, y=0.005, 'Amino Acids', cex=2.5)

segments(x0=0.31, y0=0.017, x1=0.37, y1=0.017, lwd=7) # carbs
text(x=0.34, y=0.005, 'Carbohydrates', cex=2.5)

segments(x0=0.375, y0=0.017, x1=0.42, y1=0.017, lwd=7) # vit
text(x=0.4, y=0.005, 'Vitamins', cex=2.5)

segments(x0=0.425, y0=0.017, x1=0.445, y1=0.017, lwd=7) # Energy 
text(x=0.4375, y=0.005, 'Energy', cex=2.5)

segments(x0=0.45, y0=0.017, x1=0.74, y1=0.017, lwd=7) # Lipid 
text(x=0.59, y=0.005, 'Lipid', cex=2.5)

segments(x0=0.745, y0=0.017, x1=0.809, y1=0.017, lwd=7) # Nucleotide 
text(x=0.784, y=0.005, 'Nucleotide', cex=2.5)

segments(x0=0.815, y0=0.017, x1=0.854, y1=0.017, lwd=7) # Peptide 
text(x=0.836, y=0.005, 'Peptide', cex=2.5)

segments(x0=0.86, y0=0.017, x1=0.94, y1=0.017, lwd=7) # Xenobiotics  
text(x=0.9, y=0.005, 'Xenobiotics ', cex=2.2)
dev.off()


#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only=TRUE)}
setwd(starting_dir)
rm(list=ls())
gc()
