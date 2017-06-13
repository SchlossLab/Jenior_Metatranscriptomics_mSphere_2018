
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
plot_4 <- 'results/supplement/figures/figure_S4.pdf'
plot_5a <- 'results/supplement/figures/figure_S5a.pdf'
plot_5b <- 'results/supplement/figures/figure_S5b.pdf'
plot_5c <- 'results/supplement/figures/figure_S5c.pdf'

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
#metabolome <- metabolome[metabolome$SUPER_PATHWAY %in% c('Carbohydrate','Amino_Acid','Lipid'),] # Subset to amino acids and carbohydrates
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
rownames(noabx_metabolome) <- 'No_Antibiotics-Mock'
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

# Create and write supplementary table
supp_metabolome <- t(rbind(noabx_metabolome, abx_metabolome))
supp_metabolome <- merge(metabolome_annotation, supp_metabolome, by='row.names')
colnames(supp_metabolome)[1] <- 'BIOCHEMICAL'
supp_metabolome$BIOCHEMICAL <- gsub(' ', '_', supp_metabolome$BIOCHEMICAL)
write.table(supp_metabolome, file=supp_table3, sep='\t', row.names=FALSE, quote=FALSE)
rm(supp_metabolome)

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
rownames(metabolome) <- c('Cefoperazone - Infected 1','Cefoperazone - Infected 2','Cefoperazone - Infected 3',
                          'Cefoperazone - Infected 4','Cefoperazone - Infected 5','Cefoperazone - Infected 6',
                          'Cefoperazone - Infected 7','Cefoperazone - Infected 8','Cefoperazone - Infected 9',
                          'Cefoperazone - Mock 1','Cefoperazone - Mock 2','Cefoperazone - Mock 3',
                          'Cefoperazone - Mock 4','Cefoperazone - Mock 5','Cefoperazone - Mock 6',
                          'Cefoperazone - Mock 7','Cefoperazone - Mock 8','Cefoperazone - Mock 9',
                          'Clindamycin - Infected 1','Clindamycin - Infected 2','Clindamycin - Infected 3',
                          'Clindamycin - Infected 4','Clindamycin - Infected 5','Clindamycin - Infected 6',
                          'Clindamycin - Infected 7','Clindamycin - Infected 8','Clindamycin - Infected 9',
                          'Clindamycin - Mock 1','Clindamycin - Mock 2','Clindamycin - Mock 3',
                          'Clindamycin - Mock 4','Clindamycin - Mock 5','Clindamycin - Mock 6',
                          'Clindamycin - Mock 7','Clindamycin - Mock 8','Clindamycin - Mock 9',
                          'No Antibiotics - Mock 1','No Antibiotics - Mock 2','No Antibiotics - Mock 3',
                          'No Antibiotics - Mock 4','No Antibiotics - Mock 5','No Antibiotics - Mock 6',
                          'No Antibiotics - Mock 7','No Antibiotics - Mock 8','No Antibiotics - Mock 9',
                          'Streptomycin - Infected 1','Streptomycin - Infected 2','Streptomycin - Infected 3',
                          'Streptomycin - Infected 4','Streptomycin - Infected 5','Streptomycin - Infected 6',
                          'Streptomycin - Infected 7','Streptomycin - Infected 8','Streptomycin - Infected 9',
                          'Streptomycin - Mock 1','Streptomycin - Mock 2','Streptomycin - Mock 3',
                          'Streptomycin - Mock 4','Streptomycin - Mock 5','Streptomycin - Mock 6',
                          'Streptomycin - Mock 7','Streptomycin - Mock 8','Streptomycin - Mock 9')
metabolome <- metabolome[c('No Antibiotics - Mock 1','No Antibiotics - Mock 2','No Antibiotics - Mock 3',
                            'No Antibiotics - Mock 4','No Antibiotics - Mock 5','No Antibiotics - Mock 6',
                            'No Antibiotics - Mock 7','No Antibiotics - Mock 8','No Antibiotics - Mock 9',
                            'Streptomycin - Infected 1','Streptomycin - Infected 2','Streptomycin - Infected 3',
                            'Streptomycin - Infected 4','Streptomycin - Infected 5','Streptomycin - Infected 6',
                            'Streptomycin - Infected 7','Streptomycin - Infected 8','Streptomycin - Infected 9',
                            'Streptomycin - Mock 1','Streptomycin - Mock 2','Streptomycin - Mock 3',
                            'Streptomycin - Mock 4','Streptomycin - Mock 5','Streptomycin - Mock 6',
                            'Streptomycin - Mock 7','Streptomycin - Mock 8','Streptomycin - Mock 9',
                            'Cefoperazone - Infected 1','Cefoperazone - Infected 2','Cefoperazone - Infected 3',
                            'Cefoperazone - Infected 4','Cefoperazone - Infected 5','Cefoperazone - Infected 6',
                            'Cefoperazone - Infected 7','Cefoperazone - Infected 8','Cefoperazone - Infected 9',
                            'Cefoperazone - Mock 1','Cefoperazone - Mock 2','Cefoperazone - Mock 3',
                            'Cefoperazone - Mock 4','Cefoperazone - Mock 5','Cefoperazone - Mock 6',
                            'Cefoperazone - Mock 7','Cefoperazone - Mock 8','Cefoperazone - Mock 9',
                            'Clindamycin - Infected 1','Clindamycin - Infected 2','Clindamycin - Infected 3',
                            'Clindamycin - Infected 4','Clindamycin - Infected 5','Clindamycin - Infected 6',
                            'Clindamycin - Infected 7','Clindamycin - Infected 8','Clindamycin - Infected 9',
                            'Clindamycin - Mock 1','Clindamycin - Mock 2','Clindamycin - Mock 3',
                            'Clindamycin - Mock 4','Clindamycin - Mock 5','Clindamycin - Mock 6',
                            'Clindamycin - Mock 7','Clindamycin - Mock 8','Clindamycin - Mock 9'),]
metabolome <- merge(metabolome_annotation, t(metabolome), by='row.names')
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome <- metabolome[order(metabolome$SUPER_PATHWAY),]
pathways <- table(metabolome$SUPER_PATHWAY)
metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.matrix(t(metabolome))


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
rm(carb, amino, lipid, abx_metabolome_path, metabolome_annotation)

#-------------------------------------------------------------------------------------------------------------------------#

# Set plotting up parameters
heat_palette <- viridis(n=200)

# Generate figures
pdf(file=plot_4, width=50, height=30)
heatmap.2( metabolome,
           col=heat_palette,
           trace='none',
           scale='none',
           symm=FALSE,
           symbreaks=FALSE,
           dendrogram='none',
           margins=c(10, 20),
           cexRow=1.8, 
           Colv=FALSE,
           Rowv=FALSE,
           labCol=FALSE,
           keysize=1,
           density.info='none',
           symkey=FALSE,
           key.xlab='Scaled Intensity',
           key.par=list(cex=1.5))
segments(x0=0.182, y0=0.816, x1=0.182, y1=0.707, lwd=7) # Resistant
text(x=0.15, y=0.76, 'Resistant', cex=3)
segments(x0=0.182, y0=0.698, x1=0.182, y1=0.025, lwd=7) # Susceptible
text(x=0.15, y=0.36, 'Susceptible', cex=3)

segments(x0=0.189, y0=0.015, x1=0.326, y1=0.015, lwd=7) # amino acids
text(x=0.2585, y=0.005, 'Amino Acids', cex=2.2)
segments(x0=0.331, y0=0.015, x1=0.364, y1=0.015, lwd=7) # carbs
text(x=0.348, y=0.005, 'Carbohydrates', cex=2.2)
segments(x0=0.37, y0=0.015, x1=0.409, y1=0.015, lwd=7) # vit
text(x=0.39, y=0.005, 'Vitamins', cex=2.2)
segments(x0=0.413, y0=0.015, x1=0.423, y1=0.015, lwd=7) # Energy 
text(x=0.418, y=0.005, 'Energy', cex=2.2)
segments(x0=0.427, y0=0.015, x1=0.75, y1=0.015, lwd=7) # Lipid 
text(x=0.59, y=0.005, 'Lipid', cex=2.2)
segments(x0=0.756, y0=0.015, x1=0.809, y1=0.015, lwd=7) # Nucleotide 
text(x=0.784, y=0.005, 'Nucleotide', cex=2.2)
segments(x0=0.815, y0=0.015, x1=0.854, y1=0.015, lwd=7) # Peptide 
text(x=0.836, y=0.005, 'Peptide', cex=2.2)
segments(x0=0.86, y0=0.015, x1=0.94, y1=0.015, lwd=7) # Xenobiotics  
text(x=0.9, y=0.005, 'Xenobiotics ', cex=2.2)
dev.off()

pdf(file=plot_5a, width=20, height=20)
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
text(x=0.34, y=0.01, 'Amino Acids', cex=1.7)
segments(x0=0.46, y0=0.03, x1=0.845, y1=0.03, lwd=5) # Lipids
text(x=0.6525, y=0.01, 'Lipids', cex=1.7)
text(x=0.15, y=0.82, 'a', cex=4, font=2)
text(x=0.5, y=0.845, 'Susceptible Metabolomes Only', cex=3, font=2)
dev.off()

pdf(file=plot_5b, width=20, height=20)
heatmap.2( abx_metabolome_carb,
           col=heat_palette,
           trace='none',
           scale='none',
           symm=FALSE,
           symbreaks=FALSE,
           dendrogram='none',
           margins=c(35, 20),
           cexRow=2, 
           cexCol=2,
           Rowv=FALSE,
           Colv=FALSE,
           keysize=1,
           density.info='none',
           symkey=FALSE,
           key.xlab='Scaled Intensity',
           key.par=list(cex=1.5))
text(x=0.15, y=0.82, 'b', cex=4, font=2)
text(x=0.5, y=0.845, 'Carbohydrates', cex=3.5, font=2)
dev.off()

pdf(file=plot_5c, width=40, height=20)
heatmap.2( abx_metabolome_amino,
           col=heat_palette,
           trace='none',
           scale='none',
           symm=FALSE,
           symbreaks=FALSE,
           dendrogram='none',
           margins=c(30, 20),
           cexRow=2, 
           cexCol=1.8, 
           Rowv=FALSE,
           Colv=FALSE,
           keysize=1,
           density.info='none',
           symkey=FALSE,
           key.xlab='Scaled Intensity',
           key.par=list(cex=1.5))
text(x=0.17, y=0.82, 'c', cex=4, font=2)
text(x=0.5, y=0.845, 'Amino Acids', cex=3.5, font=2)
dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only=TRUE)}
setwd(starting_dir)
rm(list=ls())
gc()
