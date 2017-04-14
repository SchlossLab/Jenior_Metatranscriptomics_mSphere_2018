# Set up environment
rm(list=ls())
gc()

# Load in functions
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Output plot name
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/supplement/figures/figure_S2.pdf'

# Input 0.03 OTU shared file
shared_otu_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.unique_list.conventional.shared'

# Input Metabolomes
metabolome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolome/metabolomics.tsv'

# Input Metadata
metadata_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metadata.tsv'

#----------------#

# Read in data

# 16S data
shared_otu <- read.delim(shared_otu_file, sep='\t', header=TRUE, row.names=2)
rm(shared_otu_file)

# Metabolomes
metabolome <- read.delim(metabolome_file, sep='\t', header=TRUE)
rm(metabolome_file)

# Metadata
metadata <- read.delim(metadata_file, sep='\t', header=T, row.names=1)
rm(metadata_file)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data

# Metadata
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2'), ] # Remove contaminated sample
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$susceptibility <- NULL

# Metabolomes
metabolome <- metabolome[,!colnames(metabolome) %in% c('CefC5M2')] # Contaminated sample
metabolome <- metabolome[,!colnames(metabolome) %in% c('GfC1M1','GfC1M2','GfC1M3',
                                                      'GfC2M1','GfC2M2','GfC2M3',
                                                      'GfC3M1','GfC3M2','GfC3M3',
                                                      'GfC4M1','GfC4M2','GfC4M3',
                                                      'GfC5M1','GfC5M2','GfC5M3',
                                                      'GfC6M1','GfC6M2','GfC6M3')] # Germfree samples
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('StrepC4M1','StrepC4M2','StrepC4M3'), ] 
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
bile_metabolome <- rbind(subset(metabolome, SUB_PATHWAY == 'Primary_Bile_Acid_Metabolism'),
                         subset(metabolome, SUB_PATHWAY == 'Secondary_Bile_Acid_Metabolism'))
bile_metabolome$SUB_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
carb_metabolome <- subset(metabolome, SUPER_PATHWAY == 'Carbohydrate')
carb_metabolome$SUPER_PATHWAY <- NULL
carb_metabolome <- t(carb_metabolome)
aa_metabolome <- subset(metabolome, SUPER_PATHWAY == 'Amino_Acid')
aa_metabolome$SUPER_PATHWAY <- NULL
aa_metabolome <- t(aa_metabolome)
metabolome$SUPER_PATHWAY <- NULL
metabolome <- t(metabolome)

# 16S data
shared_otu$label <- NULL
shared_otu$numOtus <- NULL
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('CefC5M2'), ] # Contaminated sample
shared_otu <- shared_otu[,!(names(shared_otu) %in% c('Otu0004','Otu0308'))] # Remove residual C. difficile OTUs
shared_otu <- rrarefy(shared_otu, ceiling(min(rowSums(shared_otu)) * 0.9))
shared_otu <- filter_table(shared_otu)
metabolome_16s <- clean_merge(metabolome, shared_otu)
metabolome_16s <- clean_merge(metadata, metabolome_16s)

#-------------------------------------------------------------------------------------------------------------------------#

# Calculate axes

# Metabolome
metabolome_nmds <- metaMDS(metabolome, k=2, trymax=100)$points
metabolome_nmds[,1] <- metabolome_nmds[,1] + 0.15
metabolome_nmds[,2] <- metabolome_nmds[,2] * -1
metabolome_nmds <- clean_merge(metadata, metabolome_nmds)
metabolome_cefoperazone <- subset(metabolome_nmds, abx == 'cefoperazone')
metabolome_clindamycin <- subset(metabolome_nmds, abx == 'clindamycin')
metabolome_streptomycin <- subset(metabolome_nmds, abx == 'streptomycin')
metabolome_untreated <- subset(metabolome_nmds, abx == 'none')

# Community structure
otu_nmds <- metaMDS(shared_otu, k=2, trymax=100)$points
otu_nmds[,1] <- otu_nmds[,1] - 0.2
otu_nmds[,2] <- otu_nmds[,2] - 0.2
otu_nmds <- clean_merge(metadata, otu_nmds)
otu_cefoperazone <- subset(otu_nmds, abx == 'cefoperazone')
otu_clindamycin <- subset(otu_nmds, abx == 'clindamycin')
otu_streptomycin <- subset(otu_nmds, abx == 'streptomycin')
otu_untreated <- subset(otu_nmds, abx == 'none')

#----------------#

# Test differences between untreated and all other groups
shared_otu <- clean_merge(metadata, shared_otu)
shared_otu$infection <- NULL

# Subset for testing
cef_otu <- rbind(subset(shared_otu, abx == 'cefoperazone'), subset(shared_otu, abx == 'none'))
cef_abx <- cef_otu$abx
cef_otu$abx <- NULL
clinda_otu <- rbind(subset(shared_otu, abx == 'clindamycin'), subset(shared_otu, abx == 'none'))
clinda_abx <- clinda_otu$abx
clinda_otu$abx <- NULL
strep_otu <- rbind(subset(shared_otu, abx == 'streptomycin'), subset(shared_otu, abx == 'none'))
strep_abx <- strep_otu$abx
strep_otu$abx <- NULL

# Calculate significant differences
cef_p <- anosim(cef_otu, cef_abx, permutations=999, distance='bray')$signif
clinda_p <- anosim(clinda_otu, clinda_abx, permutations=999, distance='bray')$signif
strep_p <- anosim(strep_otu, strep_abx, permutations=999, distance='bray')$signif

rm(shared_otu, cef_otu, cef_abx, clinda_otu, clinda_abx, strep_otu, strep_abx)

#----------------#

# Subset for feature selection
all_infection <- metabolome_16s
all_infection$abx <- NULL
conv_infection <- subset(metabolome_16s, abx != 'germfree')
abx_infection <- subset(metabolome_16s, abx != 'none')
conv_infection$abx <- NULL
abx_infection$abx <- NULL
cef_infection <- subset(metabolome_16s, abx == 'cefoperazone')
cef_infection$abx <- NULL
strep_infection <- subset(metabolome_16s, abx == 'streptomycin')
strep_infection$abx <- NULL
clinda_infection <- subset(metabolome_16s, abx == 'clindamycin')
clinda_infection$abx <- NULL

rm(all_infection, conv_infection, cef_infection, strep_infection, clinda_infection, abx_infection)

#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_file, width=12, height=6)
layout(matrix(c(1,2),
              nrow=1, ncol=2, byrow=TRUE))

#-------------------#

# OTUs

par(mar=c(3,4,1,1), las=1, mgp=c(2,0.75,0), xaxs='i', yaxs='i')
plot(x=otu_nmds$MDS1, y=otu_nmds$MDS2, xlim=c(-1.6,1.6), ylim=c(-1.5,1.5),
     xlab='NMDS axis 1', ylab='NMDS axis 2', xaxt='n', yaxt='n', pch=19, cex=0.2)
axis(side=1, at=seq(-1.6,1.6,0.4), labels=seq(-1.6,1.6,0.4))
axis(side=2, at=seq(-1.5,1.5,0.3), labels=seq(-1.5,1.5,0.3))
points(x=otu_cefoperazone$MDS1, y=otu_cefoperazone$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=otu_clindamycin$MDS1, y=otu_clindamycin$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=otu_streptomycin$MDS1, y=otu_streptomycin$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=otu_untreated$MDS1, y=otu_untreated$MDS2, bg='azure2', pch=21, cex=2, lwd=1.2)
legend('bottomright', legend=c('No Antibiotics','Streptomycin-pretreated','Cefoperzone-pretreated','Clindamycin-pretreated'), 
       pt.bg=c('azure2',strep_col,cef_col,clinda_col), 
       pch=21, cex=1.1, pt.cex=2.2, bty='n')
text(x=-0.73, y=1.35, labels='16S rRNA Gene Sequencing', cex=1.2)
mtext('a', side=2, line=2, las=2, adj=2, padj=-16, cex=1.3, font=2)

#-------------------#

# Metabolomics

par(mar=c(3,4,1,1), las=1, mgp=c(2,0.75,0), xaxs='i', yaxs='i')
plot(x=metabolome_nmds$MDS1, y=metabolome_nmds$MDS2, xlim=c(-0.3,0.3), ylim=c(-0.15,0.1),
     xlab='NMDS axis 1', ylab='NMDS axis 2', xaxt='n', yaxt='n', pch=19, cex=0.2)
axis(side=1, at=seq(-0.4,0.4,0.1), labels=seq(-0.4,0.4,0.1))
axis(side=2, at=seq(-0.3,0.3,0.1), labels=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3))
points(x=metabolome_cefoperazone$MDS1, y=metabolome_cefoperazone$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_clindamycin$MDS1, y=metabolome_clindamycin$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_streptomycin$MDS1, y=metabolome_streptomycin$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_germfree$MDS1, y=metabolome_germfree$MDS2, bg='gold1', pch=21, cex=2, lwd=1.2)
points(x=metabolome_untreated$MDS1, y=metabolome_untreated$MDS2, bg='azure2', pch=21, cex=2, lwd=1.2)
legend('bottomright', legend=c('No Antibiotics','Streptomycin-pretreated','Cefoperzone-pretreated','Clindamycin-pretreated'), 
       pt.bg=c('azure2',strep_col,cef_col,clinda_col), 
       pch=21, cex=1.1, pt.cex=2.2, bty='n')
text(x=-0.15, y=0.085, labels='Untargeted Metabolomics', cex=1.2)
mtext('b', side=2, line=2, las=2, adj=2, padj=-16, cex=1.3, font=2)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()

