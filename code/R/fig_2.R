# Set up environment
rm(list=ls())
gc()

# Load dependencies
deps <- c('vegan', 'plotrix', 'reshape2')
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Set seed for RNG
set.seed(9861)

# Load in functions
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Output plot name
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/figures/figure_2.pdf'

# Input 0.03 OTU shared file
shared_otu_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'

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
conv_metabolome <- metabolome[,colnames(metabolome) %in% c('ConvC1M1','ConvC1M2','ConvC1M3','ConvC1M4',
                                                       'ConvC2M1','ConvC2M2','ConvC2M3','ConvC2M4','ConvC2M5')] # Untreated SPF samples
metabolome <- metabolome[,!colnames(metabolome) %in% c('ConvC1M1','ConvC1M2','ConvC1M3','ConvC1M4',
                                                       'ConvC2M1','ConvC2M2','ConvC2M3','ConvC2M4','ConvC2M5')] # Untreated SPF samples
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
shared_otu <- shared_otu[,!(names(shared_otu) == 'Otu0004')] # Remove residual C. difficile OTU
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('StrepC4M1','StrepC4M2','StrepC4M3'), ]
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('GfC1M1','GfC1M2','GfC1M3',
                                                      'GfC2M1','GfC2M2','GfC2M3',
                                                      'GfC3M1','GfC3M2','GfC3M3',
                                                      'GfC4M1','GfC4M2','GfC4M3',
                                                      'GfC5M1','GfC5M2','GfC5M3',
                                                      'GfC6M1','GfC6M2','GfC6M3'), ] # Germfree samples
conv_otu <- shared_otu[rownames(shared_otu) %in% c('ConvC1M1','ConvC1M2','ConvC1M3','ConvC1M4',
                                                    'ConvC2M1','ConvC2M2','ConvC2M3','ConvC2M4','ConvC2M5'), ]
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('ConvC1M1','ConvC1M2','ConvC1M3','ConvC1M4',
                                                       'ConvC2M1','ConvC2M2','ConvC2M3','ConvC2M4','ConvC2M5'), ] # Untreated SPF samples
metabolome_16s <- clean_merge(metabolome, shared_otu)
metabolome_16s <- clean_merge(metadata, metabolome_16s)

#-------------------------------------------------------------------------------------------------------------------------#

# Calculate axes

# Metabolome
metabolome_nmds <- metaMDS(metabolome, k=2, trymax=100)$points
metabolome_nmds[,1] <- metabolome_nmds[,1] - 0.05
metabolome_nmds[,2] <- metabolome_nmds[,2] * -1
metabolome_nmds <- clean_merge(metadata, metabolome_nmds)

# Community structure
otu_nmds <- metaMDS(shared_otu, k=2, trymax=100)$points
otu_nmds[,1] <- otu_nmds[,1] - 0.2
otu_nmds[,2] <- otu_nmds[,2] + 0.1
otu_nmds <- clean_merge(metadata, otu_nmds)

# Subset to color points
metabolome_cefoperazone <- subset(metabolome_nmds, abx == 'cefoperazone')
metabolome_cefoperazone_630 <- subset(metabolome_cefoperazone, infection == '630')
metabolome_cefoperazone_mock <- subset(metabolome_cefoperazone, infection == 'mock')
rm(metabolome_cefoperazone)
metabolome_clindamycin <- subset(metabolome_nmds, abx == 'clindamycin')
metabolome_clindamycin_630 <- subset(metabolome_clindamycin, infection == '630')
metabolome_clindamycin_mock <- subset(metabolome_clindamycin, infection == 'mock')
rm(metabolome_clindamycin)
metabolome_streptomycin <- subset(metabolome_nmds, abx == 'streptomycin')
metabolome_streptomycin_630 <- subset(metabolome_streptomycin, infection == '630')
metabolome_streptomycin_mock <- subset(metabolome_streptomycin, infection == 'mock')
rm(metabolome_streptomycin)
otu_cefoperazone <- subset(otu_nmds, abx == 'cefoperazone')
otu_cefoperazone_630 <- subset(otu_cefoperazone, infection == '630')
otu_cefoperazone_mock <- subset(otu_cefoperazone, infection == 'mock')
rm(otu_cefoperazone)
otu_clindamycin <- subset(otu_nmds, abx == 'clindamycin')
otu_clindamycin_630 <- subset(otu_clindamycin, infection == '630')
otu_clindamycin_mock <- subset(otu_clindamycin, infection == 'mock')
rm(otu_clindamycin)
otu_streptomycin <- subset(otu_nmds, abx == 'streptomycin')
otu_streptomycin_630 <- subset(otu_streptomycin, infection == '630')
otu_streptomycin_mock <- subset(otu_streptomycin, infection == 'mock')
rm(otu_streptomycin)

#----------------#

# Combine 16S and metabolome for heatmaps

metabolome <- clean_merge(metadata, metabolome)
shared_otu <- clean_merge(metadata, shared_otu)




rm(metabolome, shared_otu)

#----------------#

# Calculate correlations for heatmap

#cormat <- round(cor(mydata),2)
#melted_cormat <- melt(cormat)


#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_file, width=8.5, height=11)
layout(matrix(c(1,2,
                3,3,
                3,3),
              nrow=3, ncol=2, byrow=TRUE))
par(mar=c(4,4,1,1), las=1, mgp=c(3,0.75,0), xaxs='i', yaxs='i')

#-------------------#

# OTUs alone
plot(x=otu_nmds$MDS1, y=otu_nmds$MDS2, xlim=c(-1.2,1.2), ylim=c(-1.5,1.5),
     xlab='NMDS axis 1', ylab='NMDS axis 2', xaxt='n', yaxt='n', pch=19, cex=0.2)
axis(side=1, at=seq(-1.2,1.2,0.4), labels=c(-1.2,-0.8,-0.4,0,0.4,0.8,1.2))
axis(side=2, at=seq(-1.5,1.5,0.5), labels=seq(-1.5,1.5,0.5))
mtext('a', side=2, line=2, las=2, adj=2, padj=-10, cex=1.2, font=2)
legend('topright', legend=c('16S rRNA Gene Sequencing'), pch=21, cex=1.4, pt.cex=0, bty='n')
legend('bottomright', legend=c('Streptomycin-pretreated','Cefoperzone-pretreated','Clindamycin-pretreated'), 
       pt.bg=c(strep_col,cef_col,clinda_col), 
       pch=22, cex=1.2, pt.cex=2.4)
legend('bottomleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(16,17), cex=1.2, pt.cex=c(2.2,2))

points(x=otu_cefoperazone_630$MDS1, y=otu_cefoperazone_630$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=otu_clindamycin_630$MDS1, y=otu_clindamycin_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=otu_streptomycin_630$MDS1, y=otu_streptomycin_630$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=otu_cefoperazone_mock$MDS1, y=otu_cefoperazone_mock$MDS2, bg=cef_col, pch=24, cex=1.8, lwd=1.2)
points(x=otu_clindamycin_mock$MDS1, y=otu_clindamycin_mock$MDS2, bg=clinda_col, pch=24, cex=1.8, lwd=1.2)
points(x=otu_streptomycin_mock$MDS1, y=otu_streptomycin_mock$MDS2, bg=strep_col, pch=24, cex=1.8, lwd=1.2)

#-------------------#

# Metabolomics alone
plot(x=metabolome_nmds$MDS1, y=metabolome_nmds$MDS2, xlim=c(-0.25,0.2), ylim=c(-0.25,0.2),
     xlab='NMDS axis 1', ylab='NMDS axis 2', xaxt='n', yaxt='n', pch=19, cex=0.2)
axis(side=1, at=seq(-0.25,0.2,0.05), labels=seq(-0.25,0.2,0.05))
axis(side=2, at=seq(-0.25,0.2,0.05), labels=seq(-0.25,0.2,0.05))
mtext('b', side=2, line=2, las=2, adj=2, padj=-10, cex=1.2, font=2)
legend('topright', legend=c('Untargeted Metabolomics'), pch=21, cex=1.4, pt.cex=0, bty='n')
legend('bottomright', legend=c('Streptomycin-pretreated','Cefoperzone-pretreated','Clindamycin-pretreated'), 
       pt.bg=c(strep_col,cef_col,clinda_col), 
       pch=22, cex=1.2, pt.cex=2.4)
legend('bottomleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(16,17), cex=1.2, pt.cex=c(2.2,2))

points(x=metabolome_cefoperazone_630$MDS1, y=metabolome_cefoperazone_630$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_clindamycin_630$MDS1, y=metabolome_clindamycin_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_streptomycin_630$MDS1, y=metabolome_streptomycin_630$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_cefoperazone_mock$MDS1, y=metabolome_cefoperazone_mock$MDS2, bg=cef_col, pch=24, cex=1.8, lwd=1.2)
points(x=metabolome_clindamycin_mock$MDS1, y=metabolome_clindamycin_mock$MDS2, bg=clinda_col, pch=24, cex=1.8, lwd=1.2)
points(x=metabolome_streptomycin_mock$MDS1, y=metabolome_streptomycin_mock$MDS2, bg=strep_col, pch=24, cex=1.8, lwd=1.2)

#-------------------#

par(mar=c(4,4,1,1), las=1, mgp=c(3,0.75,0), xaxs='i', yaxs='i')
plot(0, type='n', axes=FALSE, xlab='', ylab='')
mtext('c', side=2, line=2, las=2, adj=2, padj=-23, cex=1.2, font=2)

# Heatmap or correlation somehow...

dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up

for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()

