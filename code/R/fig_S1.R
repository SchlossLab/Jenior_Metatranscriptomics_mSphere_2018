
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
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/supplement/figures/figure_S1.pdf'

# Input 0.03 OTU shared file
shared_otu_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'

# Input Metadata
metadata_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metadata.tsv'

#----------------#

# Read in data

# 16S data
shared_otu <- read.delim(shared_otu_file, sep='\t', header=TRUE, row.names=2)
rm(shared_otu_file)

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

# 16S data
shared_otu$label <- NULL
shared_otu$numOtus <- NULL
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('CefC5M2'), ] # Contaminated sample
shared_otu <- shared_otu[,!(names(shared_otu) == 'Otu0004')] # Remove residual C. difficile OTU
#shared_otu <- shared_otu[!rownames(shared_otu) %in% c('StrepC4M1','StrepC4M2','StrepC4M3'), ]

#-------------------------------------------------------------------------------------------------------------------------#

# Subset antibiotic treatments
cef_otu <- shared_otu[rownames(shared_otu) %in% c('CefC1M1','CefC1M2','CefC1M3',
                                                  'CefC2M1','CefC2M2','CefC2M3',
                                                  'CefC3M1','CefC3M2','CefC3M3',
                                                  'CefC4M1','CefC4M2','CefC4M3',
                                                  'CefC5M1','CefC5M2','CefC5M3',
                                                  'CefC6M1','CefC6M2','CefC6M3'), ]
clinda_otu <- shared_otu[rownames(shared_otu) %in% c('ClindaC1M1','ClindaC1M2','ClindaC1M3',
                                                     'ClindaC2M1','ClindaC2M2','ClindaC2M3',
                                                     'ClindaC3M1','ClindaC3M2','ClindaC3M3',
                                                     'ClindaC4M1','ClindaC4M2','ClindaC4M3',
                                                     'ClindaC5M1','ClindaC5M2','ClindaC5M3',
                                                     'ClindaC6M1','ClindaC6M2','ClindaC6M3'), ]
strep_otu <- shared_otu[rownames(shared_otu) %in% c('StrepC1M1','StrepC1M2','StrepC1M3',
                                                    'StrepC2M1','StrepC2M2','StrepC2M3',
                                                    'StrepC3M1','StrepC3M2','StrepC3M3',
                                                    'StrepC4M1','StrepC4M2','StrepC4M3',
                                                    'StrepC5M1','ClindaC5M2','ClindaC5M3',
                                                    'StrepC6M1','StrepC6M2','StrepC6M3'), ]
rm(shared_otu)

#-------------------------------------------------------------------------------------------------------------------------#

# Calculate axes and merge with metadata
cef_nmds <- metaMDS(cef_otu, k=2, trymax=100)$points
cef_nmds <- clean_merge(metadata, cef_nmds)
clinda_nmds <- metaMDS(clinda_otu, k=2, trymax=100)$points
clinda_nmds <- clean_merge(metadata, clinda_nmds)
strep_nmds <- metaMDS(strep_otu, k=2, trymax=100)$points
strep_nmds <- clean_merge(metadata, strep_nmds)
rm(metadata)

# Calculate significant differences
cef_p <- as.character(anosim(cef_otu, cef_nmds$infection, permutations=999, distance='bray')$signif)
cef_r <- as.character(round(anosim(cef_otu, cef_nmds$infection, permutations=999, distance='bray')$statistic, 3))
clinda_p <- as.character(anosim(clinda_otu, clinda_nmds$infection, permutations=999, distance='bray')$signif)
clinda_r <- as.character(round(anosim(clinda_otu, clinda_nmds$infection, permutations=999, distance='bray')$statistic, 3))
strep_p <- as.character(anosim(strep_otu, strep_nmds$infection, permutations=999, distance='bray')$signif)
strep_r <- as.character(round(anosim(strep_otu, strep_nmds$infection, permutations=999, distance='bray')$statistic, 3))
rm(cef_otu, clinda_otu, strep_otu)

# Subset to points for plot
cef_nmds_630 <- subset(cef_nmds, infection == '630')
cef_nmds_mock <- subset(cef_nmds, infection == 'mock')
clinda_nmds_630 <- subset(clinda_nmds, infection == '630')
clinda_nmds_mock <- subset(clinda_nmds, infection == 'mock')
strep_nmds_630 <- subset(strep_nmds, infection == '630')
strep_nmds_mock <- subset(strep_nmds, infection == 'mock')

# Calculate centroids
cef_centoids <- aggregate(cbind(cef_nmds$MDS1,cef_nmds$MDS2)~cef_nmds$infection, data=cef_nmds, mean)
clinda_centoids <- aggregate(cbind(clinda_nmds$MDS1,clinda_nmds$MDS2)~clinda_nmds$infection, data=clinda_nmds, mean)
strep_centoids <- aggregate(cbind(strep_nmds$MDS1,strep_nmds$MDS2)~strep_nmds$infection, data=strep_nmds, mean)

#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_file, width=4, height=10)
layout(matrix(c(1,
                2,
                3),
              nrow=3, ncol=1, byrow=TRUE))
par(mar=c(4,4,1,1), las=1, mgp=c(2.5,0.75,0), xaxs='i', yaxs='i')

#-------------------#

# Streptomycin
plot(x=strep_nmds$MDS1, y=strep_nmds$MDS2, xlim=c(-1.2,1.2), ylim=c(-1.2,1.2),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('a', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.2, font=2)
legend('topright', legend='Streptomycin-pretreated', pch=1, cex=1.4, pt.cex=0, bty='n')
legend('bottomleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(16,17), cex=1.2, pt.cex=c(2.2,2))
segments(x0=strep_nmds_630$MDS1, y0=strep_nmds_630$MDS2, x1=strep_centoids[1,2], y1=strep_centoids[1,3], col='gray30')
segments(x0=strep_nmds_mock$MDS1, y0=strep_nmds_mock$MDS2, x1=strep_centoids[2,2], y1=strep_centoids[2,3], col='gray30')
points(x=strep_nmds_630$MDS1, y=strep_nmds_630$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=strep_nmds_mock$MDS1, y=strep_nmds_mock$MDS2, bg=strep_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c(as.expression(bquote(paste(italic('R'),' = ',.(strep_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(strep_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')

#-------------------#

# Cefoperazone
plot(x=cef_nmds$MDS1, y=cef_nmds$MDS2, xlim=c(-0.8,0.8), ylim=c(-0.8,0.8),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('b', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.2, font=2)
legend('topright', legend='Cefoperazone-pretreated', pch=1, cex=1.4, pt.cex=0, bty='n')
legend('bottomleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(16,17), cex=1.2, pt.cex=c(2.2,2))
segments(x0=cef_nmds_630$MDS1, y0=cef_nmds_630$MDS2, x1=cef_centoids[1,2], y1=cef_centoids[1,3], col='gray30')
segments(x0=cef_nmds_mock$MDS1, y0=cef_nmds_mock$MDS2, x1=cef_centoids[2,2], y1=cef_centoids[2,3], col='gray30')
points(x=cef_nmds_630$MDS1, y=cef_nmds_630$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=cef_nmds_mock$MDS1, y=cef_nmds_mock$MDS2, bg=cef_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c(as.expression(bquote(paste(italic('R'),' = ',.(cef_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(cef_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')

#-------------------#

# Clindamycin
plot(x=clinda_nmds$MDS1, y=clinda_nmds$MDS2, xlim=c(-0.9,0.9), ylim=c(-0.8,0.8),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('c', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.2, font=2)
legend('topright', legend='Clindamycin-pretreated', pch=1, cex=1.4, pt.cex=0, bty='n')
legend('bottomleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(16,17), cex=1.2, pt.cex=c(2.2,2))
segments(x0=clinda_nmds_630$MDS1, y0=clinda_nmds_630$MDS2, x1=clinda_centoids[1,2], y1=clinda_centoids[1,3], col='gray30')
segments(x0=clinda_nmds_mock$MDS1, y0=clinda_nmds_mock$MDS2, x1=clinda_centoids[2,2], y1=clinda_centoids[2,3], col='gray30')
points(x=clinda_nmds_630$MDS1, y=clinda_nmds_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=clinda_nmds_mock$MDS1, y=clinda_nmds_mock$MDS2, bg=clinda_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c(as.expression(bquote(paste(italic('R'),' = ',.(clinda_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(clinda_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')

dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()

