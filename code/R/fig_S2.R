# Set up environment
rm(list=ls())
gc()

# Load in functions
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Output plot name
plot_file <- 'results/supplement/figures/figure_S2.pdf'

# Input 0.03 OTU shared file
shared_otu_file <- 'data/16S_analysis/all_treatments.0.03.unique_list.conventional.shared'

# Input Metabolomes
metabolome_file <- 'data/metabolome/metabolomics.tsv'

# Input Metadata
metadata_file <- 'data/metadata.tsv'

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
metabolome <- metabolome[,!colnames(metabolome) %in% c('ConvC1M1','ConvC1M2','ConvC1M3','ConvC1M4',
                                                       'ConvC2M1','ConvC2M2','ConvC2M3','ConvC2M4','ConvC2M5')] # Untreated SPF samples
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome_annotation <- metabolome[,1:4]
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome <- t(metabolome)

# 16S data
shared_otu$label <- NULL
shared_otu$numOtus <- NULL
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('CefC5M2'), ] # Contaminated sample
shared_otu <- shared_otu[,!(names(shared_otu) %in% c('Otu0004','Otu0308'))] # Remove residual C. difficile OTUs
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('ConvC1M1','ConvC1M2','ConvC1M3','ConvC1M4',
                                                       'ConvC2M1','ConvC2M2','ConvC2M3','ConvC2M4','ConvC2M5'), ] # Untreated SPF samples
sub_sample <- ceiling(min(rowSums(shared_otu)) * 0.9) # calculate rarefaction level
shared_otu <- rrarefy(shared_otu, sub_sample) # subsample shared file
shared_otu <- filter_table(shared_otu)
rm(sub_sample)

#-------------------------------------------------------------------------------------------------------------------------#

# Calculate axes

# Metabolome - all
metabolome_nmds <- metaMDS(metabolome, k=2, trymax=100, distance='bray')$points
metabolome_nmds[,1] <- metabolome_nmds[,1] - 0.05
metabolome_nmds[,2] <- metabolome_nmds[,2] * -1
metabolome_nmds[,2] <- metabolome_nmds[,2] + 0.05
metabolome_nmds <- clean_merge(metadata, metabolome_nmds)

# Community structure - all
otu_nmds <- metaMDS(shared_otu, k=2, trymax=100, distance='bray')$points
otu_nmds[,1] <- otu_nmds[,1] - 0.2
otu_nmds[,2] <- otu_nmds[,2] + 0.1
otu_nmds <- clean_merge(metadata, otu_nmds)

# Subset NMDS axes to color points - all
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

# Separate plots
cef_metabolome <- metabolome[rownames(metabolome) %in% c('CefC1M1','CefC1M2','CefC1M3',
                                                  'CefC2M1','CefC2M2','CefC2M3',
                                                  'CefC3M1','CefC3M2','CefC3M3',
                                                  'CefC4M1','CefC4M2','CefC4M3',
                                                  'CefC5M1','CefC5M2','CefC5M3',
                                                  'CefC6M1','CefC6M2','CefC6M3'), ]
clinda_metabolome <- metabolome[rownames(metabolome) %in% c('ClindaC1M1','ClindaC1M2','ClindaC1M3',
                                                     'ClindaC2M1','ClindaC2M2','ClindaC2M3',
                                                     'ClindaC3M1','ClindaC3M2','ClindaC3M3',
                                                     'ClindaC4M1','ClindaC4M2','ClindaC4M3',
                                                     'ClindaC5M1','ClindaC5M2','ClindaC5M3',
                                                     'ClindaC6M1','ClindaC6M2','ClindaC6M3'), ]
strep_metabolome <- metabolome[rownames(metabolome) %in% c('StrepC1M1','StrepC1M2','StrepC1M3',
                                                    'StrepC2M1','StrepC2M2','StrepC2M3',
                                                    'StrepC3M1','StrepC3M2','StrepC3M3',
                                                    'StrepC4M1','StrepC4M2','StrepC4M3',
                                                    'StrepC5M1','ClindaC5M2','ClindaC5M3',
                                                    'StrepC6M1','StrepC6M2','StrepC6M3'), ]
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
rm(metabolome, shared_otu)

# Filter out 0 sum columns
cef_otu <- filter_table(cef_otu)
clinda_otu <- filter_table(clinda_otu)
strep_otu <- filter_table(strep_otu)

# Calculate axes and merge with metadata
cef_otu_nmds <- metaMDS(cef_otu, k=2, trymax=100)$points
cef_otu_nmds <- clean_merge(metadata, cef_otu_nmds)
clinda_otu_nmds <- metaMDS(clinda_otu, k=2, trymax=100)$points
clinda_otu_nmds <- clean_merge(metadata, clinda_otu_nmds)
strep_otu_nmds <- metaMDS(strep_otu, k=2, trymax=100)$points
strep_otu_nmds <- clean_merge(metadata, strep_otu_nmds)
cef_metabolome_nmds <- metaMDS(cef_metabolome, k=2, trymax=100)$points
cef_metabolome_nmds <- clean_merge(metadata, cef_metabolome_nmds)
clinda_metabolome_nmds <- metaMDS(clinda_metabolome, k=2, trymax=100)$points
clinda_metabolome_nmds <- clean_merge(metadata, clinda_metabolome_nmds)
strep_metabolome_nmds <- metaMDS(strep_metabolome, k=2, trymax=100)$points
strep_metabolome_nmds <- clean_merge(metadata, strep_metabolome_nmds)
rm(metadata)


# Calculate significant differences
cef_otu_p <- as.character(anosim(cef_otu, cef_otu_nmds$infection, permutations=999, distance='bray')$signif)
cef_otu_r <- as.character(round(anosim(cef_otu, cef_otu_nmds$infection, permutations=999, distance='bray')$statistic, 3))
clinda_otu_p <- as.character(anosim(clinda_otu, clinda_otu_nmds$infection, permutations=999, distance='bray')$signif)
clinda_otu_r <- as.character(round(anosim(clinda_otu, clinda_otu_nmds$infection, permutations=999, distance='bray')$statistic, 3))
strep_otu_p <- as.character(anosim(strep_otu, strep_otu_nmds$infection, permutations=999, distance='bray')$signif)
strep_otu_r <- as.character(round(anosim(strep_otu, strep_otu_nmds$infection, permutations=999, distance='bray')$statistic, 3))
rm(cef_otu, clinda_otu, strep_otu)
cef_metabolome_p <- as.character(anosim(cef_metabolome, cef_metabolome_nmds$infection, permutations=999, distance='bray')$signif)
cef_metabolome_r <- as.character(round(anosim(cef_metabolome, cef_metabolome_nmds$infection, permutations=999, distance='bray')$statistic, 3))
clinda_metabolome_p <- as.character(anosim(clinda_metabolome, clinda_metabolome_nmds$infection, permutations=999, distance='bray')$signif)
clinda_metabolome_r <- as.character(round(anosim(clinda_metabolome, clinda_metabolome_nmds$infection, permutations=999, distance='bray')$statistic, 3))
strep_metabolome_p <- as.character(anosim(strep_metabolome, strep_metabolome_nmds$infection, permutations=999, distance='bray')$signif)
strep_metabolome_r <- as.character(round(anosim(strep_metabolome, strep_metabolome_nmds$infection, permutations=999, distance='bray')$statistic, 3))
rm(cef_metabolome, clinda_metabolome, strep_metabolome)

# Move points for easier viewing
strep_otu_nmds$MDS2 <- strep_otu_nmds$MDS2 + 0.2
cef_metabolome_nmds$MDS2 <- cef_metabolome_nmds$MDS2 + 0.025

# Subset to points for plot
cef_otu_nmds_630 <- subset(cef_otu_nmds, infection == '630')
cef_otu_nmds_mock <- subset(cef_otu_nmds, infection == 'mock')
clinda_otu_nmds_630 <- subset(clinda_otu_nmds, infection == '630')
clinda_otu_nmds_mock <- subset(clinda_otu_nmds, infection == 'mock')
strep_otu_nmds_630 <- subset(strep_otu_nmds, infection == '630')
strep_otu_nmds_mock <- subset(strep_otu_nmds, infection == 'mock')
cef_metabolome_nmds_630 <- subset(cef_metabolome_nmds, infection == '630')
cef_metabolome_nmds_mock <- subset(cef_metabolome_nmds, infection == 'mock')
clinda_metabolome_nmds_630 <- subset(clinda_metabolome_nmds, infection == '630')
clinda_metabolome_nmds_mock <- subset(clinda_metabolome_nmds, infection == 'mock')
strep_metabolome_nmds_630 <- subset(strep_metabolome_nmds, infection == '630')
strep_metabolome_nmds_mock <- subset(strep_metabolome_nmds, infection == 'mock')

# Calculate centroids
cef_otu_centoids <- aggregate(cbind(cef_otu_nmds$MDS1,cef_otu_nmds$MDS2)~cef_otu_nmds$infection, data=cef_otu_nmds, mean)
clinda_otu_centoids <- aggregate(cbind(clinda_otu_nmds$MDS1,clinda_otu_nmds$MDS2)~clinda_otu_nmds$infection, data=clinda_otu_nmds, mean)
strep_otu_centoids <- aggregate(cbind(strep_otu_nmds$MDS1,strep_otu_nmds$MDS2)~strep_otu_nmds$infection, data=strep_otu_nmds, mean)
cef_metabolome_centoids <- aggregate(cbind(cef_metabolome_nmds$MDS1,cef_metabolome_nmds$MDS2)~cef_metabolome_nmds$infection, data=cef_metabolome_nmds, mean)
clinda_metabolome_centoids <- aggregate(cbind(clinda_metabolome_nmds$MDS1,clinda_metabolome_nmds$MDS2)~clinda_metabolome_nmds$infection, data=clinda_metabolome_nmds, mean)
strep_metabolome_centoids <- aggregate(cbind(strep_metabolome_nmds$MDS1,strep_metabolome_nmds$MDS2)~strep_metabolome_nmds$infection, data=strep_metabolome_nmds, mean)

#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_file, width=8.5, height=14)
layout(matrix(c(1,2,
                3,4,
                5,6,
                7,8),
              nrow=4, ncol=2, byrow=TRUE))

#-------------------#

# OTUs all
par(mar=c(4,4,1,1), las=1, mgp=c(3,0.75,0), xaxs='i', yaxs='i')
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

# Metabolomics all
par(mar=c(4,4,1,1), las=1, mgp=c(3,0.75,0), xaxs='i', yaxs='i')
plot(x=metabolome_nmds$MDS1, y=metabolome_nmds$MDS2, xlim=c(-0.2,0.2), ylim=c(-0.2,0.2),
     xlab='NMDS axis 1', ylab='NMDS axis 2', xaxt='n', yaxt='n', pch=19, cex=0.2)
axis(side=1, at=seq(-0.25,0.2,0.05), labels=seq(-0.25,0.2,0.05))
axis(side=2, at=seq(-0.25,0.2,0.05), labels=seq(-0.25,0.2,0.05))
mtext('e', side=2, line=2, las=2, adj=2, padj=-10, cex=1.2, font=2)
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

#-------------------------------------------------------------#

# Streptomycin - otu
par(mar=c(4,4,1,1), las=1, mgp=c(2.5,0.75,0), xaxs='i', yaxs='i')
plot(x=strep_otu_nmds$MDS1, y=strep_otu_nmds$MDS2, xlim=c(-1.3,1.3), ylim=c(-1.8,1.8),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('b', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.2, font=2)
legend('topright', legend='Streptomycin-pretreated', pch=1, cex=1.4, pt.cex=0, bty='n')
legend('bottomleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(16,17), cex=1.2, pt.cex=c(2.2,2))
segments(x0=strep_otu_nmds_630$MDS1, y0=strep_otu_nmds_630$MDS2, x1=strep_otu_centoids[1,2], y1=strep_otu_centoids[1,3], col='gray30')
segments(x0=strep_otu_nmds_mock$MDS1, y0=strep_otu_nmds_mock$MDS2, x1=strep_otu_centoids[2,2], y1=strep_otu_centoids[2,3], col='gray30')
points(x=strep_otu_nmds_630$MDS1, y=strep_otu_nmds_630$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=strep_otu_nmds_mock$MDS1, y=strep_otu_nmds_mock$MDS2, bg=strep_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c(as.expression(bquote(paste(italic('R'),' = ',.(strep_otu_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(strep_otu_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')
legend('topleft', legend='16S', pch=1, pt.cex=0, bty='n')

# Streptomycin - metabolome
par(mar=c(4,4,1,1), las=1, mgp=c(2.5,0.75,0), xaxs='i', yaxs='i')
plot(x=strep_metabolome_nmds$MDS1, y=strep_metabolome_nmds$MDS2, xlim=c(-0.25,0.25), ylim=c(-0.25,0.25),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('f', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.2, font=2)
legend('topright', legend='Streptomycin-pretreated', pch=1, cex=1.4, pt.cex=0, bty='n')
legend('bottomleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(16,17), cex=1.2, pt.cex=c(2.2,2))
segments(x0=strep_metabolome_nmds_630$MDS1, y0=strep_metabolome_nmds_630$MDS2, x1=strep_metabolome_centoids[1,2], y1=strep_metabolome_centoids[1,3], col='gray30')
segments(x0=strep_metabolome_nmds_mock$MDS1, y0=strep_metabolome_nmds_mock$MDS2, x1=strep_metabolome_centoids[2,2], y1=strep_metabolome_centoids[2,3], col='gray30')
points(x=strep_metabolome_nmds_630$MDS1, y=strep_metabolome_nmds_630$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=strep_metabolome_nmds_mock$MDS1, y=strep_metabolome_nmds_mock$MDS2, bg=strep_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c(as.expression(bquote(paste(italic('R'),' = ',.(strep_metabolome_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(strep_metabolome_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')
legend('topleft', legend='Metabolome', pch=1, pt.cex=0, bty='n')

#-------------------#

# Cefoperazone - otu
par(mar=c(4,4,1,1), las=1, mgp=c(2.5,0.75,0), xaxs='i', yaxs='i')
plot(x=cef_otu_nmds$MDS1, y=cef_otu_nmds$MDS2, xlim=c(-1.2,1.2), ylim=c(-0.8,0.8),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('c', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.2, font=2)
legend('topright', legend='Cefoperazone-pretreated', pch=1, cex=1.4, pt.cex=0, bty='n')
legend('bottomleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(16,17), cex=1.2, pt.cex=c(2.2,2))
segments(x0=cef_otu_nmds_630$MDS1, y0=cef_otu_nmds_630$MDS2, x1=cef_otu_centoids[1,2], y1=cef_otu_centoids[1,3], col='gray30')
segments(x0=cef_otu_nmds_mock$MDS1, y0=cef_otu_nmds_mock$MDS2, x1=cef_otu_centoids[2,2], y1=cef_otu_centoids[2,3], col='gray30')
points(x=cef_otu_nmds_630$MDS1, y=cef_otu_nmds_630$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=cef_otu_nmds_mock$MDS1, y=cef_otu_nmds_mock$MDS2, bg=cef_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c(as.expression(bquote(paste(italic('R'),' = ',.(cef_otu_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(cef_otu_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')
legend('topleft', legend='16S', pch=1, pt.cex=0, bty='n')

# Cefoperazone - metabolome
par(mar=c(4,4,1,1), las=1, mgp=c(2.5,0.75,0), xaxs='i', yaxs='i')
plot(x=cef_metabolome_nmds$MDS1, y=cef_metabolome_nmds$MDS2, xlim=c(-0.15,0.15), ylim=c(-0.1,0.1),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('g', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.2, font=2)
legend('topright', legend='Cefoperazone-pretreated', pch=1, cex=1.4, pt.cex=0, bty='n')
legend('bottomleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(16,17), cex=1.2, pt.cex=c(2.2,2))
segments(x0=cef_metabolome_nmds_630$MDS1, y0=cef_metabolome_nmds_630$MDS2, x1=cef_metabolome_centoids[1,2], y1=cef_metabolome_centoids[1,3], col='gray30')
segments(x0=cef_metabolome_nmds_mock$MDS1, y0=cef_metabolome_nmds_mock$MDS2, x1=cef_metabolome_centoids[2,2], y1=cef_metabolome_centoids[2,3], col='gray30')
points(x=cef_metabolome_nmds_630$MDS1, y=cef_metabolome_nmds_630$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=cef_metabolome_nmds_mock$MDS1, y=cef_metabolome_nmds_mock$MDS2, bg=cef_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c(as.expression(bquote(paste(italic('R'),' = ',.(cef_metabolome_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(cef_metabolome_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')
legend('topleft', legend='Metabolome', pch=1, pt.cex=0, bty='n')

#-------------------#

# Clindamycin - otu
par(mar=c(4,4,1,1), las=1, mgp=c(2.5,0.75,0), xaxs='i', yaxs='i')
plot(x=clinda_otu_nmds$MDS1, y=clinda_otu_nmds$MDS2, xlim=c(-1.1,1.1), ylim=c(-1.2,1.2),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('d', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.2, font=2)
legend('topright', legend='Clindamycin-pretreated', pch=1, cex=1.4, pt.cex=0, bty='n')
legend('bottomleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(16,17), cex=1.2, pt.cex=c(2.2,2))
segments(x0=clinda_otu_nmds_630$MDS1, y0=clinda_otu_nmds_630$MDS2, x1=clinda_otu_centoids[1,2], y1=clinda_otu_centoids[1,3], col='gray30')
segments(x0=clinda_otu_nmds_mock$MDS1, y0=clinda_otu_nmds_mock$MDS2, x1=clinda_otu_centoids[2,2], y1=clinda_otu_centoids[2,3], col='gray30')
points(x=clinda_otu_nmds_630$MDS1, y=clinda_otu_nmds_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=clinda_otu_nmds_mock$MDS1, y=clinda_otu_nmds_mock$MDS2, bg=clinda_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c(as.expression(bquote(paste(italic('R'),' = ',.(clinda_otu_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(clinda_otu_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')
legend('topleft', legend='16S', pch=1, pt.cex=0, bty='n')

# Clindamycin - metabolome
par(mar=c(4,4,1,1), las=1, mgp=c(2.5,0.75,0), xaxs='i', yaxs='i')
plot(x=clinda_metabolome_nmds$MDS1, y=clinda_metabolome_nmds$MDS2, xlim=c(-0.2,0.1), ylim=c(-0.1,0.1),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('h', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.2, font=2)
legend('topright', legend='Clindamycin-pretreated', pch=1, cex=1.4, pt.cex=0, bty='n')
legend('bottomleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(16,17), cex=1.2, pt.cex=c(2.2,2))
segments(x0=clinda_metabolome_nmds_630$MDS1, y0=clinda_metabolome_nmds_630$MDS2, x1=clinda_metabolome_centoids[1,2], y1=clinda_metabolome_centoids[1,3], col='gray30')
segments(x0=clinda_metabolome_nmds_mock$MDS1, y0=clinda_metabolome_nmds_mock$MDS2, x1=clinda_metabolome_centoids[2,2], y1=clinda_metabolome_centoids[2,3], col='gray30')
points(x=clinda_metabolome_nmds_630$MDS1, y=clinda_metabolome_nmds_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=clinda_metabolome_nmds_mock$MDS1, y=clinda_metabolome_nmds_mock$MDS2, bg=clinda_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c(as.expression(bquote(paste(italic('R'),' = ',.(clinda_metabolome_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(clinda_metabolome_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')
legend('topleft', legend='Metabolome', pch=1, pt.cex=0, bty='n')

dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()

