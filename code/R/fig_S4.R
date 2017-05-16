
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Output plot name
plot_file <- 'results/supplement/figures/figure_S4.pdf'

# Input Metabolomes
metabolome_file <- 'data/metabolome/metabolomics.tsv'

# Input Metadata
metadata_file <- 'data/metadata.tsv'

#----------------#

# Read in data

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

# Metabolomes
metabolome <- metabolome[,!colnames(metabolome) %in% c('CefC5M2')] # Contaminated sample
metabolome <- metabolome[,!colnames(metabolome) %in% c('GfC1M1','GfC1M2','GfC1M3', 
                                                       'GfC2M1','GfC2M2','GfC2M3', 
                                                       'GfC3M1','GfC3M2','GfC3M3', 
                                                       'GfC4M1','GfC4M2','GfC4M3', 
                                                       'GfC5M1','GfC5M2','GfC5M3', 
                                                       'GfC6M1','GfC6M2','GfC6M3')] # Germfree samples 
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome_annotation <- metabolome[,1:4]
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome <- t(metabolome)

#-------------------------------------------------------------------------------------------------------------------------#

# Calculate axes
# Metabolome - all
metabolome_nmds <- metaMDS(metabolome, k=2, trymax=100, distance='bray')$points
metabolome_nmds <- clean_merge(metadata, metabolome_nmds)
metabolome_nmds$MDS1 <- metabolome_nmds$MDS1 + 0.15

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
metabolome_noantibiotics <- subset(metabolome_nmds, abx == 'none')

#----------------#

# Separate plots
abx_metabolome <- metabolome[,!colnames(metabolome) %in% c('ConvC1M1','ConvC1M2','ConvC1M3','ConvC1M4',
                                                           'ConvC2M1','ConvC2M2','ConvC2M3','ConvC2M4','ConvC2M5')] # Untreated SPF samples 
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

# Calculate axes and merge with metadata
abx_metabolome_nmds <- metaMDS(abx_metabolome, k=2, trymax=100)$points
abx_metabolome_nmds <- clean_merge(metadata, abx_metabolome_nmds)
abx_metabolome_nmds$MDS1 <- abx_metabolome_nmds$MDS1 - 0.05
abx_metabolome_nmds$MDS2 <- abx_metabolome_nmds$MDS2 - 0.01
cef_metabolome_nmds <- metaMDS(cef_metabolome, k=2, trymax=100)$points
cef_metabolome_nmds <- clean_merge(metadata, cef_metabolome_nmds)
clinda_metabolome_nmds <- metaMDS(clinda_metabolome, k=2, trymax=100)$points
clinda_metabolome_nmds <- clean_merge(metadata, clinda_metabolome_nmds)
strep_metabolome_nmds <- metaMDS(strep_metabolome, k=2, trymax=100)$points
strep_metabolome_nmds <- clean_merge(metadata, strep_metabolome_nmds)
rm(metadata)

# Calculate significant differences
metabolome_p <- as.character(anosim(metabolome, metabolome_nmds$susceptibility, permutations=999, distance='bray')$signif)
metabolome_r <- as.character(round(anosim(metabolome, metabolome_nmds$susceptibility, permutations=999, distance='bray')$statistic, 3))
abx_metabolome_p <- as.character(anosim(abx_metabolome, abx_metabolome_nmds$abx, permutations=999, distance='bray')$signif)
abx_metabolome_r <- as.character(round(anosim(abx_metabolome, abx_metabolome_nmds$abx, permutations=999, distance='bray')$statistic, 3))
cef_metabolome_p <- as.character(anosim(cef_metabolome, cef_metabolome_nmds$infection, permutations=999, distance='bray')$signif)
cef_metabolome_r <- as.character(round(anosim(cef_metabolome, cef_metabolome_nmds$infection, permutations=999, distance='bray')$statistic, 3))
clinda_metabolome_p <- as.character(anosim(clinda_metabolome, clinda_metabolome_nmds$infection, permutations=999, distance='bray')$signif)
clinda_metabolome_r <- as.character(round(anosim(clinda_metabolome, clinda_metabolome_nmds$infection, permutations=999, distance='bray')$statistic, 3))
strep_metabolome_p <- as.character(anosim(strep_metabolome, strep_metabolome_nmds$infection, permutations=999, distance='bray')$signif)
strep_metabolome_r <- as.character(round(anosim(strep_metabolome, strep_metabolome_nmds$infection, permutations=999, distance='bray')$statistic, 3))
rm(cef_metabolome, clinda_metabolome, strep_metabolome)

# Subset to points for plot
metabolome_abx <- subset(abx_metabolome_nmds, abx == 'cefoperazone')
metabolome_abx_cef_630 <- subset(metabolome_abx, infection == '630')
metabolome_abx_cef_mock <- subset(metabolome_abx, infection == 'mock')
metabolome_abx <- subset(abx_metabolome_nmds, abx == 'clindamycin')
metabolome_abx_clinda_630 <- subset(metabolome_abx, infection == '630')
metabolome_abx_clinda_mock <- subset(metabolome_abx, infection == 'mock')
metabolome_abx <- subset(abx_metabolome_nmds, abx == 'streptomycin')
metabolome_abx_strep_630 <- subset(metabolome_abx, infection == '630')
metabolome_abx_strep_mock <- subset(metabolome_abx, infection == 'mock')
rm(metabolome_abx)
cef_metabolome_nmds$MDS2 <- cef_metabolome_nmds$MDS2 + 0.024
cef_metabolome_nmds_630 <- subset(cef_metabolome_nmds, infection == '630')
cef_metabolome_nmds_mock <- subset(cef_metabolome_nmds, infection == 'mock')
clinda_metabolome_nmds_630 <- subset(clinda_metabolome_nmds, infection == '630')
clinda_metabolome_nmds_mock <- subset(clinda_metabolome_nmds, infection == 'mock')
strep_metabolome_nmds_630 <- subset(strep_metabolome_nmds, infection == '630')
strep_metabolome_nmds_mock <- subset(strep_metabolome_nmds, infection == 'mock')

# Calculate centroids
cef_metabolome_centoids <- aggregate(cbind(cef_metabolome_nmds$MDS1,cef_metabolome_nmds$MDS2)~cef_metabolome_nmds$infection, data=cef_metabolome_nmds, mean)
clinda_metabolome_centoids <- aggregate(cbind(clinda_metabolome_nmds$MDS1,clinda_metabolome_nmds$MDS2)~clinda_metabolome_nmds$infection, data=clinda_metabolome_nmds, mean)
strep_metabolome_centoids <- aggregate(cbind(strep_metabolome_nmds$MDS1,strep_metabolome_nmds$MDS2)~strep_metabolome_nmds$infection, data=strep_metabolome_nmds, mean)

# sort metabolome tables equally
metabolome <- metabolome[match(rownames(metabolome_nmds), rownames(metabolome)),]

#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_file, width=12, height=8)
layout(matrix(c(1,2,3,
                4,5,6),
              nrow=2, ncol=3, byrow=TRUE))
par(mar=c(4,4,1,1), las=1, mgp=c(2.5,0.75,0))

#-------------------#

# All conventional mice
plot(x=metabolome_nmds$MDS1, y=metabolome_nmds$MDS2, xlim=c(-0.25,0.25), ylim=c(-0.15,0.15),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
legend('topleft', legend='All groups', pch=1, cex=1.4, pt.cex=0, bty='n')
mtext('a', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.4, font=2)
points(x=metabolome_cefoperazone_630$MDS1, y=metabolome_cefoperazone_630$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_clindamycin_630$MDS1, y=metabolome_clindamycin_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_streptomycin_630$MDS1, y=metabolome_streptomycin_630$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_cefoperazone_mock$MDS1, y=metabolome_cefoperazone_mock$MDS2, bg=cef_col, pch=24, cex=1.8, lwd=1.2)
points(x=metabolome_clindamycin_mock$MDS1, y=metabolome_clindamycin_mock$MDS2, bg=clinda_col, pch=24, cex=1.8, lwd=1.2)
points(x=metabolome_streptomycin_mock$MDS1, y=metabolome_streptomycin_mock$MDS2, bg=strep_col, pch=24, cex=1.8, lwd=1.2)
points(x=metabolome_noantibiotics$MDS1, y=metabolome_noantibiotics$MDS2, bg=noabx_col, pch=24, cex=1.8, lwd=1.2)
legend('bottomright', legend=c('Resistant vs Susceptible:',
                               as.expression(bquote(paste(italic('R'),' = ',.(metabolome_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(metabolome_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')

#-------------------------------------------------------------#

# Antibiotics only
plot(x=abx_metabolome_nmds$MDS1, y=abx_metabolome_nmds$MDS2, xlim=c(-0.075,0.075), ylim=c(-0.075,0.075),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
legend('topleft', legend='Antibiotic pretreatments only', pch=1, cex=1.4, pt.cex=0, bty='n')
mtext('b', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.4, font=2)
segments(x0=abx_metabolome_nmds_630$MDS1, y0=abx_metabolome_nmds_630$MDS2, x1=abx_metabolome_centoids[1,2], y1=abx_metabolome_centoids[1,3], col='gray30')
segments(x0=abx_metabolome_nmds_mock$MDS1, y0=abx_metabolome_nmds_mock$MDS2, x1=abx_metabolome_centoids[2,2], y1=abx_metabolome_centoids[2,3], col='gray30')
points(x=abx_metabolome_nmds_630$MDS1, y=abx_metabolome_nmds_630$MDS2, bg=abx_col, pch=21, cex=2, lwd=1.2)
points(x=abx_metabolome_nmds_mock$MDS1, y=abx_metabolome_nmds_mock$MDS2, bg=abx_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c('Strep vs Cef vs Clinda:',
                               as.expression(bquote(paste(italic('R'),' = ',.(abx_metabolome_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(abx_metabolome_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')
points(x=metabolome_abx_cef_630$MDS1, y=metabolome_abx_cef_630$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_abx_cef_mock$MDS1, y=metabolome_abx_cef_mock$MDS2, bg=cef_col, pch=24, cex=2, lwd=1.2)
points(x=metabolome_abx_clinda_630$MDS1, y=metabolome_abx_clinda_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_abx_clinda_mock$MDS1, y=metabolome_abx_clinda_mock$MDS2, bg=clinda_col, pch=24, cex=2, lwd=1.2)
points(x=metabolome_abx_strep_630$MDS1, y=metabolome_abx_strep_630$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_abx_strep_mock$MDS1, y=metabolome_abx_strep_mock$MDS2, bg=strep_col, pch=24, cex=2, lwd=1.2)

#-------------------------------------------------------------#

#Streptomycin
plot(x=strep_metabolome_nmds$MDS1, y=strep_metabolome_nmds$MDS2, xlim=c(-0.25,0.25), ylim=c(-0.25,0.25),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('c', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.4, font=2)
legend('topleft', legend='Streptomycin-pretreated', pch=1, cex=1.4, pt.cex=0, bty='n')
segments(x0=strep_metabolome_nmds_630$MDS1, y0=strep_metabolome_nmds_630$MDS2, x1=strep_metabolome_centoids[1,2], y1=strep_metabolome_centoids[1,3], col='gray30')
segments(x0=strep_metabolome_nmds_mock$MDS1, y0=strep_metabolome_nmds_mock$MDS2, x1=strep_metabolome_centoids[2,2], y1=strep_metabolome_centoids[2,3], col='gray30')
points(x=strep_metabolome_nmds_630$MDS1, y=strep_metabolome_nmds_630$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=strep_metabolome_nmds_mock$MDS1, y=strep_metabolome_nmds_mock$MDS2, bg=strep_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c('Mock vs Infected:',
                               as.expression(bquote(paste(italic('R'),' = ',.(strep_metabolome_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(strep_metabolome_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')

#-------------------#

# Cefoperazone
plot(x=cef_metabolome_nmds$MDS1, y=cef_metabolome_nmds$MDS2, xlim=c(-0.15,0.15), ylim=c(-0.1,0.1),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('d', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.4, font=2)
legend('topleft', legend='Cefoperazone-pretreated', pch=1, cex=1.4, pt.cex=0, bty='n')
segments(x0=cef_metabolome_nmds_630$MDS1, y0=cef_metabolome_nmds_630$MDS2, x1=cef_metabolome_centoids[1,2], y1=cef_metabolome_centoids[1,3], col='gray30')
segments(x0=cef_metabolome_nmds_mock$MDS1, y0=cef_metabolome_nmds_mock$MDS2, x1=cef_metabolome_centoids[2,2], y1=cef_metabolome_centoids[2,3], col='gray30')
points(x=cef_metabolome_nmds_630$MDS1, y=cef_metabolome_nmds_630$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=cef_metabolome_nmds_mock$MDS1, y=cef_metabolome_nmds_mock$MDS2, bg=cef_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c('Mock vs Infected:',
                               as.expression(bquote(paste(italic('R'),' = ',.(cef_metabolome_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(cef_metabolome_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')

#-------------------#

# Clindamycin
plot(x=clinda_metabolome_nmds$MDS1, y=clinda_metabolome_nmds$MDS2, xlim=c(-0.2,0.1), ylim=c(-0.1,0.1),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('e', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.4, font=2)
legend('topleft', legend='Clindamycin-pretreated', pch=1, cex=1.4, pt.cex=0, bty='n')
segments(x0=clinda_metabolome_nmds_630$MDS1, y0=clinda_metabolome_nmds_630$MDS2, x1=clinda_metabolome_centoids[1,2], y1=clinda_metabolome_centoids[1,3], col='gray30')
segments(x0=clinda_metabolome_nmds_mock$MDS1, y0=clinda_metabolome_nmds_mock$MDS2, x1=clinda_metabolome_centoids[2,2], y1=clinda_metabolome_centoids[2,3], col='gray30')
points(x=clinda_metabolome_nmds_630$MDS1, y=clinda_metabolome_nmds_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=clinda_metabolome_nmds_mock$MDS1, y=clinda_metabolome_nmds_mock$MDS2, bg=clinda_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c('Mock vs Infected:',
                               as.expression(bquote(paste(italic('R'),' = ',.(clinda_metabolome_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(clinda_metabolome_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')

#-------------------#

# Legends
par(mar=c(8,1,6,1))
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-5,5), ylim=c(-5,5))
legend('top', legend=c('Streptomycin-pretreated','Cefoperzone-pretreated','Clindamycin-pretreated','No Antibiotics'), 
       pt.bg=c(strep_col,cef_col,clinda_col, noabx_col), pch=22, cex=1.5, pt.cex=2.7)
legend('bottom', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(16,17), cex=1.5, pt.cex=c(2.5,2.2))

dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()
