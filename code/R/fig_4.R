
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Metabolomes
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'

# Input Metadata
metadata <- 'data/metadata.tsv'

# Output files
plot_4abc <- 'results/figures/figure_4abc.pdf'
plot_4d <- 'results/figures/figure_4d.pdf'
plot_4e <- 'results/figures/figure_4e.pdf'
plot_4f <- 'results/figures/figure_4f.pdf'

#----------------#

# Read in data
# Metabolomes
metabolome <- read.delim(metabolome, sep='\t', header=TRUE)

# Metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data
# Metadata
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL

# Metabolomes
metabolome <- metabolome[,!colnames(metabolome) %in% c('GfC1M1','GfC1M2','GfC1M3', 
                                                       'GfC2M1','GfC2M2','GfC2M3', 
                                                       'GfC3M1','GfC3M2','GfC3M3', 
                                                       'GfC4M1','GfC4M2','GfC4M3', 
                                                       'GfC5M1','GfC5M2','GfC5M3', 
                                                       'GfC6M1','GfC6M2','GfC6M3')] # Germfree
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.data.frame(t(metabolome))
colnames(metabolome) <- gsub('_', ' ', colnames(metabolome))
substr(colnames(metabolome), 1, 1) <- toupper(substr(colnames(metabolome), 1, 1))

#-------------------------------------------------------------------------------------------------------------------------#

# Ordination
# Prep data
metabolome_subset <- clean_merge(metadata, metabolome)
cef_subset <- subset(metabolome_subset, abx == 'cefoperazone')
cef_subset$abx <- NULL
cef_subset$susceptibility <- NULL
clinda_subset <- subset(metabolome_subset, abx == 'clindamycin')
clinda_subset$abx <- NULL
clinda_subset$susceptibility <- NULL
strep_subset <- subset(metabolome_subset, abx == 'streptomycin')
strep_subset$abx <- NULL
strep_subset$susceptibility <- NULL

# Calculate significant differences
strep_p <- round(adonis(strep_subset[,2:ncol(strep_subset)]~factor(strep_subset$infection), data=strep_subset, permutations=10000, method='jaccard')$aov.tab[[6]][1], 3)
cef_p <- round(adonis(cef_subset[,2:ncol(cef_subset)]~factor(cef_subset$infection), data=cef_subset, permutations=10000, method='jaccard')$aov.tab[[6]][1], 3)
clinda_p <- round(adonis(clinda_subset[,2:ncol(clinda_subset)]~factor(clinda_subset$infection), data=clinda_subset, permutations=10000, method='jaccard')$aov.tab[[6]][1], 3)

# Calculate axes
# Separate analyses
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
                                                           'StrepC5M1','StrepC5M2','StrepC5M3',
                                                           'StrepC6M1','StrepC6M2','StrepC6M3'), ]

# Calculate axes and merge with metadata
cef_metabolome_nmds <- metaMDS(cef_metabolome, k=2, trymax=100)$points
cef_metabolome_nmds <- clean_merge(metadata, cef_metabolome_nmds)
clinda_metabolome_nmds <- metaMDS(clinda_metabolome, k=2, trymax=100)$points
clinda_metabolome_nmds <- clean_merge(metadata, clinda_metabolome_nmds)
strep_metabolome_nmds <- metaMDS(strep_metabolome, k=2, trymax=100)$points
strep_metabolome_nmds <- clean_merge(metadata, strep_metabolome_nmds)
rm(abx_metabolome, cef_metabolome, clinda_metabolome, strep_metabolome)

# Subset to points for plot
cef_metabolome_nmds$MDS2 <- cef_metabolome_nmds$MDS2 - 0.02
clinda_metabolome_nmds$MDS1 <- clinda_metabolome_nmds$MDS1 + 0.1
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

#-------------------------------------------------------------------------------------------------------------------------#

# Feature selection
# Separate groups
metabolome <- clean_merge(metadata, metabolome)
cef_metabolome <- subset(metabolome, abx == 'cefoperazone')
cef_metabolome$abx <- NULL
cef_metabolome$susceptibility <- NULL
cef_metabolome$infection <- factor(cef_metabolome$infection)
clinda_metabolome <- subset(metabolome, abx == 'clindamycin')
clinda_metabolome$abx <- NULL
clinda_metabolome$susceptibility <- NULL
clinda_metabolome$infection <- factor(clinda_metabolome$infection)
strep_metabolome <- subset(metabolome, abx == 'streptomycin')
strep_metabolome$abx <- NULL
strep_metabolome$susceptibility <- NULL
strep_metabolome$infection <- factor(strep_metabolome$infection)
rm(metadata, metabolome)

# Random Forest
cef_rf <- featureselect_RF(cef_metabolome, 'infection')
clinda_rf <- featureselect_RF(clinda_metabolome, 'infection')
strep_rf <- featureselect_RF(strep_metabolome, 'infection')

# Sort and subset top hits
cef_rf <- cef_rf[order(-cef_rf$MDA),][1:5,]
clinda_rf <- clinda_rf[order(-clinda_rf$MDA),][1:5,]
strep_rf <- strep_rf[order(-strep_rf$MDA),][1:5,]

# Subset concentrations
inf_cef_metabolome <- subset(cef_metabolome, infection == '630')[, cef_rf$feature]
inf_cef_metabolome$infection<- NULL
mock_cef_metabolome <- subset(cef_metabolome, infection == 'mock')[, cef_rf$feature]
mock_cef_metabolome$infection<- NULL
rm(cef_metabolome)

inf_clinda_metabolome <- subset(clinda_metabolome, infection == '630')[, clinda_rf$feature]
inf_clinda_metabolome$infection<- NULL
mock_clinda_metabolome <- subset(clinda_metabolome, infection == 'mock')[, clinda_rf$feature]
mock_clinda_metabolome$infection<- NULL
rm(clinda_metabolome)

inf_strep_metabolome <- subset(strep_metabolome, infection == '630')[, strep_rf$feature]
inf_strep_metabolome$infection<- NULL
mock_strep_metabolome <- subset(strep_metabolome, infection == 'mock')[, strep_rf$feature]
mock_strep_metabolome$infection<- NULL
rm(strep_metabolome)

# Find significant differences
cef_pvalues <- c()
for (i in 1:ncol(inf_cef_metabolome)){cef_pvalues[i] <- wilcox.test(inf_cef_metabolome[,i], mock_cef_metabolome[,i], exact=FALSE)$p.value}
cef_pvalues <- p.adjust(cef_pvalues, method='BH')
clinda_pvalues <- c()
for (i in 1:ncol(inf_clinda_metabolome)){clinda_pvalues[i] <- wilcox.test(inf_clinda_metabolome[,i], mock_clinda_metabolome[,i], exact=FALSE)$p.value}
pvalues <- p.adjust(clinda_pvalues, method='BH')
strep_pvalues <- c()
for (i in 1:ncol(inf_strep_metabolome)){strep_pvalues[i] <- wilcox.test(inf_strep_metabolome[,i], mock_strep_metabolome[,i], exact=FALSE)$p.value}
strep_pvalues <- p.adjust(strep_pvalues, method='BH')

#-------------------------------------------------------------------------------------------------------------------------#

# Generate figure panels
# Ordinations
pdf(file=plot_4abc, width=12, height=4)
layout(matrix(c(1,2,3),
              nrow=1, ncol=3, byrow=TRUE))
par(mar=c(4,4,1,1), las=1, mgp=c(2.8,0.75,0))
#Streptomycin - Fig. 4a
plot(x=strep_metabolome_nmds$MDS1, y=strep_metabolome_nmds$MDS2, xlim=c(-0.25,0.25), ylim=c(-0.15,0.15),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('a', side=2, line=2, las=2, adj=1.8, padj=-8.5, cex=1.6, font=2)
legend('topleft', legend='Streptomycin-pretreated', pch=1, cex=1.5, pt.cex=0, bty='n')
segments(x0=strep_metabolome_nmds_630$MDS1, y0=strep_metabolome_nmds_630$MDS2, x1=strep_metabolome_centoids[1,2], y1=strep_metabolome_centoids[1,3], col='gray30')
segments(x0=strep_metabolome_nmds_mock$MDS1, y0=strep_metabolome_nmds_mock$MDS2, x1=strep_metabolome_centoids[2,2], y1=strep_metabolome_centoids[2,3], col='gray30')
points(x=strep_metabolome_nmds_630$MDS1, y=strep_metabolome_nmds_630$MDS2, bg=strep_col, pch=22, cex=2, lwd=1.2)
points(x=strep_metabolome_nmds_mock$MDS1, y=strep_metabolome_nmds_mock$MDS2, bg=strep_col, pch=24, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Mock vs Infected', as.expression(bquote(paste(italic('p'),' = 0.039 *')))), 
       pch=1, cex=1.4, pt.cex=0, bty='n')
legend('bottomright', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(15,17), cex=1.2, pt.cex=2)
# Cefoperazone - Fig. 4b
plot(x=cef_metabolome_nmds$MDS1, y=cef_metabolome_nmds$MDS2, xlim=c(-0.15,0.15), ylim=c(-0.15,0.15),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('b', side=2, line=2, las=2, adj=1.8, padj=-8.5, cex=1.6, font=2)
legend('topleft', legend='Cefoperazone-pretreated', pch=1, cex=1.5, pt.cex=0, bty='n')
segments(x0=cef_metabolome_nmds_630$MDS1, y0=cef_metabolome_nmds_630$MDS2, x1=cef_metabolome_centoids[1,2], y1=cef_metabolome_centoids[1,3], col='gray30')
segments(x0=cef_metabolome_nmds_mock$MDS1, y0=cef_metabolome_nmds_mock$MDS2, x1=cef_metabolome_centoids[2,2], y1=cef_metabolome_centoids[2,3], col='gray30')
points(x=cef_metabolome_nmds_630$MDS1, y=cef_metabolome_nmds_630$MDS2, bg=cef_col, pch=22, cex=2, lwd=1.2)
points(x=cef_metabolome_nmds_mock$MDS1, y=cef_metabolome_nmds_mock$MDS2, bg=cef_col, pch=24, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Mock vs Infected', as.expression(bquote(paste(italic('p'),' = 0.016 *')))), 
       pch=1, cex=1.4, pt.cex=0, bty='n')
legend('bottomright', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(15,17), cex=1.2, pt.cex=2)
# Clindamycin - Fig. 4c
plot(x=clinda_metabolome_nmds$MDS1-0.075, y=clinda_metabolome_nmds$MDS2, xlim=c(-0.15,0.15), ylim=c(-0.1,0.1),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('c', side=2, line=2, las=2, adj=1.8, padj=-8.5, cex=1.6, font=2)
legend('topleft', legend='Clindamycin-pretreated', pch=1, cex=1.5, pt.cex=0, bty='n')
segments(x0=clinda_metabolome_nmds_630$MDS1-0.075, y0=clinda_metabolome_nmds_630$MDS2, x1=clinda_metabolome_centoids[1,2]-0.075, y1=clinda_metabolome_centoids[1,3], col='gray30')
segments(x0=clinda_metabolome_nmds_mock$MDS1-0.0751, y0=clinda_metabolome_nmds_mock$MDS2, x1=clinda_metabolome_centoids[2,2]-0.075, y1=clinda_metabolome_centoids[2,3], col='gray30')
points(x=clinda_metabolome_nmds_630$MDS1-0.075, y=clinda_metabolome_nmds_630$MDS2, bg=clinda_col, pch=22, cex=2, lwd=1.2)
points(x=clinda_metabolome_nmds_mock$MDS1-0.075, y=clinda_metabolome_nmds_mock$MDS2, bg=clinda_col, pch=24, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Mock vs Infected',as.expression(bquote(paste(italic('p'),' = 0.127 n.s.')))), 
       pch=1, cex=1.4, pt.cex=0, bty='n')
legend('bottomright', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(15,17), cex=1.2, pt.cex=2)
dev.off()

# Feature Selection
# Strep Infected - Fig. 4b
metabolite_stripchart(plot_4d, inf_strep_metabolome, mock_strep_metabolome, strep_pvalues, strep_rf$MDA, 0, 'Infected', 'Mock', 'Streptomycin-pretreated', strep_col, 'd')
# Cef Infected - Fig. 4c
metabolite_stripchart(plot_4e, inf_cef_metabolome, mock_cef_metabolome, cef_pvalues, cef_rf$MDA, 0, 'Infected', 'Mock', 'Cefoperazone-pretreated', cef_col, 'e')
# Clinda Infected - Fig. 4d
metabolite_stripchart(plot_4f, inf_clinda_metabolome, mock_clinda_metabolome, clinda_pvalues, clinda_rf$MDA, 33.33, 'Infected', 'Mock', 'Clindamycin-pretreated', clinda_col, 'f')

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()

