
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')

# Output plot name
plot_ac <- 'results/supplement/figures/figure_S1ac.pdf'
metabolome_plot <- 'results/supplement/figures/figure_S1d.pdf'
otu_plot <- 'results/supplement/figures/figure_S1b.pdf'

# Metabolomes
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'

# 16S
shared_otu <- 'data/16S_analysis/all_treatments.0.03.unique_list.conventional.shared'
otu_tax <- 'data/16S_analysis/formatted.all_treatments.0.03.cons.taxonomy'

# Metadata
metadata <- 'data/metadata.tsv'

#----------------#

# Read in data

# Metabolomes
metabolome <- read.delim(metabolome, sep='\t', header=TRUE)

# 16S
shared_otu <- read.delim(shared_otu, sep='\t', header=T, row.names=2)
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('CefC5M2'), ]  # Remove possible contaminated sample
shared_otu <- shared_otu[,!(names(shared_otu) %in% c('Otu0004','Otu0308'))] # Remove residual C. difficile OTUs
shared_otu$numOtus <- NULL
shared_otu$label <- NULL
otu_tax <- read.delim(otu_tax, sep='\t', header=T, row.names=1)

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
metabolome$BIOCHEMICAL <- gsub('_', ' ', metabolome$BIOCHEMICAL)
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.data.frame(t(metabolome))
metabolome <- clean_merge(metadata, metabolome)
metabolome <- subset(metabolome, abx != 'germfree')
metabolome <- subset(metabolome, susceptibility == 'susceptible')
metabolome$abx <- NULL
metabolome$susceptibility <- NULL

# 16S
otu_tax$genus <- gsub('_', ' ', otu_tax$genus)
otu_tax$genus <- gsub('Ruminococcus2', 'Ruminococcus', otu_tax$genus)
otu_tax$taxon <- paste(otu_tax$genus, otu_tax$OTU.1, sep='_')
otu_tax$phylum <- NULL
otu_tax$genus <- NULL
otu_tax$OTU.1 <- NULL
shared_otu <- t(shared_otu)
shared_otu <- clean_merge(otu_tax, shared_otu)
rownames(shared_otu) <- shared_otu$taxon
shared_otu$taxon <- NULL
shared_otu <- t(shared_otu)
shared_otu <- clean_merge(metadata, shared_otu)
shared_otu <- subset(shared_otu, abx != 'germfree')
shared_otu <- subset(shared_otu, susceptibility == 'susceptible')
shared_otu$abx <- NULL
shared_otu$susceptibility <- NULL
rm(otu_tax)

#-------------------------------------------------------------------------------------------------------------------------#

# Ordination analysis

# Metabolome
#metabolome_dist <- designdist(metabolome[,2:ncol(metabolome)], method='1-(J/(A+B-J))', terms='quadratic', abcd=FALSE) # Theta-YC
metabolome_dist <- vegdist(metabolome[,2:ncol(metabolome)], method='bray') # Bray-Curtis
metabolome_nmds <- as.data.frame(metaMDS(metabolome_dist, k=2, trymax=100)$points)
# Subset NMDS axes to color points
rownames(metabolome_nmds) <- rownames(metabolome)
metabolome_nmds <- clean_merge(metadata, metabolome_nmds)
metabolome_nmds_mock <- subset(metabolome_nmds, infection == 'mock')
metabolome_nmds_infected <- subset(metabolome_nmds, infection == '630')
# Calculate centroids
metabolome_mock_centroids <- aggregate(cbind(metabolome_nmds_mock$MDS1,metabolome_nmds_mock$MDS2)~metabolome_nmds_mock$infection, data=metabolome_nmds_mock, mean)
metabolome_infected_centroids <- aggregate(cbind(metabolome_nmds_infected$MDS1,metabolome_nmds_infected$MDS2)~metabolome_nmds_infected$infection, data=metabolome_nmds_infected, mean)
# permANOVA
metabolome_permANOVA_pval <- adonis(metabolome_dist ~ metabolome$infection, metabolome, perm=999)$aov.tab
metabolome_permANOVA_pval <- round(metabolome_permANOVA_pval[1,6], 3)
# Average within group distances
rm(metabolome_dist)

# 16S
#otu_dist <- designdist(shared_otu[,2:ncol(shared_otu)], method='1-(J/(A+B-J))', terms='quadratic', abcd=FALSE) # Theta-YC
otu_dist <- vegdist(shared_otu[,2:ncol(shared_otu)], method='bray') # Bray-Curtis
otu_nmds <- as.data.frame(metaMDS(otu_dist, k=2, trymax=100)$points)
# Subset NMDS axes to color points
rownames(otu_nmds) <- rownames(shared_otu)
otu_nmds <- clean_merge(metadata, otu_nmds)
otu_nmds_mock <- subset(otu_nmds, infection == 'mock')
otu_nmds_infected <- subset(otu_nmds, infection == '630')
# Calculate centroids
otu_mock_centroids <- aggregate(cbind(otu_nmds_mock$MDS1,otu_nmds_mock$MDS2)~otu_nmds_mock$infection, data=otu_nmds_mock, mean)
otu_infected_centroids <- aggregate(cbind(otu_nmds_infected$MDS1,otu_nmds_infected$MDS2)~otu_nmds_infected$infection, data=otu_nmds_infected, mean)
# permANOVA
otu_permANOVA_pval <- adonis(otu_dist ~ shared_otu$infection, shared_otu, perm=999)$aov.tab
otu_permANOVA_pval <- round(otu_permANOVA_pval[1,6], 3)
# Average within group distances
rm(otu_dist)

rm(metadata)

#-------------------------------------------------------------------------------------------------------------------------#

# Feature selection

# Metabolome
# AUCRF feature selection/reduction (to 0% OOB)
colnames(metabolome) <- make.names(colnames(metabolome))
metabolome_aucrf <- aucrfInfection(metabolome)
# Get OOB
metabolome_aucrf_oob <- metabolome_aucrf$RFopt
metabolome_aucrf_oob <- metabolome_aucrf_oob$err.rate
metabolome_aucrf_oob <- as.character(round(median(metabolome_aucrf_oob[,1]) * 100, 2))
# Get features
metabolome_aucrf <- as.data.frame(OptimalSet(metabolome_aucrf))
print(nrow(metabolome_aucrf))
metabolome_aucrf <- metabolome_aucrf[order(-metabolome_aucrf$Importance), ][1:10,]
metabolome_aucrf <- as.character(metabolome_aucrf$Name)
mock_metabolome <- subset(metabolome, infection == 'mock')[, metabolome_aucrf]
mock_metabolome$infection <- NULL
infected_metabolome <- subset(metabolome, infection == '630')[, metabolome_aucrf]
infected_metabolome$sinfection <- NULL
rm(metabolome, metabolome_aucrf)
# Find significant differences
metabolome_pval <- c()
for (i in 1:ncol(mock_metabolome)){metabolome_pval[i] <- wilcox.test(mock_metabolome[,i], infected_metabolome[,i], exact=FALSE)$p.value}
metabolome_pval <- as.character(round(p.adjust(metabolome_pval, method='BH'), 4))

# 16S
# AUCRF feature selection/reduction (to 0% OOB)
colnames(shared_otu) <- make.names(colnames(shared_otu))
shared_otu_aucrf <- aucrfInfection(shared_otu)
# Get OOB
shared_otu_aucrf_oob <- shared_otu_aucrf$RFopt
shared_otu_aucrf_oob <- shared_otu_aucrf_oob$err.rate
shared_otu_aucrf_oob <- as.character(round(median(shared_otu_aucrf_oob[,1]) * 100, 2))
# Get features
shared_otu_aucrf <- as.data.frame(OptimalSet(shared_otu_aucrf))
print(nrow(shared_otu_aucrf))
shared_otu_aucrf <- shared_otu_aucrf[order(-shared_otu_aucrf$Importance), ][1:10,]
shared_otu_aucrf <- as.character(shared_otu_aucrf$Name)
mock_shared_otu <- subset(shared_otu, infection == 'mock')[, shared_otu_aucrf]
mock_shared_otu$infection <- NULL
infected_shared_otu <- subset(shared_otu, infection == '630')[, shared_otu_aucrf]
infected_shared_otu$sinfection <- NULL
rm(shared_otu, shared_otu_aucrf)
# Find significant differences
shared_otu_pval <- c()
for (i in 1:ncol(mock_shared_otu)){shared_otu_pval[i] <- wilcox.test(mock_shared_otu[,i], infected_shared_otu[,i], exact=FALSE)$p.value}
shared_otu_pval <- as.character(round(p.adjust(shared_otu_pval, method='BH'), 4))
# Log transform values
mock_shared_otu <- log10(mock_shared_otu + 1)
infected_shared_otu <- log10(infected_shared_otu + 1)

#-------------#

# Reformat names to be more human readable
# Metabolome
colnames(mock_metabolome) <- c("5-aminovalerate","proline","N2-dimethylguanine","trans-4-hydroxyproline","pro-hydroxy-pro",
                               "arabonate/xylonate","thioproline","N-acetylthreonine","N-acetylserine","adenine")
colnames(infected_metabolome) <- c("5-aminovalerate","proline","N2-dimethylguanine","trans-4-hydroxyproline","pro-hydroxy-pro",
                                   "arabonate/xylonate","thioproline","N-acetylthreonine","N-acetylserine","adenine")
# 16S
colnames(mock_shared_otu) <- c("Porphyromonadaceae (OTU32)","Prevotella (OTU24)","Olsenella (OTU20)","Enterococcus (OTU68)","Porphyromonadaceae (OTU5)",
                               "Barnesiella (OTU147)","Pyramidobacter (OTU145)","Bacteroides (OTU3)","Fusobacterium (OTU100)","Parabacteroides (OTU69)")
colnames(infected_shared_otu) <- c("Porphyromonadaceae (OTU32)","Prevotella (OTU24)","Olsenella (OTU20)","Enterococcus (OTU68)","Porphyromonadaceae (OTU5)",
                                   "Barnesiella (OTU147)","Pyramidobacter (OTU145)","Bacteroides (OTU3)","Fusobacterium (OTU100)","Parabacteroides (OTU69)")

#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_ac, width=12, height=6)
layout(matrix(c(1,2),
              nrow=1, ncol=2, byrow=TRUE))
par(mar=c(5,4,1,1), las=1, mgp=c(2.8,0.75,0))

# 16S
plot(x=otu_nmds$MDS1, y=otu_nmds$MDS2, xlim=c(-0.6,0.6), ylim=c(-0.6,0.6),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('A', side=2, line=2, las=2, adj=1.5, padj=-9, cex=2, font=2)
mtext('B', side=2, line=2, las=2, adj=1.5, padj=14, cex=2, font=2)
segments(x0=otu_nmds_infected$MDS1, y0=otu_nmds_infected$MDS2, x1=otu_infected_centroids[1,2], y1=otu_infected_centroids[1,3], col='gray30')
points(x=otu_nmds_infected$MDS1, y=otu_nmds_infected$MDS2, bg='chocolate2', pch=21, cex=2, lwd=1.2)
segments(x0=otu_nmds_mock$MDS1, y0=otu_nmds_mock$MDS2, x1=otu_mock_centroids[1,2], y1=otu_mock_centroids[1,3], col='gray30')
points(x=otu_nmds_mock$MDS1, y=otu_nmds_mock$MDS2, bg='darkblue', pch=21, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Mock vs Infected', as.expression(bquote(paste(italic('p'),' = 0.191')))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
legend('bottomright', legend=c('Mock-infected','C. difficile-infected'), 
       pt.bg=c('darkblue', 'chocolate2'), pch=21, cex=1.2, pt.cex=2.5)
legend('topleft', legend='Community Structure', bty='n', cex=1.4, pt.cex=0)
box()

# Metabolome
plot(x=metabolome_nmds$MDS1, y=metabolome_nmds$MDS2, xlim=c(-0.6,0.6), ylim=c(-0.6,0.6),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('C', side=2, line=2, las=2, adj=1.5, padj=-9, cex=2, font=2)
mtext('D', side=2, line=2, las=2, adj=1.5, padj=14, cex=2, font=2)
segments(x0=metabolome_nmds_infected$MDS1, y0=metabolome_nmds_infected$MDS2, x1=metabolome_infected_centroids[1,2], y1=metabolome_infected_centroids[1,3], col='gray30')
points(x=metabolome_nmds_infected$MDS1, y=metabolome_nmds_infected$MDS2, bg='chocolate2', pch=21, cex=2, lwd=1.2)
segments(x0=metabolome_nmds_mock$MDS1, y0=metabolome_nmds_mock$MDS2, x1=metabolome_mock_centroids[1,2], y1=metabolome_mock_centroids[1,3], col='gray30')
points(x=metabolome_nmds_mock$MDS1, y=metabolome_nmds_mock$MDS2, bg='darkblue', pch=21, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Mock vs Infected', as.expression(bquote(paste(italic('p'),' = 0.08')))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
legend('bottomright', legend=c('Mock-infected','C. difficile-infected'), 
       pt.bg=c('darkblue', 'chocolate2'), pch=21, cex=1.2, pt.cex=2.5)
legend('topleft', legend='Metabolome', bty='n', cex=1.4, pt.cex=0)
box()

dev.off()

#---------------#

# Feature selection results
# 16S
multiStripchart(otu_plot, mock_shared_otu, infected_shared_otu, shared_otu_pval, shared_otu_aucrf_oob, 
                'Mock-infected', 'C. difficile-infected', 'darkblue', 'chocolate2', '', 'black', 
                as.list(colnames(mock_shared_otu)), expression(paste('Relative Abundance (',log[10],')')))
# Metabolome
multiStripchart(metabolome_plot, mock_metabolome, infected_metabolome, metabolome_pval, metabolome_aucrf_oob, 
                'Mock-infected', 'C. difficile-infected', 'darkblue', 'chocolate2', '', 'black',
                as.list(colnames(mock_metabolome)), expression(paste('Scaled Intensity (',log[10],')')))

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()
