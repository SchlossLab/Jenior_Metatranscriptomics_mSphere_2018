
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')

# Output plot name
plot_ab <- 'results/supplement/figures/figure_S1ac.pdf'
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
metabolome <- subset(metabolome, infection == 'mock')
metabolome$abx <- NULL
metabolome$infection <- NULL

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
shared_otu <- subset(shared_otu, infection == 'mock')
shared_otu$abx <- NULL
shared_otu$infection <- NULL
rm(otu_tax)

#-------------------------------------------------------------------------------------------------------------------------#

# Ordination analysis

# Metabolome - 728 metabolites
#metabolome_dist <- designdist(metabolome[,2:ncol(metabolome)], method='1-(J/(A+B-J))', terms='quadratic', abcd=FALSE) # Theta-YC
metabolome_dist <- vegdist(metabolome[,2:ncol(metabolome)], method='bray') # Bray-Curtis
metabolome_nmds <- as.data.frame(metaMDS(metabolome_dist, k=2, trymax=100)$points)
metabolome_nmds$MDS1 <- metabolome_nmds$MDS1 - 0.15
metabolome_nmds$MDS2 <- metabolome_nmds$MDS2 + 0.1
# Subset NMDS axes to color points
rownames(metabolome_nmds) <- rownames(metabolome)
metabolome_nmds <- clean_merge(metadata, metabolome_nmds)
metabolome_nmds_susceptible <- subset(metabolome_nmds, abx != 'none')
metabolome_nmds_strep <- subset(metabolome_nmds, abx == 'streptomycin')
metabolome_nmds_cef <- subset(metabolome_nmds, abx == 'cefoperazone')
metabolome_nmds_clinda <- subset(metabolome_nmds, abx == 'clindamycin')
metabolome_nmds_resistant <- subset(metabolome_nmds, abx == 'none')
# Calculate centroids
metabolome_res_centroids <- aggregate(cbind(metabolome_nmds_resistant$MDS1,metabolome_nmds_resistant$MDS2)~metabolome_nmds_resistant$susceptibility, data=metabolome_nmds_resistant, mean)
metabolome_sus_centroids <- aggregate(cbind(metabolome_nmds_susceptible$MDS1,metabolome_nmds_susceptible$MDS2)~metabolome_nmds_susceptible$susceptibility, data=metabolome_nmds_susceptible, mean)
# permANOVA
metabolome_PA_pval <- adonis(metabolome_dist ~ metabolome$susceptibility, metabolome, perm=999)$aov.tab
metabolome_PA_pval <- round(metabolome_pval[1,6], 3)
# Average within group distances
rm(metabolome_dist)

# 16S - 810 OTUs
#otu_dist <- designdist(shared_otu[,2:ncol(shared_otu)], method='1-(J/(A+B-J))', terms='quadratic', abcd=FALSE) # Theta-YC
otu_dist <- vegdist(shared_otu[,2:ncol(shared_otu)], method='bray') # Bray-Curtis
otu_nmds <- as.data.frame(metaMDS(otu_dist, k=2, trymax=100)$points)
otu_nmds$MDS1 <- otu_nmds$MDS1 - 0.1
# Subset NMDS axes to color points
rownames(otu_nmds) <- rownames(shared_otu)
otu_nmds <- clean_merge(metadata, otu_nmds)
otu_nmds_susceptible <- subset(otu_nmds, abx != 'none')
otu_nmds_strep <- subset(otu_nmds, abx == 'streptomycin')
otu_nmds_cef <- subset(otu_nmds, abx == 'cefoperazone')
otu_nmds_clinda <- subset(otu_nmds, abx == 'clindamycin')
otu_nmds_resistant <- subset(otu_nmds, abx == 'none')
# Calculate centroids
otu_res_centroids <- aggregate(cbind(otu_nmds_resistant$MDS1,otu_nmds_resistant$MDS2)~otu_nmds_resistant$infection, data=otu_nmds_resistant, mean)
otu_sus_centroids <- aggregate(cbind(otu_nmds_susceptible$MDS1,otu_nmds_susceptible$MDS2)~otu_nmds_susceptible$infection, data=otu_nmds_susceptible, mean)
# permANOVA
otu_PA_pval <- adonis(otu_dist ~ shared_otu$susceptibility, shared_otu, perm=999)$aov.tab
otu_PA_pval <- round(otu_pval[1,6], 3)
# Average within group distances
rm(otu_dist)

rm(metadata)

#-------------------------------------------------------------------------------------------------------------------------#

# Feature selection

# Metabolome
# AUCRF feature selection/reduction (to 0% OOB)
colnames(metabolome) <- make.names(colnames(metabolome))
metabolome_aucrf <- aucrfSusceptibility(metabolome)
# Get OOB
metabolome_aucrf_oob <- metabolome_aucrf$RFopt
metabolome_aucrf_oob <- metabolome_aucrf_oob$err.rate
metabolome_aucrf_oob <- as.character(round(median(metabolome_aucrf_oob[,1]) * 100, 3))
# Get features
metabolome_aucrf <- as.character(OptimalSet(metabolome_aucrf)$Name)
res_metabolome <- subset(metabolome, susceptibility == 'resistant')[, metabolome_aucrf]
res_metabolome$susceptibility <- NULL
sus_metabolome <- subset(metabolome, susceptibility == 'susceptible')[, metabolome_aucrf]
sus_metabolome$susceptibility <- NULL
rm(metabolome, metabolome_aucrf)
# Find significant differences
metabolome_pval <- c()
for (i in 1:ncol(res_metabolome)){metabolome_pval[i] <- wilcox.test(res_metabolome[,i], sus_metabolome[,i], exact=FALSE)$p.value}
metabolome_pval <- round(p.adjust(metabolome_pval, method='BH'), 4)

# 16S
# AUCRF feature selection/reduction (to 0% OOB)
colnames(shared_otu) <- make.names(colnames(shared_otu))
shared_otu_aucrf <- aucrfSusceptibility(shared_otu)
# Get OOB
shared_otu_aucrf_oob <- shared_otu_aucrf$RFopt
shared_otu_aucrf_oob <- shared_otu_aucrf_oob$err.rate
shared_otu_aucrf_oob <- as.character(round(median(shared_otu_aucrf_oob[,1]) * 100, 3))
# Get features
shared_otu_aucrf <- as.character(OptimalSet(shared_otu_aucrf)$Name)
res_shared_otu <- subset(shared_otu, susceptibility == 'resistant')[, shared_otu_aucrf]
res_shared_otu$susceptibility <- NULL
sus_shared_otu <- subset(shared_otu, susceptibility == 'susceptible')[, shared_otu_aucrf]
sus_shared_otu$susceptibility <- NULL
rm(shared_otu, shared_otu_aucrf)
# Find significant differences
shared_otu_pval <- c()
for (i in 1:ncol(res_shared_otu)){shared_otu_pval[i] <- wilcox.test(res_shared_otu[,i], sus_shared_otu[,i], exact=FALSE)$p.value}
shared_otu_pval <- round(p.adjust(shared_otu_pval, method='BH'), 4)
# Log transform values
res_shared_otu <- log10(res_shared_otu + 1)
sus_shared_otu <- log10(sus_shared_otu + 1)

#-------------#

# Reformat names to be more human readable
# Metabolome
colnames(res_metabolome) <- c("Nudifloramide","N-Acetylproline","Sebacate / Decanedioate","Hyodeoxycholate","Murideoxycholate")
colnames(sus_metabolome) <- c("Nudifloramide","N-Acetylproline","Sebacate / Decanedioate","Hyodeoxycholate","Murideoxycholate")
# 16S
# Top BLAST results
species <- c("Clostridium aerotolerans","Dialister succinatiphilus",
             "Clostridium saccharolyticum","Clostridium saccharolyticum",
             "Ruthenibacterium lactatiformans")
otu <- c("(OTU42)","(OTU102)","(OTU111)","(OTU115)","(OTU128)")
formatted_names <- lapply(1:length(species), function(i) bquote(paste(italic(.(species[i])), ' ', .(otu[i]), sep='')))
rm(species, otu)
# c("Lachnospiraceae (OTU42)","Clostridia unclassified (OTU102)","Lachnospiraceae (OTU111)","Clostridiales unclassified (OTU115)","Lachnospiraceae (OTU128)")

#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_ab, width=6.5, height=12)
layout(matrix(c(1,
                2),
              nrow=2, ncol=1, byrow=TRUE))
par(mar=c(4,4,1,3), las=1, mgp=c(2.8,0.75,0))

# 16S
plot(x=otu_nmds$MDS1, y=otu_nmds$MDS2, xlim=c(-0.6,0.6), ylim=c(-0.6,0.6),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('A', side=2, line=2, las=2, adj=1.5, padj=-10, cex=2, font=2)
segments(x0=otu_nmds_susceptible$MDS1, y0=otu_nmds_susceptible$MDS2, x1=otu_sus_centroids[1,2], y1=otu_sus_centroids[1,3], col='gray30')
points(x=otu_nmds_strep$MDS1, y=otu_nmds_strep$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=otu_nmds_cef$MDS1, y=otu_nmds_cef$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=otu_nmds_clinda$MDS1, y=otu_nmds_clinda$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
segments(x0=otu_nmds_resistant$MDS1, y0=otu_nmds_resistant$MDS2, x1=otu_res_centroids[1,2], y1=otu_res_centroids[1,3], col='gray30')
points(x=otu_nmds_resistant$MDS1, y=otu_nmds_resistant$MDS2, bg=noabx_col, pch=21, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Resistant vs Susceptible', as.expression(bquote(paste(italic('p'),' << 0.001')))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
legend('bottomright', legend=c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin'), 
       pt.bg=c(noabx_col, strep_col, cef_col, clinda_col), pch=21, cex=1.2, pt.cex=2.5)
legend('topleft', legend='Community Structure', bty='n', cex=1.4, pt.cex=0)
box()
mtext('B', side=4, line=2, las=2, adj=0.5, padj=-10, cex=2, font=2)

# Metabolome
plot(x=metabolome_nmds$MDS1, y=metabolome_nmds$MDS2, xlim=c(-0.6,0.6), ylim=c(-0.6,0.6),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('C', side=2, line=2, las=2, adj=1.5, padj=-10, cex=2, font=2)
segments(x0=metabolome_nmds_susceptible$MDS1, y0=metabolome_nmds_susceptible$MDS2, x1=metabolome_sus_centroids[1,2], y1=metabolome_sus_centroids[1,3], col='gray30')
points(x=metabolome_nmds_strep$MDS1, y=metabolome_nmds_strep$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_nmds_cef$MDS1, y=metabolome_nmds_cef$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_nmds_clinda$MDS1, y=metabolome_nmds_clinda$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
segments(x0=metabolome_nmds_resistant$MDS1, y0=metabolome_nmds_resistant$MDS2, x1=metabolome_res_centroids[1,2], y1=metabolome_res_centroids[1,3], col='gray30')
points(x=metabolome_nmds_resistant$MDS1, y=metabolome_nmds_resistant$MDS2, bg=noabx_col, pch=21, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Resistant vs Susceptible', as.expression(bquote(paste(italic('p'),' << 0.001')))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
legend('bottomright', legend=c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin'), 
       pt.bg=c(noabx_col, strep_col, cef_col, clinda_col), pch=21, cex=1.2, pt.cex=2.5)
legend('topleft', legend='Metabolome', bty='n', cex=1.4, pt.cex=0)
box()
mtext('D', side=4, line=2, las=2, adj=0.5, padj=-10, cex=2, font=2)

dev.off()

#---------------#

# Feature selection results
# 16S
multiStripchart(otu_plot, res_shared_otu, sus_shared_otu, shared_otu_pval, shared_otu_aucrf_oob, 
                'Resistant', 'Susceptible', noabx_col, 'forestgreen', '', 'black', 
                formatted_names, expression(paste('Relative Abundance (',log[10],')')))
# Metabolome
multiStripchart(metabolome_plot, res_metabolome, sus_metabolome, metabolome_pval, metabolome_aucrf_oob, 
                'Resistant', 'Susceptible', noabx_col, 'forestgreen', '', 'black', 
                as.list(colnames(res_metabolome)), expression(paste('Scaled Intensity (',log[10],')')))

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()
