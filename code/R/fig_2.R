
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')

# Output plot name
plot_ab <- 'results/figures/figure_2ac.pdf'
metabolome_plot <- 'results/figures/figure_2d.pdf'
otu_plot <- 'results/figures/figure_2b.pdf'

# Metabolomes
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'

# 16S
shared_otu <- 'data/16S_analysis/all_treatments.0.03.unique_list.conventional.shared'
otu_tax <- 'data/16S_analysis/formatted.all_treatments.0.03.cons.taxonomy'

# Metadata
metadata <- 'data/metadata.tsv'

aucrfClearance <- function(training_data){
  # Format levels of susceptibility for AUCRF
  colnames(training_data) <- make.names(colnames(training_data))
  levels <- as.vector(unique(training_data$clearance))
  training_data$clearance <- as.character(training_data$clearance)
  training_data$clearance[which(training_data$clearance==levels[1])] <- 0
  training_data$clearance[which(training_data$clearance==levels[2])] <- 1
  training_data$clearance <- as.factor(as.numeric(training_data$clearance))
  rm(levels)
  # Run AUCRF with reproduceable parameters
  set.seed(906801)
  data_RF <- AUCRF(clearance ~ ., data=training_data, pdel=0.05, k0=5, ranking='MDA')
  return(data_RF)
}

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
metadata$clearance <- c(rep('colonized',18), rep('cleared',18), rep('colonized',18),
                        rep('mock',12), rep('colonized',18))

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
metabolome <- subset(metabolome, abx != 'none')
metabolome <- subset(metabolome, infection == 'mock')
metabolome$abx <- NULL
metabolome$susceptibility <- NULL
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
shared_otu <- subset(shared_otu, abx != 'none')
shared_otu <- subset(shared_otu, infection == 'mock')
shared_otu$abx <- NULL
shared_otu$susceptibility <- NULL
shared_otu$infection <- NULL
rm(otu_tax)

#-------------------------------------------------------------------------------------------------------------------------#

# Ordination analysis

# Metabolome - 728 metabolites
#metabolome_dist <- designdist(metabolome[,2:ncol(metabolome)], method='1-(J/(A+B-J))', terms='quadratic', abcd=FALSE) # Theta-YC
metabolome_dist <- vegdist(metabolome[,2:ncol(metabolome)], method='bray') # Bray-Curtis
metabolome_nmds <- as.data.frame(metaMDS(metabolome_dist, k=2, trymax=100)$points)
# Subset NMDS axes to color points
rownames(metabolome_nmds) <- rownames(metabolome)
metabolome_nmds <- clean_merge(metadata, metabolome_nmds)
metabolome_nmds_strep <- subset(metabolome_nmds, abx == 'streptomycin')
metabolome_nmds_cef <- subset(metabolome_nmds, abx == 'cefoperazone')
metabolome_nmds_clinda <- subset(metabolome_nmds, abx == 'clindamycin')
metabolome_nmds_cleared <- subset(metabolome_nmds, clearance == 'cleared')
metabolome_nmds_colonized <- subset(metabolome_nmds, clearance == 'colonized')
# Calculate centroids
metabolome_cleared_centroids <- aggregate(cbind(metabolome_nmds_cleared$MDS1, metabolome_nmds_cleared$MDS2)~metabolome_nmds_cleared$clearance, data=metabolome_nmds_cleared, mean)
metabolome_colonized_centroids <- aggregate(cbind(metabolome_nmds_colonized$MDS1, metabolome_nmds_colonized$MDS2)~metabolome_nmds_colonized$clearance, data=metabolome_nmds_colonized, mean)
# permANOVA
metabolome_permANOVA_pval <- adonis(metabolome_dist ~ metabolome$clearance, metabolome, perm=999)$aov.tab
metabolome_permANOVA_pval <- round(metabolome_permANOVA_pval[1,6], 3)
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
otu_nmds_strep <- subset(otu_nmds, abx == 'streptomycin')
otu_nmds_cef <- subset(otu_nmds, abx == 'cefoperazone')
otu_nmds_clinda <- subset(otu_nmds, abx == 'clindamycin')
otu_nmds_cleared <- subset(otu_nmds, clearance == 'cleared')
otu_nmds_colonized <- subset(otu_nmds, clearance == 'colonized')
# Calculate centroids
otu_cleared_centroids <- aggregate(cbind(otu_nmds_cleared$MDS1, otu_nmds_cleared$MDS2)~otu_nmds_cleared$clearance, data=otu_nmds_cleared, mean)
otu_colonized_centroids <- aggregate(cbind(otu_nmds_colonized$MDS1, otu_nmds_colonized$MDS2)~otu_nmds_colonized$clearance, data=otu_nmds_colonized, mean)
# permANOVA
otu_permANOVA_pval <- adonis(otu_dist ~ shared_otu$clearance, shared_otu, perm=999)$aov.tab
otu_permANOVA_pval <- round(otu_permANOVA_pval[1,6], 3)
# Average within group distances
rm(otu_dist)

rm(metadata)

#-------------------------------------------------------------------------------------------------------------------------#

# Feature selection

# Metabolome
# AUCRF feature selection/reduction (to 0% OOB)
colnames(metabolome) <- make.names(colnames(metabolome))
metabolome_aucrf <- aucrfClearance(metabolome)
# Get OOB
metabolome_aucrf_oob <- metabolome_aucrf$RFopt
metabolome_aucrf_oob <- metabolome_aucrf_oob$err.rate
metabolome_aucrf_oob <- as.character(round(median(metabolome_aucrf_oob[,1]) * 100, 3))
# Get features
metabolome_aucrf <- as.character(OptimalSet(metabolome_aucrf)$Name)
cleared_metabolome <- subset(metabolome, clearance == 'cleared')[, metabolome_aucrf]
cleared_metabolome$clearance <- NULL
colonized_metabolome <- subset(metabolome, clearance == 'colonized')[, metabolome_aucrf]
colonized_metabolome$clearance <- NULL
rm(metabolome, metabolome_aucrf)
# Find significant differences
metabolome_pval <- c()
for (i in 1:ncol(cleared_metabolome)){metabolome_pval[i] <- wilcox.test(cleared_metabolome[,i], colonized_metabolome[,i], exact=FALSE)$p.value}
metabolome_pval <- round(p.adjust(metabolome_pval, method='BH'), 4)

# 16S
# AUCRF feature selection/reduction (to 0% OOB)
colnames(shared_otu) <- make.names(colnames(shared_otu))
shared_otu_aucrf <- aucrfClearance(shared_otu)
# Get OOB
shared_otu_aucrf_oob <- shared_otu_aucrf$RFopt
shared_otu_aucrf_oob <- shared_otu_aucrf_oob$err.rate
shared_otu_aucrf_oob <- as.character(round(median(shared_otu_aucrf_oob[,1]) * 100, 3))
# Get features
shared_otu_aucrf <- as.character(OptimalSet(shared_otu_aucrf)$Name)
cleared_shared_otu <- subset(shared_otu, clearance == 'cleared')[, shared_otu_aucrf]
cleared_shared_otu$clearance <- NULL
colonized_shared_otu <- subset(shared_otu, clearance == 'colonized')[, shared_otu_aucrf]
colonized_shared_otu$clearance <- NULL
rm(shared_otu, shared_otu_aucrf)
# Find significant differences
shared_otu_pval <- c()
for (i in 1:ncol(cleared_shared_otu)){shared_otu_pval[i] <- wilcox.test(cleared_shared_otu[,i], colonized_shared_otu[,i], exact=FALSE)$p.value}
shared_otu_pval <- round(p.adjust(shared_otu_pval, method='BH'), 4)
# Log transform values
cleared_shared_otu <- log10(cleared_shared_otu + 1)
colonized_shared_otu <- log10(colonized_shared_otu + 1)

#-------------#

# Reformat names to be more human readable
# Metabolome
colnames(cleared_metabolome) <- c("cis-4-Decenoylcarnitine","Sucrose","4-Imidazoleacetate","Nicotinamide-ribonucleotide","Palmitoyl-Dihydrosphingomyelin")
colnames(colonized_metabolome) <- c("cis-4-Decenoylcarnitine","Sucrose","4-Imidazoleacetate","Nicotinamide-ribonucleotide","Palmitoyl-Dihydrosphingomyelin")
# 16S
# Top BLAST results
otu <- c("(OTU2)","(OTU6)","(OTU17)","(OTU3)","(OTU16)")
# "Escherichia.Shigella_.OTU2." "Turicibacter_.OTU6." "Arthrobacter_.OTU17." "Bacteroides_.OTU3." "Lactobacillus_.OTU16."  
species <- c("Escherichia coli",
             "Turicibacter sanguinis",
             "Arthrobacter spp.",
             "Bacteroides caecimuris",
             "Lactobacillus spp.")
formatted_names <- lapply(1:length(species), function(i) bquote(paste(italic(.(species[i])), ' ', .(otu[i]), sep='')))
rm(species, otu)

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
segments(x0=otu_nmds_colonized$MDS1, y0=otu_nmds_colonized$MDS2, x1=otu_colonized_centroids[1,2], y1=otu_colonized_centroids[1,3], col='gray30')
points(x=otu_nmds_strep$MDS1, y=otu_nmds_strep$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=otu_nmds_cef$MDS1, y=otu_nmds_cef$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
segments(x0=otu_nmds_cleared$MDS1, y0=otu_nmds_cleared$MDS2, x1=otu_cleared_centroids[1,2], y1=otu_cleared_centroids[1,3], col='gray30')
points(x=otu_nmds_clinda$MDS1, y=otu_nmds_clinda$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Cleared vs Colonized', as.expression(bquote(paste(italic('p'),' = 0.001')))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
legend('bottomright', legend=c('Streptomycin','Cefoperazone','Clindamycin'), 
       pt.bg=c(strep_col, cef_col, clinda_col), pch=21, cex=1.2, pt.cex=2.5)
legend('topleft', legend='Community Structure', bty='n', cex=1.4, pt.cex=0)
box()
mtext('B', side=4, line=2, las=2, adj=0.5, padj=-10, cex=2, font=2)

# Metabolome
plot(x=metabolome_nmds$MDS1, y=metabolome_nmds$MDS2, xlim=c(-0.5,0.5), ylim=c(-0.5,0.5),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('C', side=2, line=2, las=2, adj=1.5, padj=-10, cex=2, font=2)
segments(x0=metabolome_nmds_colonized$MDS1, y0=metabolome_nmds_colonized$MDS2, x1=metabolome_colonized_centroids[1,2], y1=metabolome_colonized_centroids[1,3], col='gray30')
points(x=metabolome_nmds_strep$MDS1, y=metabolome_nmds_strep$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_nmds_cef$MDS1, y=metabolome_nmds_cef$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
segments(x0=metabolome_nmds_cleared$MDS1, y0=metabolome_nmds_cleared$MDS2, x1=metabolome_cleared_centroids[1,2], y1=metabolome_cleared_centroids[1,3], col='gray30')
points(x=metabolome_nmds_clinda$MDS1, y=metabolome_nmds_clinda$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Cleared vs Colonized', as.expression(bquote(paste(italic('p'),' << 0.001')))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
legend('bottomright', legend=c('Streptomycin','Cefoperazone','Clindamycin'), 
       pt.bg=c(strep_col, cef_col, clinda_col), pch=21, cex=1.2, pt.cex=2.5)
legend('topleft', legend='Metabolome', bty='n', cex=1.4, pt.cex=0)
box()
mtext('D', side=4, line=2, las=2, adj=0.5, padj=-10, cex=2, font=2)

dev.off()

#---------------#

# Feature selection results
# 16S
multiStripchart(otu_plot, cleared_shared_otu, colonized_shared_otu, shared_otu_pval, shared_otu_aucrf_oob, 
                'Cleared', 'Colonized', 'white', 'darkorchid4', '', 'black', 
                formatted_names, expression(paste('Relative Abundance (',log[10],')')))
# Metabolome
multiStripchart(metabolome_plot, cleared_metabolome, colonized_metabolome, metabolome_pval, metabolome_aucrf_oob, 
                'Cleared', 'Colonized', 'white', 'darkorchid4', '', 'black', 
                as.list(colnames(cleared_metabolome)), expression(paste('Scaled Intensity (',log[10],')')))

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()
