# Set up environment

# Load dependencies
deps <- c('vegan', 'plotrix', 'reshape2', 'randomForest')
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Set seed for RNG
set.seed(6189)

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
metabolome <- metabolome[!colnames(metabolome) %in% c('CefC5M2'), ] # Contaminated sample
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
shared_otu <- shared_otu[ , !(names(shared_otu) == 'Otu0004')] # Remove residual C. difficile OTU
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('StrepC4M1','StrepC4M2','StrepC4M3'), ]

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
metabolome_germfree <- subset(metabolome_nmds, abx == 'germfree')
metabolome_untreated <- subset(metabolome_nmds, abx == 'none')

# Community structure
otu_nmds <- metaMDS(shared_otu, k=2, trymax=100)$points
otu_nmds[,1] <- otu_nmds[,1] - 0.2
otu_nmds[,2] <- otu_nmds[,2] - 0.2
otu_nmds <- clean_merge(metadata, otu_nmds)
otu_nmds <- subset(otu_nmds, abx != 'germfree')
otu_cefoperazone <- subset(otu_nmds, abx == 'cefoperazone')
otu_clindamycin <- subset(otu_nmds, abx == 'clindamycin')
otu_streptomycin <- subset(otu_nmds, abx == 'streptomycin')
otu_untreated <- subset(otu_nmds, abx == 'none')

#----------------#

# Test differences between groups

#adonis()




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

# Filter subset tables




# Run random forest
all_infection_rf <- featureselect_RF(all_infection, 'infection') # OOB = 7.69%
conv_infection_rf <- featureselect_RF(conv_infection, 'infection') # OOB = 10.71%
abx_infection_rf <- featureselect_RF(abx_infection, 'infection') # OOB = 11.86%
cef_infection_rf <- featureselect_RF(cef_infection, 'infection') # OOB = 0%    why???
strep_infection_rf <- featureselect_RF(strep_infection, 'infection') # OOB = 6.67%
clinda_infection_rf <- featureselect_RF(clinda_infection, 'infection') # OOB = 33.33%
rm(all_infection, conv_infection, cef_infection, strep_infection, clinda_infection, abx_infection)


#----------------#


# Calculate correlations for heatmap

cormat <- round(cor(mydata),2)
melted_cormat <- melt(cormat)








#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_file, width=8.5, height=11)
layout(matrix(c(1,2,
                3,4,
                3,4),
              nrow=3, ncol=2, byrow=TRUE))

#-------------------#

# OTUs alone

par(mar=c(3,4,1,1), las=1, mgp=c(2,0.75,0), xaxs='i', yaxs='i')
plot(x=otu_nmds$MDS1, y=otu_nmds$MDS2, xlim=c(-1.6,1.6), ylim=c(-1.5,1.5),
     xlab='NMDS axis 1', ylab='NMDS axis 2', xaxt='n', yaxt='n', pch=19, cex=0.2)
axis(side=1, at=seq(-1.6,1.6,0.4), labels=seq(-1.6,1.6,0.4))
axis(side=2, at=seq(-1.5,1.5,0.3), labels=seq(-1.5,1.5,0.3))
points(x=otu_cefoperazone$MDS1, y=otu_cefoperazone$MDS2, bg='chartreuse3', pch=21, cex=1.7, lwd=1.2)
points(x=otu_clindamycin$MDS1, y=otu_clindamycin$MDS2, bg='blue2', pch=21, cex=1.7, lwd=1.2)
points(x=otu_streptomycin$MDS1, y=otu_streptomycin$MDS2, bg='firebrick1', pch=21, cex=1.7, lwd=1.2)
points(x=otu_untreated$MDS1, y=otu_untreated$MDS2, bg='azure2', pch=21, cex=1.7, lwd=1.2)
legend('bottomright', legend=c('Untreated (SPF)','Streptomycin-treated (SPF)','Cefoperzone-treated (SPF)','Clindamycin-treated (SPF)'), 
       pt.bg=c('azure2','firebrick1','chartreuse3','blue2'), 
       pch=21, pt.cex=1.7, bty='n')


draw.ellipse(x=-0.27, y=-0.13, a=0.15, b=0.1, angle=-100, lty=2, lwd=2) # resistant
text(x=-0.27, y=0.04, labels='Resistant', font=2)


mtext('A', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.3)

#-------------------#

# Metabolomics alone

par(mar=c(3,4,1,1), las=1, mgp=c(2,0.75,0), xaxs='i', yaxs='i')
plot(x=metabolome_nmds$MDS1, y=metabolome_nmds$MDS2, xlim=c(-0.4,0.4), ylim=c(-0.3,0.3),
     xlab='NMDS axis 1', ylab='NMDS axis 2', xaxt='n', yaxt='n', pch=19, cex=0.2)
axis(side=1, at=seq(-0.4,0.4,0.1), labels=seq(-0.4,0.4,0.1))
axis(side=2, at=seq(-0.3,0.3,0.1), labels=c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3))
points(x=metabolome_cefoperazone$MDS1, y=metabolome_cefoperazone$MDS2, bg='chartreuse3', pch=21, cex=2, lwd=1.2)
points(x=metabolome_clindamycin$MDS1, y=metabolome_clindamycin$MDS2, bg='blue2', pch=21, cex=2, lwd=1.2)
points(x=metabolome_streptomycin$MDS1, y=metabolome_streptomycin$MDS2, bg='firebrick1', pch=21, cex=2, lwd=1.2)
points(x=metabolome_germfree$MDS1, y=metabolome_germfree$MDS2, bg='gold1', pch=21, cex=2, lwd=1.2)
points(x=metabolome_untreated$MDS1, y=metabolome_untreated$MDS2, bg='azure2', pch=21, cex=2, lwd=1.2)
draw.ellipse(x=0.19, y=-0.01, a=0.28, b=0.17, angle=-60, lty=2, lwd=2) # susceptible
text(x=0.25, y=0.25, labels='Susceptible', font=2)
draw.ellipse(x=-0.27, y=-0.13, a=0.15, b=0.1, angle=-100, lty=2, lwd=2) # resistant
text(x=-0.27, y=0.04, labels='Resistant', font=2)
legend('topleft', legend=c('Untreated (SPF)','Streptomycin-treated (SPF)','Cefoperzone-treated (SPF)','Clindamycin-treated (SPF)','Germfree'), 
       pt.bg=c('azure2','firebrick1','chartreuse3','blue2','gold1'), 
       pch=21, pt.cex=2, bty='n')

mtext('B', side=2, line=2, las=2, adj=1.5, padj=-9, cex=1.3)

#-------------------#

plot(0, type='n', axes=FALSE, xlab='', ylab='')
# Feature selection
mtext('C', side=2, line=2, las=2, adj=1.5, padj=-22, cex=1.3)

#-------------------#

plot(0, type='n', axes=FALSE, xlab='', ylab='')
# Heatmap or correlation somehow...
mtext('D', side=2, line=2, las=2, adj=1.5, padj=-22, cex=1.3)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up

#for (dep in deps){
#  pkg <- paste('package:', dep, sep='')
#  detach(pkg, character.only = TRUE)
#}
#rm(list=ls())
#gc()

