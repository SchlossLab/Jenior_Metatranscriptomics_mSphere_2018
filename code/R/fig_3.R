
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Output plot name
plot_ab <- 'results/figures/figure_3ab.pdf'
plot_c <- 'results/figures/figure_3c.pdf'
plot_d <- 'results/figures/figure_3d.pdf'

# Input Metabolomes
metabolome_file <- 'data/metabolome/scaled_intensities.log10.tsv'

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

#-------------------------------------------------------------------------------------------------------------------------#

# Stats
# Prep data
metabolome <- clean_merge(metadata, metabolome)
metabolome <- subset(metabolome, abx != 'germfree')
metabolome <- subset(metabolome, infection == 'mock')
abx_metabolome <- subset(metabolome, abx != 'none')
abx_metabolome$infection <- NULL
abx_metabolome$susceptibility <- NULL
metabolome$abx <- NULL
metabolome$infection <- NULL

# Calculate significant differences
noabx_p <- round(adonis(metabolome[,2:ncol(metabolome)]~factor(metabolome$susceptibility), data=metabolome, permutations=10000, method='jaccard')$aov.tab[[6]][1], 3)
abx_p <- round(adonis(abx_metabolome[,2:ncol(abx_metabolome)]~factor(abx_metabolome$abx), data=abx_metabolome, permutations=10000, method='jaccard')$aov.tab[[6]][1], 3)
metabolome$susceptibility <- NULL
abx_metabolome$abx <-  NULL

#-------------------------------------------------------------------------------------------------------------------------#

# Calculate axes
# Metabolome - all
metabolome_nmds <- metaMDS(metabolome, k=2, trymax=100, distance='bray')$points
metabolome_nmds <- clean_merge(metadata, metabolome_nmds)
metabolome_nmds$MDS1 <- metabolome_nmds$MDS1 - 0.1

# Subset NMDS axes to color points - all
metabolome_cefoperazone <- subset(metabolome_nmds, abx == 'cefoperazone')
metabolome_cefoperazone_mock <- subset(metabolome_cefoperazone, infection == 'mock')
rm(metabolome_cefoperazone)
metabolome_clindamycin <- subset(metabolome_nmds, abx == 'clindamycin')
metabolome_clindamycin_mock <- subset(metabolome_clindamycin, infection == 'mock')
rm(metabolome_clindamycin)
metabolome_streptomycin <- subset(metabolome_nmds, abx == 'streptomycin')
metabolome_streptomycin_mock <- subset(metabolome_streptomycin, infection == 'mock')
rm(metabolome_streptomycin)
metabolome_noantibiotics <- subset(metabolome_nmds, abx == 'none')

# Separate analysis for abx only
# Calculate axes and merge with metadata
abx_metabolome_nmds <- metaMDS(abx_metabolome, k=2, trymax=100)$points
abx_metabolome_nmds <- clean_merge(metadata, abx_metabolome_nmds)
abx_metabolome_nmds$MDS1 <- abx_metabolome_nmds$MDS1 - 0.06

# Subset to points for plot
metabolome_abx <- subset(abx_metabolome_nmds, abx == 'cefoperazone')
metabolome_abx_cef_mock <- subset(metabolome_abx, infection == 'mock')
metabolome_abx <- subset(abx_metabolome_nmds, abx == 'clindamycin')
metabolome_abx_clinda_mock <- subset(metabolome_abx, infection == 'mock')
metabolome_abx <- subset(abx_metabolome_nmds, abx == 'streptomycin')
metabolome_abx_strep_mock <- subset(metabolome_abx, infection == 'mock')

#-------------------------------------------------------------------------------------------------------------------------#

# Feature selection
# Separate groups
metabolome <- clean_merge(metadata, metabolome)
abx_metabolome <- subset(metabolome, abx != 'none')
abx_metabolome$infection <- NULL
abx_metabolome$susceptibility <- NULL
abx_metabolome$abx <- factor(abx_metabolome$abx)
metabolome$infection <- NULL
metabolome$abx <- NULL
metabolome$abx <- factor(metabolome$susceptibility)
rm(metadata)

# Random Forest
all_rf <- featureselect_RF(metabolome, 'susceptibility')
abx_rf <- featureselect_RF(abx_metabolome, 'abx')

# Sort and subset top hits
all_rf <- all_rf[order(-all_rf$MDA),][1:7,]
abx_rf <- abx_rf[order(-abx_rf$MDA),][1:7,]

# Subset concentrations
res_metabolome <- subset(metabolome, susceptibility == 'resistant')[, all_rf$feature]
res_metabolome$susceptibility <- NULL
sus_metabolome <- subset(metabolome, susceptibility == 'susceptible')[, all_rf$feature]
sus_metabolome$susceptibility <- NULL
rm(metabolome)

cef_abx_metabolome <- subset(abx_metabolome, abx == 'cefoperazone')[, abx_rf$feature]
cef_abx_metabolome$abx <- NULL
clinda_abx_metabolome <- subset(abx_metabolome, abx == 'clindamycin')[, abx_rf$feature]
clinda_abx_metabolome$abx <- NULL
strep_abx_metabolome <- subset(abx_metabolome, abx == 'streptomycin')[, abx_rf$feature]
strep_abx_metabolome$abx <- NULL
rm(abx_metabolome)

# Find significant differences
resistant_pvalues <- c()
for (i in 1:ncol(res_metabolome)){resistant_pvalues[i] <- wilcox.test(res_metabolome[,i], sus_metabolome[,i], exact=FALSE)$p.value}
resistant_pvalues <- p.adjust(resistant_pvalues, method='BH')

abx_pvalues1 <- c()
for (i in 1:ncol(cef_abx_metabolome)){abx_pvalues1[i] <- wilcox.test(cef_abx_metabolome[,i], clinda_abx_metabolome[,i], exact=FALSE)$p.value}
abx_pvalues2 <- c()
for (i in 1:ncol(cef_abx_metabolome)){abx_pvalues2[i] <- wilcox.test(cef_abx_metabolome[,i], strep_abx_metabolome[,i], exact=FALSE)$p.value}
abx_pvalues3 <- c()
for (i in 1:ncol(clinda_abx_metabolome)){abx_pvalues3[i] <- wilcox.test(clinda_abx_metabolome[,i], strep_abx_metabolome[,i], exact=FALSE)$p.value}
abx_pvalues1 <- p.adjust(abx_pvalues1, method='BH')
abx_pvalues2 <- p.adjust(abx_pvalues2, method='BH')
abx_pvalues3 <- p.adjust(abx_pvalues3, method='BH')

#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_ab, width=12, height=6)
layout(matrix(c(1,2),
              nrow=1, ncol=2, byrow=TRUE))
par(mar=c(4,4,1,1), las=1, mgp=c(2.8,0.75,0))

# All conventional mice
plot(x=metabolome_nmds$MDS1, y=metabolome_nmds$MDS2, xlim=c(-0.25,0.25), ylim=c(-0.15,0.15),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('a', side=2, line=2, las=2, adj=1.9, padj=-10, cex=2, font=2)
points(x=metabolome_cefoperazone_mock$MDS1, y=metabolome_cefoperazone_mock$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_clindamycin_mock$MDS1, y=metabolome_clindamycin_mock$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_streptomycin_mock$MDS1, y=metabolome_streptomycin_mock$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_noantibiotics$MDS1, y=metabolome_noantibiotics$MDS2, bg=noabx_col, pch=21, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Resistant vs Susceptible', as.expression(bquote(paste(italic('p'),' < 0.001 ***')))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
legend('topright', legend=c('No Antibiotic','Streptomycin','Cefoperzone','Clindamycin'), 
       pt.bg=c(noabx_col,strep_col,cef_col,clinda_col), pch=21, cex=1.2, pt.cex=2.5)

# Antibiotics individually
plot(x=abx_metabolome_nmds$MDS1, y=abx_metabolome_nmds$MDS2, xlim=c(-0.15,0.15), ylim=c(-0.15,0.15),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('b', side=2, line=2, las=2, adj=1.9, padj=-10, cex=2, font=2)
points(x=metabolome_abx_cef_mock$MDS1, y=metabolome_abx_cef_mock$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_abx_clinda_mock$MDS1, y=metabolome_abx_clinda_mock$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_abx_strep_mock$MDS1, y=metabolome_abx_strep_mock$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Between Pretreatments', as.expression(bquote(paste(italic('p'),' < 0.001 ***')))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
legend('topright', legend=c('Streptomycin','Cefoperzone','Clindamycin'), 
       pt.bg=c(strep_col,cef_col,clinda_col), pch=21, cex=1.2, pt.cex=2.5)

dev.off()

#---------------#

# Feature Selection
# All abx vs Untreated
metabolite_stripchart(plot_c, res_metabolome, sus_metabolome, resistant_pvalues, all_rf$MDA, 
                      0, 'Resistant', 'Susceptible', '', 'white', 'c', 'gray80', 'salmon2')

# Each antibiotic group
pdf(file=plot_d, width=4, height=ncol(strep_abx_metabolome)*1.5)
layout(matrix(c(1:(ncol(strep_abx_metabolome)+2)), nrow=(ncol(strep_abx_metabolome)+2), ncol=1, byrow = TRUE))

par(mar=c(0.2, 0, 0, 2), mgp=c(2.3, 0.75, 0), xpd=FALSE)
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
text(x=-10.2, y=-3, labels='d', cex=2.4, font=2, xpd=TRUE)
legend('bottomright', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin'), bty='n',
       pt.bg=c(strep_col,cef_col,clinda_col), pch=21, cex=1.2, pt.cex=2, ncol=3)

par(mar=c(0.2, 2, 0.2, 2.5), mgp=c(2.3, 0.75, 0), xpd=FALSE, yaxs='i')
for(i in c(1:(ncol(strep_abx_metabolome)))){
  xmax <- ceiling(max(c(max(strep_abx_metabolome[,i]), max(cef_abx_metabolome[,i]), max(clinda_abx_metabolome[,i]))))
  while(xmax %% 5 != 0 ){xmax <- xmax + 1}
  plot(0, type='n', xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,xmax), ylim=c(0.3,2.25))
  stripchart(at=1.65, jitter(strep_abx_metabolome[,i], amount=1e-5),
             pch=21, bg=strep_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  stripchart(at=1.2, jitter(cef_abx_metabolome[,i], amount=1e-5), 
             pch=21, bg=cef_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  stripchart(at=0.66, jitter(clinda_abx_metabolome[,i], amount=1e-5), 
             pch=21, bg=clinda_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  metabolite <- paste(colnames(strep_abx_metabolome)[i], ' [',as.character(round(abx_rf$MDA[i],3)),']', sep='')
  legend('topright', legend=metabolite, pch=1, cex=1.3, pt.cex=0, bty='n')
  
  mtext('|', side=4, cex=1.6, las=1, adj=-0.3, padj=-0.5)
  mtext('|', side=4, cex=1.6, las=1, adj=-0.3, padj=-0.3)
  mtext('|', side=4, cex=1.6, las=1, adj=-0.3, padj=1.7)
  mtext('|', side=4, cex=1.6, las=1, adj=-0.3, padj=1.5)
  mtext('|', side=4, cex=1.6, las=1, adj=-2.5, padj=-0.5)
  mtext('|', side=4, cex=1.6, las=1, adj=-2.5, padj=-0.3)
  mtext('|', side=4, cex=1.6, las=1, adj=-2.5, padj=0)
  mtext('|', side=4, cex=1.6, las=1, adj=-2.5, padj=0.5)
  mtext('|', side=4, cex=1.6, las=1, adj=-2.5, padj=1.5)
  mtext('|', side=4, cex=1.6, las=1, adj=-2.5, padj=1.7)
  if (abx_pvalues2[i] < 0.001){
    mtext('***', side=4, font=2, cex=1, padj=0.8, adj=0.64)
  } else if (abx_pvalues2[i] <= 0.01){
    mtext('**', side=4, font=2, cex=1, padj=0.8, adj=0.64)
  } else if (abx_pvalues2[i] <= 0.05){
    mtext('*', side=4, font=2, cex=1, padj=0.8, adj=0.64)
  } else {
    mtext('n.s.', cex=0.7, side=4, padj=0.5, adj=0.64)
  }
  if (abx_pvalues1[i] < 0.001){
    mtext('***', side=4, font=2, cex=1, padj=0.8, adj=0.27)
  } else if (abx_pvalues1[i] <= 0.01){
    mtext('**', side=4, font=2, cex=1, padj=0.8, adj=0.27)
  } else if (abx_pvalues1[i] <= 0.05){
    mtext('*', side=4, font=2, cex=1, padj=0.8, adj=0.27)
  } else {
    mtext('n.s.', cex=0.7, side=4, padj=0.5, adj=0.27)
  }
  if (abx_pvalues3[i] < 0.001){
    mtext('***', side=4, font=2, cex=1, padj=2.2, adj=0.45)
  } else if (abx_pvalues3[i] <= 0.01){
    mtext('**', side=4, font=2, cex=1, padj=2.2, adj=0.45)
  } else if (abx_pvalues3[i] <= 0.05){
    mtext('*', side=4, font=2, cex=1, padj=2.2, adj=0.45)
  } else {
    mtext('n.s.', cex=0.7, side=4, padj=2.4, adj=0.45)
  }
  
  if (xmax <= 10) {
    text(x=seq(0,xmax,1), y=0.42, labels=seq(0,xmax,1), cex=1)
    axis(1, at=seq(0,5,1), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 1000){
    text(x=seq(0,xmax,200), y=0.42, labels=seq(0,xmax,200), cex=1)
    axis(1, at=seq(0,xmax,200), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 500){
    text(x=seq(0,xmax,100), y=0.42, labels=seq(0,xmax,100), cex=1)
    axis(1, at=seq(0,xmax,100), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 100){
    text(x=seq(0,xmax,50), y=0.42, labels=seq(0,xmax,50), cex=1)
    axis(1, at=seq(0,xmax,50), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 50){
    text(x=seq(0,xmax,10), y=0.42, labels=seq(0,xmax,10), cex=1)
    axis(1, at=seq(0,xmax,10), NA, cex.axis=0.8, tck=0.015)
  } else {
    text(x=seq(0,xmax,5), y=0.42, labels=seq(0,xmax,5), cex=1)
    axis(1, at=seq(0,xmax,5), NA, cex.axis=0.8, tck=0.015)
  }
  segments(median(strep_abx_metabolome[,i]), 1.48, median(strep_abx_metabolome[,i]), 1.82, lwd=2.5)
  segments(median(cef_abx_metabolome[,i]), 1.03, median(cef_abx_metabolome[,i]), 1.37, lwd=2.5)
  segments(median(clinda_abx_metabolome[,i]), 0.49, median(clinda_abx_metabolome[,i]), 0.83, lwd=2.5)
  
}
par(mar=c(0, 0, 0, 0))
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
text(x=0, y=4, labels=expression(paste('Scaled Intensity (',log[10],')')), cex=1.4)
text(x=7, y=4.5, labels='OOB Error = 1.85%')

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
