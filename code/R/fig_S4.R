
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Output plot name
plot_abc <- 'results/supplement/figures/figure_S4abc.pdf'
plot_d <- 'results/supplement/figures/figure_S4d.pdf'
plot_e <- 'results/supplement/figures/figure_S4e.pdf'
plot_f <- 'results/supplement/figures/figure_S4f.pdf'

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

#-------------------------------------------------------------------------------------------------------------------------#

# Stats
# Prep data
metabolome_subset <- clean_merge(metadata, metabolome)
inf_subset <- subset(metabolome_subset, abx != 'none')
inf_subset$susceptibility <- NULL
inf_subset$abx <- NULL
abx_subset <- subset(metabolome_subset, abx != 'none')
abx_subset$infection <- NULL
abx_subset$susceptibility <- NULL
metabolome_subset$abx <- NULL
metabolome_subset$infection <- NULL

# Calculate significant differences
noabx_p <- round(adonis(metabolome_subset[,2:ncol(metabolome_subset)]~factor(metabolome_subset$susceptibility), data=metabolome_subset, permutations=10000, method='jaccard')$aov.tab[[6]][1], 3)
inf_p <- round(adonis(abx_subset[,2:ncol(inf_subset)]~factor(inf_subset$infection), data=inf_subset, permutations=10000, method='jaccard')$aov.tab[[6]][1], 3)
abx_p <- round(adonis(abx_subset[,2:ncol(abx_subset)]~factor(abx_subset$abx), data=abx_subset, permutations=10000, method='jaccard')$aov.tab[[6]][1], 3)

#-------------------------------------------------------------------------------------------------------------------------#

# Calculate axes
# Metabolome - all
metabolome_nmds <- metaMDS(metabolome, k=2, trymax=100, distance='bray')$points
metabolome_nmds <- clean_merge(metadata, metabolome_nmds)
metabolome_nmds$MDS1 <- metabolome_nmds$MDS1 + 0.15
metabolome_nmds$MDS2 <- metabolome_nmds$MDS2 - 0.05

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


# Separate analysis for abx only
abx_metabolome <- metabolome[,!colnames(metabolome) %in% c('ConvC1M1','ConvC1M2','ConvC1M3','ConvC1M4',
                                                           'ConvC2M1','ConvC2M2','ConvC2M3','ConvC2M4','ConvC2M5')] # Untreated SPF samples 

# Calculate axes and merge with metadata
abx_metabolome_nmds <- metaMDS(abx_metabolome, k=2, trymax=100)$points
abx_metabolome_nmds <- clean_merge(metadata, abx_metabolome_nmds)
abx_metabolome_nmds$MDS1 <- abx_metabolome_nmds$MDS1 - 0.05

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
metabolome_abx_630 <- subset(abx_metabolome_nmds, infection == '630')
metabolome_abx_mock <- subset(abx_metabolome_nmds, infection == 'mock')

#-------------------------------------------------------------------------------------------------------------------------#

# Feature selection
# Separate groups
metabolome <- clean_merge(metadata, metabolome)
abx_metabolome <- subset(metabolome, abx != 'none')
abx_metabolome$infection <- NULL
abx_metabolome$susceptibility <- NULL
abx_metabolome$abx <- factor(abx_metabolome$abx)
infection_metabolome <- subset(metabolome, abx != 'none')
infection_metabolome$abx <- NULL
infection_metabolome$susceptibility <- NULL
metabolome$infection <- NULL
metabolome$abx <- NULL
metabolome$abx <- factor(metabolome$susceptibility)
rm(metadata)

# Random Forest
all_rf <- featureselect_RF(metabolome, 'susceptibility')
abx_rf <- featureselect_RF(abx_metabolome, 'abx')
inf_rf <- featureselect_RF(infection_metabolome, 'infection')

# Sort and subset top hits
all_rf <- all_rf[order(-all_rf$MDA),][1:10,]
abx_rf <- abx_rf[order(-abx_rf$MDA),][1:10,]
inf_rf <- inf_rf[order(-inf_rf$MDA),][1:10,]

# Subset concentrations
res_metabolome <- subset(metabolome, susceptibility == 'resistant')[, all_rf$feature]
res_metabolome$susceptibility <- NULL
sus_metabolome <- subset(metabolome, susceptibility == 'susceptible')[, all_rf$feature]
sus_metabolome$susceptibility <- NULL
rm(metabolome)

inf_infection_metabolome <- subset(infection_metabolome, infection == '630')[, inf_rf$feature]
inf_infection_metabolome$infection <- NULL
mock_infection_metabolome <- subset(infection_metabolome, infection == 'mock')[, inf_rf$feature]
mock_infection_metabolome$infection <- NULL
rm(infection_metabolome)

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

infection_pvalues <- c()
for (i in 1:ncol(inf_infection_metabolome)){infection_pvalues[i] <- wilcox.test(inf_infection_metabolome[,i], mock_infection_metabolome[,i], exact=FALSE)$p.value}
infection_pvalues <- p.adjust(infection_pvalues, method='BH')

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
pdf(file=plot_abc, width=12, height=4)
layout(matrix(c(1,2,3),
              nrow=1, ncol=3, byrow=TRUE))
par(mar=c(4,4,1,1), las=1, mgp=c(2.8,0.75,0))

# All conventional mice
plot(x=metabolome_nmds$MDS1, y=metabolome_nmds$MDS2, xlim=c(-0.25,0.25), ylim=c(-0.15,0.15),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('a', side=2, line=2, las=2, adj=1.8, padj=-9.5, cex=1.4, font=2)
points(x=metabolome_cefoperazone_630$MDS1, y=metabolome_cefoperazone_630$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_clindamycin_630$MDS1, y=metabolome_clindamycin_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_streptomycin_630$MDS1, y=metabolome_streptomycin_630$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_cefoperazone_mock$MDS1, y=metabolome_cefoperazone_mock$MDS2, bg=cef_col, pch=24, cex=1.8, lwd=1.2)
points(x=metabolome_clindamycin_mock$MDS1, y=metabolome_clindamycin_mock$MDS2, bg=clinda_col, pch=24, cex=1.8, lwd=1.2)
points(x=metabolome_streptomycin_mock$MDS1, y=metabolome_streptomycin_mock$MDS2, bg=strep_col, pch=24, cex=1.8, lwd=1.2)
points(x=metabolome_noantibiotics$MDS1, y=metabolome_noantibiotics$MDS2, bg=noabx_col, pch=24, cex=1.8, lwd=1.2)
legend('bottomleft', legend=c('Resistant vs Susceptible', as.expression(bquote(paste(italic('p'),' < 0.001 ***')))), 
       pch=1, cex=1.4, pt.cex=0, bty='n')
legend('topright', legend=c('No Antibiotic','Streptomycin','Cefoperzone','Clindamycin'), 
       pt.bg=c(noabx_col,strep_col,cef_col,clinda_col), pch=22, cex=1.4, pt.cex=2.5)
legend('topleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(15,17), cex=1.2, pt.cex=2)
# Infection
plot(x=abx_metabolome_nmds$MDS1, y=abx_metabolome_nmds$MDS2, xlim=c(-0.075,0.075), ylim=c(-0.075,0.075),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('b', side=2, line=2, las=2, adj=1.8, padj=-9.5, cex=1.4, font=2)
points(x=metabolome_abx_630$MDS1, y=metabolome_abx_630$MDS2, bg='mediumorchid4', pch=21, cex=2, lwd=1.2)
points(x=metabolome_abx_mock$MDS1, y=metabolome_abx_mock$MDS2, bg='chartreuse2', pch=21, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Between Antibiotic Groups', as.expression(bquote(paste(italic('p'),' = 0.075 n.s.')))), 
       pch=1, cex=1.4, pt.cex=0, bty='n')
legend('topleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pt.bg=c('mediumorchid4','chartreuse2'), pch=21, cex=1.2, pt.cex=2)
# Antibiotics individually
plot(x=abx_metabolome_nmds$MDS1, y=abx_metabolome_nmds$MDS2, xlim=c(-0.075,0.075), ylim=c(-0.075,0.075),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('c', side=2, line=2, las=2, adj=1.8, padj=-9.5, cex=1.4, font=2)
points(x=metabolome_abx_cef_630$MDS1, y=metabolome_abx_cef_630$MDS2, bg=cef_col, pch=22, cex=2, lwd=1.2)
points(x=metabolome_abx_cef_mock$MDS1, y=metabolome_abx_cef_mock$MDS2, bg=cef_col, pch=24, cex=2, lwd=1.2)
points(x=metabolome_abx_clinda_630$MDS1, y=metabolome_abx_clinda_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_abx_clinda_mock$MDS1, y=metabolome_abx_clinda_mock$MDS2, bg=clinda_col, pch=24, cex=2, lwd=1.2)
points(x=metabolome_abx_strep_630$MDS1, y=metabolome_abx_strep_630$MDS2, bg=strep_col, pch=22, cex=2, lwd=1.2)
points(x=metabolome_abx_strep_mock$MDS1, y=metabolome_abx_strep_mock$MDS2, bg=strep_col, pch=24, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Between Antibiotic Groups', as.expression(bquote(paste(italic('p'),' < 0.001 ***')))), 
       pch=1, cex=1.4, pt.cex=0, bty='n')
legend('topleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(15,17), cex=1.2, pt.cex=2)

dev.off()

#---------------#

# Feature Selection
# All abx vs Untreated
metabolite_stripchart(plot_d, res_metabolome, sus_metabolome, resistant_pvalues, all_rf$MDA, 0, 'Resistant', 'Susceptible', '', 'white', 'd')
# All Infected vs All Mock
metabolite_stripchart(plot_e, inf_infection_metabolome, mock_infection_metabolome, infection_pvalues, inf_rf$MDA, 11.11, 'Infected', 'Mock', '', 'white', 'e')

# Each antibiotic group
pdf(file=plot_f, width=4, height=ncol(strep_abx_metabolome)*1.5)
layout(matrix(c(1:(ncol(strep_abx_metabolome)+2)), nrow=(ncol(strep_abx_metabolome)+2), ncol=1, byrow = TRUE))

par(mar=c(0.2, 0, 0, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE)
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
text(x=-10.2, y=-3, labels='f', cex=2.4, font=2, xpd=TRUE)
legend('bottomright', legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin'), bty='n',
       pt.bg=c(strep_col,cef_col,clinda_col), pch=21, cex=1.2, pt.cex=2, ncol=3)

par(mar=c(0.2, 2, 0.2, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE, yaxs='i')
for(i in c(1:(ncol(strep_abx_metabolome)))){
  xmax <- ceiling(max(c(max(strep_abx_metabolome[,i]), max(cef_abx_metabolome[,i]))))
  while(xmax %% 5 != 0 ){xmax <- xmax + 1}
  if (xmax > 70){
    while(xmax %% 10 != 0 ){xmax <- xmax + 1}
  }
  plot(0, type='n', xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,xmax), ylim=c(0.3,2.25))
  stripchart(at=1.65, jitter(strep_abx_metabolome[,i], amount=1e-5),
             pch=21, bg=strep_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  stripchart(at=1.2, jitter(cef_abx_metabolome[,i], amount=1e-5), 
             pch=21, bg=cef_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  stripchart(at=0.66, jitter(clinda_abx_metabolome[,i], amount=1e-5), 
             pch=21, bg=clinda_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  metabolite <- paste(colnames(strep_abx_metabolome)[i], ' [',as.character(round(abx_rf$MDA[i],3)),']', sep='')
  legend('topright', legend=metabolite, pch=1, cex=1.3, pt.cex=0, bty='n')
  if (xmax <= 10) {
    text(x=seq(0,xmax,1), y=0.42, labels=seq(0,xmax,1), cex=1)
    axis(1, at=seq(0,5,1), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 70){
    text(x=seq(0,xmax,10), y=0.42, labels=seq(0,xmax,10), cex=1)
    axis(1, at=seq(0,xmax,10), NA, cex.axis=0.8, tck=0.015)
  } else {
    text(x=seq(0,xmax,5), y=0.42, labels=seq(0,xmax,5), cex=1)
    axis(1, at=seq(0,xmax,5), NA, cex.axis=0.8, tck=0.015)
  }
  segments(median(strep_abx_metabolome[,i]), 1.58, median(strep_abx_metabolome[,i]), 1.72, lwd=2.5)
  segments(median(cef_abx_metabolome[,i]), 1.03, median(cef_abx_metabolome[,i]), 1.37, lwd=2.5)
  segments(median(clinda_abx_metabolome[,i]), 0.49, median(clinda_abx_metabolome[,i]), 0.83, lwd=2.5)
  
}
par(mar=c(0, 0, 0, 0))
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
text(x=0, y=4, labels=expression(paste('Scaled Intensity (',log[10],')')), cex=1.4)
text(x=8, y=4.5, labels='OOB Error = 1.85%', cex=0.8)

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
