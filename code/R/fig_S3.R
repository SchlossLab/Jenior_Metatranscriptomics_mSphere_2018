
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Output plot name
plot_a <- 'results/supplement/figures/figure_S3a.pdf'
plot_b <- 'results/supplement/figures/figure_S3b.pdf'

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

# Calculate significant differences
inf_p <- round(adonis(abx_subset[,2:ncol(inf_subset)]~factor(inf_subset$infection), data=inf_subset, permutations=10000, method='jaccard')$aov.tab[[6]][1], 3)

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
metabolome_abx_630 <- subset(abx_metabolome_nmds, infection == '630')
metabolome_abx_mock <- subset(abx_metabolome_nmds, infection == 'mock')

#-------------------------------------------------------------------------------------------------------------------------#

# Feature selection
# Separate groups
metabolome <- clean_merge(metadata, metabolome)
infection_metabolome <- subset(metabolome, abx != 'none')
infection_metabolome$abx <- NULL
infection_metabolome$susceptibility <- NULL
rm(metadata)

# Random Forest
inf_rf <- featureselect_RF(infection_metabolome, 'infection')

# Sort and subset top hits
inf_rf <- inf_rf[order(-inf_rf$MDA),][1:10,]

# Subset concentrations
inf_infection_metabolome <- subset(infection_metabolome, infection == '630')[, inf_rf$feature]
inf_infection_metabolome$infection <- NULL
mock_infection_metabolome <- subset(infection_metabolome, infection == 'mock')[, inf_rf$feature]
mock_infection_metabolome$infection <- NULL
rm(infection_metabolome)

# Find significant differences
infection_pvalues <- c()
for (i in 1:ncol(inf_infection_metabolome)){infection_pvalues[i] <- wilcox.test(inf_infection_metabolome[,i], mock_infection_metabolome[,i], exact=FALSE)$p.value}
infection_pvalues <- p.adjust(infection_pvalues, method='BH')

#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_a, width=6, height=6)
par(mar=c(4,4,1,1), las=1, mgp=c(2.8,0.75,0))

# Infection
plot(x=abx_metabolome_nmds$MDS1, y=abx_metabolome_nmds$MDS2, xlim=c(-0.075,0.075), ylim=c(-0.075,0.075),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('A', side=2, line=2, las=2, adj=1.8, padj=-10, cex=2, font=2)
points(x=metabolome_abx_630$MDS1, y=metabolome_abx_630$MDS2, bg='mediumorchid4', pch=21, cex=2, lwd=1.2)
points(x=metabolome_abx_mock$MDS1, y=metabolome_abx_mock$MDS2, bg='chartreuse2', pch=21, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Between Antibiotic Groups', as.expression(bquote(paste(italic('p'),' = 0.075 n.s.')))), 
       pch=1, cex=1.4, pt.cex=0, bty='n')
legend('topleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pt.bg=c('mediumorchid4','chartreuse2'), pch=21, cex=1.2, pt.cex=2)


dev.off()

#---------------#

# Feature Selection
# All Infected vs All Mock
metabolite_stripchart(plot_b, inf_infection_metabolome, mock_infection_metabolome, infection_pvalues, 
                      inf_rf$MDA, 11.11, 'Infected', 'Mock', '', 'white', 'B', 'mediumorchid4', 'chartreuse2')

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()
