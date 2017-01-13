
# Load dependencies
deps <- c('wesanderson', 'randomForest', 'vegan');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}
set.seed(6189)

featureselect_RF <- function(training_data, feature){
  
  # Random Forest procedure based on Segal et al. (2004)
  # Set parameters
  attach(training_data)
  levels <- as.vector(unique(training_data[,feature]))
  factor_1 <- round(length(rownames(data[which(data[,feature]==levels[1]),])) * 0.623)
  factor_2 <- round(length(rownames(data[which(data[,feature]==levels[2]),])) * 0.623)
  factor <- max(c(round(factor_1 / factor_2), round(factor_2 / factor_1))) * 3
  n_tree <- round(length(colnames(training_data)) - 1) * factor
  m_try <- round(sqrt(length(colnames(training_data)) - 1))
  
  # Run random forest
  data_randomForest <- randomForest(training_data[,feature]~., data=training_data, importance=TRUE, replace=FALSE, do.trace=500, err.rate=TRUE, ntree=n_tree, mtry=m_try)
  detach(training_data)
  
  # Parse features for significance
  features_RF <- importance(data_randomForest, type=1)
  final_features_RF <- subset(features_RF, features_RF > abs(min(features_RF)))
  
  return(final_features_RF)
}

clean_merge <- function(data_1, data_2){
  
  clean_merged <- merge(data_1, data_2, by = 'row.names')
  rownames(clean_merged) <- clean_merged$Row.names
  clean_merged$Row.names <- NULL
  
  return(clean_merged)
}

#-------------------------------------------------------------------------------------------------------------------------------------#

# Define variables
shared_otu_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'
taxonomy_otu_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.cons.taxonomy'
shared_family_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.family.subsample.shared'
taxonomy_family_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.family.cons.taxonomy'
metadata_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metadata.tsv'

# Load in data
metadata <- read.delim(metadata_file, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2'), ] # Remove contaminated sample
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
taxonomy_otu <- read.delim(taxonomy_otu_file, sep='\t', header=T, row.names=1)
taxonomy_otu$Size <- NULL
shared_otu <- read.delim(shared_otu_file, sep='\t', header=T, row.names=2)
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('CefC5M2'), ]  # Remove contaminated sample
shared_otu$numOtus <- NULL
shared_otu$label <- NULL
taxonomy_family <- read.delim(taxonomy_family_file, sep='\t', header=T, row.names=1)
taxonomy_family$Size <- NULL
shared_family <- read.delim(shared_family_file, sep='\t', header=T, row.names=2)
shared_family <- shared_family[!rownames(shared_family) %in% c('CefC5M2'), ]  # Remove contaminated sample
shared_family$numOtus <- NULL
shared_family$label <- NULL

rm(shared_otu_file, taxonomy_otu_file, shared_family_file, taxonomy_family_file, metadata_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format data (previously subsampled in mothur)

# OTU shared file
metadata_shared_otu <- clean_merge(metadata, shared_otu)
conv_shared_otu <- subset(metadata_shared_otu, type == 'conventional')
conv_shared_otu$type <- NULL
conv_shared_otu$infection <- NULL

conv_tax_otu <- t(clean_merge(t(conv_shared_otu), taxonomy_otu))

cef_shared_otu <- subset(metadata_shared_otu, abx == 'cefoperazone')
cef_shared_otu$abx <- NULL
cef_shared_otu$type <- NULL
strep_shared_otu <- subset(metadata_shared_otu, abx == 'streptomycin')
strep_shared_otu$abx <- NULL
strep_shared_otu$type <- NULL
clinda_shared_otu <- subset(metadata_shared_otu, abx == 'clindamycin')
clinda_shared_otu$abx <- NULL
clinda_shared_otu$type <- NULL
rm(shared_otu, metadata_shared_otu)

# Phylotype family-level shared file
metadata_shared_family <- clean_merge(metadata, shared_family)
conv_shared_family <- subset(metadata_shared_family, type == 'conventional')
conv_shared_family$type <- NULL
conv_shared_family$infection <- NULL
conv_shared_tax_family <- t(clean_merge(t(conventional_shared_family), taxonomy_family))
rm(shared_family, taxonomy_family, metadata_shared_family, conventional_shared_family)

#----------------------------------------------------------------------------------------------------------------------#

# Analysis and stats

# NMDS axes
nmds_shared <- conv_shared_otu
nmds_shared$abx <- NULL


nmds <- metaMDS(nmds_shared, k=2, trymax=1000)
rm(nmds_shared)

# Subset nmds groups
nmds_points <- nmds$points
metadata_nmds_points <- clean_merge(metadata, nmds_points)
metadata_nmds_points$type <- NULL
cef_points <- subset(metadata_nmds_points, abx == 'cefoperazone')
cef_points$abx <- NULL
clinda_points <- subset(metadata_nmds_points, abx == 'clindamycin')
clinda_points$abx <- NULL
strep_points <- subset(metadata_nmds_points, abx == 'streptomycin')
strep_points$abx <- NULL
metadata_nmds_points$infection <- NULL

# AMOVA
adonis(metadata_nmds_points[,2:3]~abx, data=metadata_nmds_points, permutations=1000, method='bray')

adonis([,2:3]~infection, data=cef_points, permutations=1000, method='bray')
adonis([,2:3]~infection, data=clinda_points, permutations=1000, method='bray')
adonis([,2:3]~infection, data=strep_points, permutations=1000, method='bray')


# Run random forest
cef_features <- featureselect_RF(cef_shared, 'infection')
strep_features <- featureselect_RF(strep_shared, 'infection')
clinda_features <- featureselect_RF(clinda_shared, 'infection')
cef_features_tax <- clean_merge(cef_features, taxonomy_otu)
strep_features_tax <- clean_merge(strep_features, taxonomy_otu)
clinda_features_tax <- clean_merge(clinda_features, taxonomy_otu)

#----------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_2.pdf'
pdf(file=plot_file, width=12, height=7)
layout(matrix(c(1,2,2,
                4,5,6),
              nrow=2, ncol=3, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# NMDS of treatment groups
par(las=1, mar=c(5,5,1,1))
plot(metadata_axes$axis1, metadata_axes$axis2, pch=21, cex=0,
     xlim=c(-0.8,0.8), ylim=c(-0.8,0.8), cex.lab=2, cex.axis=1.7,
     xlab='NMDS Axis 1', ylab='NMDS Axis 2')

# add mock points
points(x=mock_axes$axis1, y=mock_axes$axis2, 
       col=c(wes_palette("FantasticFox")[3], 
             wes_palette("FantasticFox")[5], 
             'forestgreen', 'black', 
             wes_palette("FantasticFox")[1])[mock_axes$abx], 
       pch=1, lwd=4, cex=2.5)
# add infected points
points(x=infected_axes$axis1, y=infected_axes$axis2, 
       col=c(wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen', 'black', wes_palette("FantasticFox")[1])[infected_axes$abx], 
       pch=2, lwd=4, cex=2.5)

# Add legends
legend('topleft', legend=c('Streptomycin-treated', 'Cefoperzone-treated', 'Clindamycin-treated', 'No Antibiotics'), 
       col=c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'black'), 
       pch=15, cex=1.9, pt.cex=3, bty='n')
legend('bottomleft', legend=c('Mock Infected', '630 Infected'), 
       col='black', pch=c(16,17), cex=2, pt.cex=2.5, bty='n')

mtext('A', side=2, line=2, las=2, adj=1.7, padj=-18.1, cex=1.3)

#-----------------------#

# Family-level phylotype bar chart

par(las=1, mar=c(4,4.5,1,12), xpd=TRUE)

# When needed, use this pallete   
final_colors <- c("gold1", "orangered1", "aquamarine3", "firebrick", "forestgreen", "blue3", 
                  "mediumorchid2", "violetred4", "mediumpurple4", "dodgerblue3", "goldenrod3", "chartreuse3")


# Plot the final formatted table
barplot(t(filtered_shared), col=final_colors, yaxt='n', ylim=c(0,100), ylab='% Relative Abundance', font=2)
box()
axis(side=2, at=seq(0,100,20), tick=TRUE)
segments(x0=rep(0,4), y0=seq(20,80,20), x1=rep(5,4), y1=seq(20,80,20), lty=2)

# Create a figure legend in the margin
legend(5.025, 75, legend=taxa, pt.bg=final_colors, pch=22, pt.cex=1.3, cex=0.7)

mtext('B', side=2, line=2, las=2, adj=1.7, padj=-10.5, cex=1.3)

#-----------------------#

# Random Forest results

# Cefoperazone plot
dotchart()


# Clindamycin plot
dotchart()


# Streptomycin plot
dotchart()



#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
dev.off()


rm(plot_file, metadata_axes, mock_axes, infected_axes, metadata_summary, 
   strep_div, cef_div, clinda_div, conv_div)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()