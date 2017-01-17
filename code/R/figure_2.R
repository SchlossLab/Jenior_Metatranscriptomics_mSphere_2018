
# Load dependencies
deps <- c('wesanderson', 'randomForest', 'vegan')
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
  factor_1 <- round(length(rownames(training_data[which(training_data[,feature]==levels[1]),])) * 0.623)
  factor_2 <- round(length(rownames(training_data[which(training_data[,feature]==levels[2]),])) * 0.623)
  factor <- max(c(round(factor_1 / factor_2), round(factor_2 / factor_1))) * 3
  n_tree <- round(length(colnames(training_data)) - 1) * factor
  m_try <- round(sqrt(length(colnames(training_data)) - 1))
  
  # Run random forest
  data_randomForest <- randomForest(training_data[,feature]~., data=training_data, importance=TRUE, replace=FALSE, do.trace=500, err.rate=TRUE, ntree=n_tree, mtry=m_try)
  detach(training_data)
  
  # Parse features for significance
  features_RF <- importance(data_randomForest, type=1)
  final_features_RF <- subset(features_RF, features_RF > abs(min(features_RF)))
  final_features_RF <- final_features_RF[!(rownames(final_features_RF) == feature),]
  final_features_RF <- as.data.frame(final_features_RF)

  return(final_features_RF)
}

clean_merge <- function(data_1, data_2){
  
  clean_merged <- merge(data_1, data_2, by = 'row.names')
  rownames(clean_merged) <- clean_merged$Row.names
  clean_merged$Row.names <- NULL
  
  return(clean_merged)
}

#-------------------------------------------------------------------------------------------------------------------------------------#

# Define variabless

nmds_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.unique_list.conventional.thetayc.0.03.lt.ave.nmds.axes'
shared_otu_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'
taxonomy_otu_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.cons.genus.format.taxonomy'
shared_family_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.family.subsample.shared'
taxonomy_family_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.family.cons.family.format.taxonomy'
metadata_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metadata.tsv'

# Load in data
nmds <- read.delim(nmds_file, sep='\t', header=T, row.names=1)
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

rm(nmds_file, shared_otu_file, taxonomy_otu_file, shared_family_file, taxonomy_family_file, metadata_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format data (previously subsampled in mothur)

# OTU shared file
metadata_shared_otu <- clean_merge(metadata, shared_otu)
conv_shared_otu <- subset(metadata_shared_otu, type == 'conventional')
conv_shared_otu$type <- NULL
conv_shared_otu$infection <- NULL
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

# Run random forest
cef_features <- featureselect_RF(cef_shared_otu, 'infection')
strep_features <- featureselect_RF(strep_shared_otu, 'infection')
clinda_features <- featureselect_RF(clinda_shared_otu, 'infection')
cef_features_tax <- clean_merge(cef_features, taxonomy_otu)
cef_features_tax <- cef_features_tax[order(cef_features_tax[,1]),] 
cef_features_tax$Taxonomy <- gsub('_', ' ', cef_features_tax$Taxonomy)
cef_features_tax$Taxonomy <- gsub(',', ', ', cef_features_tax$Taxonomy)
strep_features_tax <- clean_merge(strep_features, taxonomy_otu)
strep_features_tax <- strep_features_tax[order(strep_features_tax[,1]),] 
strep_features_tax$Taxonomy <- gsub('_', ' ', strep_features_tax$Taxonomy)
strep_features_tax$Taxonomy <- gsub(',', ', ', strep_features_tax$Taxonomy)
clinda_features_tax <- clean_merge(clinda_features, taxonomy_otu)
clinda_features_tax <- clinda_features_tax[order(clinda_features_tax[,1]),] 
clinda_features_tax$Taxonomy <- gsub('_', ' ', clinda_features_tax$Taxonomy)
clinda_features_tax$Taxonomy <- gsub(',', ', ', clinda_features_tax$Taxonomy)
rm(taxonomy_otu, cef_features, strep_features, clinda_features)

# Get significant taxa from shared file
cef_shared_sig_otu <- subset(cef_shared_otu, rownames(cef_shared_otu) == )





# Phylotype family-level shared file
metadata_shared_family <- clean_merge(metadata, shared_family)
conv_shared_family <- subset(metadata_shared_family, type == 'conventional')
conv_shared_family$type <- NULL
conv_shared_family$infection <- NULL
conv_shared_family$abx <- NULL
rm(shared_family, metadata_shared_family)

# Convert to relative abundance
relabund_shared <- conv_shared_family + 1
relabund_shared <- (relabund_shared / rowSums(relabund_shared)) * 100
rm(conv_shared_family)

# Bin lowly abundant OTUs into an 'Other' category
relabund_shared[relabund_shared < 1] <- 0
relabund_shared <- relabund_shared[, colSums(relabund_shared != 0) > 0]
relabund_shared$Other <- 100 - rowSums(relabund_shared)

# NMDS axes
metadata_nmds <- clean_merge(metadata, nmds)
metadata_nmds$type <- NULL
control_axes <- subset(metadata_nmds, abx == 'none')
control_axes$abx <- NULL
control_axes$infection <- NULL
cef_axes <- subset(metadata_nmds, abx == 'cefoperazone')
cef_axes$abx <- NULL
cef_infected_axes <- subset(cef_axes, infection == '630')
cef_infected_axes$infection <- NULL
cef_mock_axes <- subset(cef_axes, infection == 'mock')
cef_mock_axes$infection <- NULL
clinda_axes <- subset(metadata_nmds, abx == 'clindamycin')
clinda_axes$abx <- NULL
clinda_infected_axes <- subset(clinda_axes, infection == '630')
clinda_infected_axes$infection <- NULL
clinda_mock_axes <- subset(clinda_axes, infection == 'mock')
clinda_mock_axes$infection <- NULL
strep_axes <- subset(metadata_nmds, abx == 'streptomycin')
strep_axes$abx <- NULL
strep_infected_axes <- subset(strep_axes, infection == '630')
strep_infected_axes$infection <- NULL
strep_mock_axes <- subset(strep_axes, infection == 'mock')
strep_mock_axes$infection <- NULL
rm(cef_axes, clinda_axes, strep_axes, metadata_nmds, metadata, nmds)

#----------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_2.pdf'
pdf(file=plot_file, width=12, height=7)
layout(matrix(c(1,2,2,
                4,5,6),
              nrow=2, ncol=3, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# NMDS of treatment groups

jpeg(filename='~/Desktop/nmds.jpeg', height = 900, width = 900)
par(las=1, mar=c(5,5,1,1))
plot(0,type='n', pch=16, cex=0,
     xlim=c(-0.7,0.7), ylim=c(-0.7,0.7), cex.lab=2, cex.axis=1.7,
     xlab='NMDS Axis 1', ylab='NMDS Axis 2')

# Add points
points(x=control_axes$axis1, y=control_axes$axis2, col='black', pch=1, lwd=4, cex=3.5)
points(x=cef_infected_axes$axis1, y=cef_infected_axes$axis2, col=wes_palette("FantasticFox")[3], pch=1, lwd=4, cex=3.5)
points(x=cef_mock_axes$axis1, y=cef_mock_axes$axis2, col=wes_palette("FantasticFox")[3], pch=6, lwd=4, cex=3)
points(x=clinda_infected_axes$axis1, y=clinda_infected_axes$axis2, col=wes_palette("FantasticFox")[5], pch=1, lwd=4, cex=3.5)
points(x=clinda_mock_axes$axis1, y=clinda_mock_axes$axis2, col=wes_palette("FantasticFox")[5], pch=6, lwd=4, cex=3)
points(x=strep_infected_axes$axis1, y=strep_infected_axes$axis2, col=wes_palette("FantasticFox")[1], pch=1, lwd=4, cex=3.5)
points(x=strep_mock_axes$axis1, y=strep_mock_axes$axis2, col=wes_palette("FantasticFox")[1], pch=6, lwd=4, cex=3)

# Add legends
legend('topleft', legend=c('Streptomycin-treated', 'Cefoperzone-treated', 'Clindamycin-treated', 'No Antibiotics'), 
       col=c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'black'), 
       pch=15, cex=1.9, pt.cex=3.5, bty='n')
legend('bottomleft', legend=c('Mock Infected', '630 Infected'), 
       col='black', pch=c(16,17), cex=2.5, pt.cex=c(3.5,3), bty='n')


dev.off()

mtext('A', side=2, line=2, las=2, adj=1.7, padj=-18.1, cex=1.3)

#-----------------------#

# Family-level phylotype bar chart

par(mar=c(4,5,1,1))

# When needed, use this pallete   
final_colors <- c("gold1", "orangered1", "aquamarine3", "firebrick", 
                  "forestgreen", "blue3", "mediumorchid2", "violetred4", 
                  "mediumpurple4", "dodgerblue3", "goldenrod3", "chartreuse3")


# Plot the final formatted table
barplot(t(relabund_shared), col=final_colors, yaxt='n', 
        ylim=c(0,100), ylab='% Relative Abundance', font=2)
box()
axis(side=2, at=seq(0,100,20), tick=TRUE)
segments(x0=rep(0,4), y0=seq(20,80,20), x1=rep(5,4), y1=seq(20,80,20), lty=2)

# Create a figure legend in the margin
legend(5.025, 75, legend=taxa, pt.bg=final_colors, pch=22, pt.cex=1.3, cex=0.7)

mtext('B', side=2, line=2, las=2, adj=1.7, padj=-10.5, cex=1.3)

#-----------------------#

# Random Forest results


pdf(file='~/Desktop/random_forest.pdf', width=20, height=8)
layout(matrix(c(1,2,3), nrow=1, ncol=3, byrow = TRUE))
# Cefoperazone plot



pdf(file='~/Desktop/cef_rf.pdf', width=8, height=8)
par(mar=c(5,3,1,1))
dotchart(cef_features_tax$final_features_RF, labels=cef_features_tax$Taxonomy, pch=19, cex=1.2,
         xlab='Mean Decrease Accuracy', xlim=c(3,7))
dev.off()


#mtext('C', side=2, line=2, las=2, adj=1.7, padj=-10.5, cex=1.3)

# Clindamycin plot
pdf(file='~/Desktop/clinda_rf.pdf', width=7.5, height=8)
par(mar=c(5,3,1,1))
dotchart(clinda_features_tax$final_features_RF, labels=clinda_features_tax$Taxonomy, pch=19, cex=1.2,
         xlab='Mean Decrease Accuracy', xlim=c(2,12))
dev.off()



# Streptomycin plot
pdf(file='~/Desktop/strep_rf.pdf', width=8, height=8)
par(mar=c(5,3,1,1))
dotchart(strep_features_tax$final_features_RF, labels=strep_features_tax$Taxonomy, pch=19, cex=1.2,
         xlab='Mean Decrease Accuracy', xlim=c(2,10))
dev.off()







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