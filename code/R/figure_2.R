
# Load dependencies
deps <- c('wesanderson', 'randomForest', 'vegan');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}
set.seed(6189)


classify_RF <- function(training_data, feature){
  
  # Set parameters
  attach(training_data)
  levels <- as.vector(unique(training_data[,feature]))
  factor_1 <- round(length(rownames(data[which(data[,feature]==levels[1]),])) * 0.623)
  factor_2 <- round(length(rownames(data[which(data[,feature]==levels[2]),])) * 0.623)
  factor <- max(c(round(factor_1 / factor_2), round(factor_2 / factor_1))) * 3 # Segal et al. (2004)
  n_tree <- round(length(colnames(training_data)) - 1) * factor
  m_try <- round(sqrt(length(colnames(training_data)) - 1))
  
  # Run random forest
  data_randomForest <- randomForest(training_data[,feature]~., data=training_data, importance=TRUE, replace=FALSE, do.trace=500, err.rate=TRUE, ntree=n_tree, mtry=m_try)
  detach(training_data)
  
  # Parse features for significance
  features_RF <- importance(data_randomForest, type=1)
  final_features_RF <- subset(features_RF, features_RF>abs(min(features_RF))) # Segal et al. (2004)
  
  return(final_features_RF)
  
}

clean_merge <- function(data_1, data_2){
  
  clean_merged <- merge(data_1, data_2, by = 'row.names')
  rownames(clean_merged) <- clean_merged$Row.names
  clean_merged$Row.names <- NULL
  
  return(clean_merged)
  
}


# Define variables
nmds_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.unique_list.thetayc.0.03.lt.ave.nmds.axes'
shared_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'
taxonomy_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.family.cons.taxonomy' # wrong taxonomy
metadata_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metadata.tsv'

# Load in data
nmds <- read.delim(nmds_file, sep='\t', header=T, row.names=1)
nmds <- nmds[!rownames(nmds) %in% c('CefC5M2'), ] # Remove contaminated sample
metadata <- read.delim(metadata_file, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2'), ] # Remove contaminated sample
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
taxonomy <- read.delim(taxonomy_file, sep='\t', header=T, row.names=1)
taxonomy$Size <- NULL
shared <- read.delim(shared_file, sep='\t', header=T, row.names=2)
shared <- shared[!rownames(shared) %in% c('CefC5M2'), ]  # Remove contaminated sample
shared$numOtus <- NULL
shared$label <- NULL

rm(nmds_file, taxonomy_file, metadata_file, shared_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format data

# NMDS axes
metadata_axes <- clean_merge(metadata, nmds)
mock_axes <- subset(metadata_axes, infection == 'mock')
mock_axes <- subset(mock_axes, type != 'germfree')
infected_axes <- subset(metadata_axes, infection == '630')
infected_axes <- subset(infected_axes, type != 'germfree')
rm(nmds, metadata_axes)

# Shared file
metadata_shared <- clean_merge(metadata, shared)
conventional_shared <- subset(metadata_shared, type == 'conventional')
conventional_shared$type <- NULL
conventional_shared$infection <- NULL
conventional_tax <- t(clean_merge(taxonomy, t(conventional_shared)))
cef_shared <- subset(metadata_shared, abx == 'cefoperazone')
cef_shared$abx <- NULL
cef_shared$type <- NULL
strep_shared <- subset(metadata_shared, abx == 'streptomycin')
strep_shared$abx <- NULL
strep_shared$type <- NULL
clinda_shared <- subset(metadata_shared, abx == 'clindamycin')
clinda_shared$abx <- NULL
clinda_shared$type <- NULL
rm(shared, metadata, metadata_shared, taxonomy)

# Run random forest
cef_features <- classify_RF(cef_shared, 'infection')
strep_features <- classify_RF(strep_shared, 'infection')
clinda_features <- classify_RF(clinda_shared, 'infection')

#-------------------------------------------------------------------------------------------------------------------------------------#

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
       col=c(wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen', 'black', wes_palette("FantasticFox")[1])[mock_axes$abx], 
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