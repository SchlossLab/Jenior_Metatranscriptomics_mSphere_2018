
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
  subfactor_1 <- round(length(rownames(training_data[which(training_data[,feature]==levels[1]),])) * 0.623)
  subfactor_2 <- round(length(rownames(training_data[which(training_data[,feature]==levels[2]),])) * 0.623)
  factor <- max(c(round(subfactor_1 / subfactor_2), round(subfactor_2 / subfactor_1))) * 3
  n_tree <- round(length(colnames(training_data)) - 1) * factor
  m_try <- round(sqrt(length(colnames(training_data)) - 1))
  
  # Run random forest
  data_randomForest <- randomForest(training_data[,feature]~., data=training_data, importance=TRUE, replace=FALSE, do.trace=500, err.rate=TRUE, ntree=n_tree, mtry=m_try)
  detach(training_data)
  
  # Parse features for significance and sort
  features_RF <- importance(data_randomForest, type=1)
  final_features_RF <- subset(features_RF, features_RF > abs(min(features_RF)))
  final_features_RF <- final_features_RF[!(rownames(final_features_RF) == feature),]
  final_features_RF <- as.data.frame(final_features_RF)

  return(final_features_RF)
}

# Merge data frames with shared row names
clean_merge <- function(data_1, data_2){
  
  clean_merged <- merge(data_1, data_2, by = 'row.names')
  rownames(clean_merged) <- clean_merged$Row.names
  clean_merged$Row.names <- NULL
  
  return(clean_merged)
}

# Select files
nmds_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.unique_list.conventional.thetayc.0.03.lt.ave.nmds.axes'
shared_otu_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'
taxonomy_otu_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.cons.genus.format.taxonomy'
shared_family_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.family.subsample.shared'
taxonomy_family_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.family.cons.family.format.taxonomy'
cfu_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/cfu.dat'
toxin_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/toxin_titer.dat'
metadata_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metadata.tsv'

#-------------------------------------------------------------------------------------------------------------------------------------#

# Load in data
metadata <- read.delim(metadata_file, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2'), ] # Remove contaminated sample
metadata <- subset(metadata, type == 'conventional') # remove germfree
metadata$type <- NULL
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

# Format data (previously subsampled and filtered in mothur)

# OTU shared file
metadata_shared_otu <- clean_merge(metadata, shared_otu)
res_shared_otu <- metadata_shared_otu
res_shared_otu$infection <- NULL
res_shared_otu$abx <- NULL
abx_shared_otu <- subset(metadata_shared_otu, abx != 'none')
abx_shared_otu$infection <- NULL
abx_shared_otu$susceptibility <- NULL
cef_shared_otu <- subset(metadata_shared_otu, abx == 'cefoperazone')
cef_shared_otu$abx <- NULL
cef_shared_otu$susceptibility <- NULL
strep_shared_otu <- subset(metadata_shared_otu, abx == 'streptomycin')
strep_shared_otu$abx <- NULL
strep_shared_otu$susceptibility <- NULL
clinda_shared_otu <- subset(metadata_shared_otu, abx == 'clindamycin')
clinda_shared_otu$abx <- NULL
clinda_shared_otu$susceptibility <- NULL
rm(shared_otu)

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
rm(cef_features, strep_features, clinda_features)

# Subset shared file to significant OTUs from random forest




# Transform OTU abundances
cef_infected_otu <- subset(cef_shared_otu, infection == '630')
cef_infected_otu$infection <- NULL
cef_infected_otu <- log10(cef_infected_otu + 1)
cef_mock_otu <- subset(cef_shared_otu, infection == 'mock')
cef_mock_otu$infection <- NULL
cef_mock_otu <- log10(cef_mock_otu + 1)
rm(cef_shared_otu)
clinda_infected_otu <- subset(clinda_shared_otu, infection == '630')
clinda_infected_otu$infection <- NULL
clinda_infected_otu <- log10(clinda_infected_otu + 1)
clinda_mock_otu <- subset(clinda_shared_otu, infection == 'mock')
clinda_mock_otu$infection <- NULL
clinda_mock_otu <- log10(clinda_mock_otu + 1)
rm(clinda_shared_otu)
strep_infected_otu <- subset(strep_shared_otu, infection == '630')
strep_infected_otu$infection <- NULL
strep_infected_otu <- log10(strep_infected_otu + 1)
strep_mock_otu <- subset(strep_shared_otu, infection == 'mock')
strep_mock_otu$infection <- NULL
strep_mock_otu <- log10(strep_mock_otu + 1)
rm(strep_shared_otu)





# Integrate taxa into OTUs
cef_infected_otu <- t(clean_merge(t(cef_infected_otu), taxonomy_otu))
cef_mock_otu <- t(clean_merge(t(cef_mock_otu), taxonomy_otu))
clinda_infected_otu <- t(clean_merge(t(clinda_infected_otu), taxonomy_otu))
clinda_mock_otu <- t(clean_merge(t(clinda_mock_otu), taxonomy_otu))
strep_infected_otu <- t(clean_merge(t(strep_infected_otu), taxonomy_otu))
strep_mock_otu <- t(clean_merge(t(strep_mock_otu), taxonomy_otu))
rm(taxonomy_otu)

# Phylotype family-level shared file
metadata_shared_family <- clean_merge(metadata, shared_family)
conv_shared_family <- subset(metadata_shared_family, type == 'conventional')
conv_shared_family$type <- NULL
conv_shared_family$infection <- NULL
conv_shared_family$abx <- NULL
rm(shared_family, metadata_shared_family)

# Convert to relative abundance
relabund_shared <- (conv_shared_family / rowSums(conv_shared_family)) * 100
rm(conv_shared_family)

# Bin lowly abundant OTUs into an 'Other' category
relabund_shared[relabund_shared < 1] <- 0
relabund_shared <- relabund_shared[, colSums(relabund_shared != 0) > 0]
relabund_shared$Other <- 100 - rowSums(relabund_shared)

# NMDS axes
#metadata_nmds <- clean_merge(metadata, nmds)
#metadata_nmds$type <- NULL
#control_axes <- subset(metadata_nmds, abx == 'none')
#control_axes$abx <- NULL
#control_axes$infection <- NULL
#cef_axes <- subset(metadata_nmds, abx == 'cefoperazone')
#cef_axes$abx <- NULL
#cef_infected_axes <- subset(cef_axes, infection == '630')
#cef_infected_axes$infection <- NULL
#cef_mock_axes <- subset(cef_axes, infection == 'mock')
#cef_mock_axes$infection <- NULL
#clinda_axes <- subset(metadata_nmds, abx == 'clindamycin')
#clinda_axes$abx <- NULL
#clinda_infected_axes <- subset(clinda_axes, infection == '630')
#clinda_infected_axes$infection <- NULL
#clinda_mock_axes <- subset(clinda_axes, infection == 'mock')
#clinda_mock_axes$infection <- NULL
#strep_axes <- subset(metadata_nmds, abx == 'streptomycin')
#strep_axes$abx <- NULL
#strep_infected_axes <- subset(strep_axes, infection == '630')
#strep_infected_axes$infection <- NULL
#strep_mock_axes <- subset(strep_axes, infection == 'mock')
#strep_mock_axes$infection <- NULL
#rm(cef_axes, clinda_axes, strep_axes, metadata_nmds, metadata, nmds)











#----------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_1.pdf'
pdf(file=plot_file, width=12, height=7)
layout(matrix(c(1,2,2,
                4,5,6),
              nrow=2, ncol=3, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# NMDS of treatment groups

#par(las=1, mar=c(5,5,1,1))
#plot(0,type='n', pch=16, cex=0,
#     xlim=c(-0.7,0.7), ylim=c(-0.7,0.7), cex.lab=2, cex.axis=1.7,
#     xlab='NMDS Axis 1', ylab='NMDS Axis 2')
#
## Add points
#points(x=control_axes$axis1, y=control_axes$axis2, col='black', pch=1, lwd=4, cex=3.5)
#points(x=cef_infected_axes$axis1, y=cef_infected_axes$axis2, col=wes_palette("FantasticFox")[3], pch=1, lwd=4, cex=3.5)
#points(x=cef_mock_axes$axis1, y=cef_mock_axes$axis2, col=wes_palette("FantasticFox")[3], pch=6, lwd=4, cex=3)
#points(x=clinda_infected_axes$axis1, y=clinda_infected_axes$axis2, col=wes_palette("FantasticFox")[5], pch=1, lwd=4, cex=3.5)
#points(x=clinda_mock_axes$axis1, y=clinda_mock_axes$axis2, col=wes_palette("FantasticFox")[5], pch=6, lwd=4, cex=3)
#points(x=strep_infected_axes$axis1, y=strep_infected_axes$axis2, col=wes_palette("FantasticFox")[1], pch=1, lwd=4, cex=3.5)
#points(x=strep_mock_axes$axis1, y=strep_mock_axes$axis2, col=wes_palette("FantasticFox")[1], pch=6, lwd=4, cex=3)
#
## Add legends
#legend('topleft', legend=c('Streptomycin-treated', 'Cefoperzone-treated', 'Clindamycin-treated', 'No Antibiotics'), 
#       col=c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'black'), 
#       pch=15, cex=1.9, pt.cex=3.5, bty='n')
#legend('bottomleft', legend=c('Mock Infected', '630 Infected'), 
#       col='black', pch=c(16,17), cex=2.5, pt.cex=c(3.5,3), bty='n')
#
#
#dev.off()



# Create an empty plot
par(mar=c(0,0,0,0))
plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-4.75,4), ylim=c(-2,5))

# Abx in drinking water timeline
rect(xleft=-4, ybottom=2.8, xright=0, ytop=3.2, col=fox[3], border='black')
Arrows(x0=-4, y0=3, x1=3.5, y1=3, lwd=4, arr.type='triangle', arr.length=0.75, arr.width=0.4)
segments(x0=c(-4,0,2,2.75), y0=c(3.5,3.5,3.5,3.5), x1=c(-4,0,2,2.75), y1=c(2.5,2.5,2.5,2.5), lwd=4)
segments(x0=c(-4,-3,-2,-1,1), y0=c(3.25,3.25,3.25,3.25,3.25), x1=c(-4,-3,-2,-1,1), y1=c(2.75,2.75,2.75,2.75,2.75), lwd=2)
points(x=c(2,2.75), y=c(4,4), pch=25, bg=c('white','black'), col='black', cex=2.5)
text(x=c(-4,0,2,2.75), y=c(2.2,2.2,2.2,2.2), c('Day -7', 'Day -2', 'Day 0', '18 hrs'), cex=0.9)
text(x=-4.5, y=3.2, 'Cefoperazone', cex=0.8)
text(x=-4.5, y=2.95, 'or', font=2)
text(x=-4.5, y=2.7, 'Streptomycin', cex=0.8)

# IP injection abx timeline
Arrows(x0=-4, y0=0, x1=-1.5, y1=0, lwd=4, arr.type='triangle', arr.length=0.75, arr.width=0.4)
segments(x0=c(-4,-3,-2.25), y0=c(-0.5,-0.5,-0.5), x1=c(-4,-3,-2.25), y1=c(0.5,0.5,0.5), lwd=4)
points(x=c(-4,-3,-2.25), y=c(1,1,1), pch=c(25,25,25), bg=c(fox[5],'white','black'), col='black', cex=2.5)
text(x=c(-4,-3,-2.25), y=c(-0.8,-0.8,-0.8), c('Day -1', 'Day 0', '18 hrs'), cex=0.9)
text(x=-4.5, y=0, 'Clindamycin', cex=0.8)

# Legend
legend(x=0, y=0.7, legend=expression('Antibiotic in Drinking Water', 'IP Injection of Antibiotic',paste(italic('C. difficile'), ' Spore Gavage'), 'Sacrifice & Necropsy'), 
       pt.bg=c(fox[3],fox[5],'white','black'), pch=c(22,25,25,25), cex=1.2, pt.cex=c(3.2,2.2,2.2,2.2), bty='n')

# Plot label
legend('topleft', legend='A', cex=2, bty='n')


#-------------------------------------------------------------------#

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

#-------------------------------------------------------------------#

# Random Forest results - Ssignificant OTUs and relative abundance changes

# Cefoperazone plot
par(mar=c(5, 15, 1, 1))
plot(1, type="n", ylim=c(0,length(cef_otus)*2), xlim=c(0,4), 
     ylab="", xlab="Normalized Abundance", xaxt="n", yaxt="n") # make blank plot
index <- 1
for(i in cef_otus){
  stripchart(at=index-0.35, jitter(cef_mock[,i], amount=1e-5), 
             pch=21, bg="royalblue1", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(cef_infected[,i], amount=1e-5), 
             pch=21, bg="red", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  segments(median(cef_mock[,i]), index-0.7, median(cef_mock[,i]), index, lwd=3) #adds line for median
  segments(median(cef_infected[,i]), index+0.7, median(cef_infected[,i]), index, lwd=3)
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=cef_otus, las=1, line=-0.5, tick=F, cex.axis=0.8) 
axis(1, at=c(0, 1, 2, 3, 4), label=c('0','10', '100', "1000", '10000'))
legend('topright', legend=c("630 Infected", "Mock Infected"), 
       pch=c(21, 21), pt.bg=c("red","royalblue1"), pt.cex=2, cex=1.2)

#-----------------#

# Clindamycin plot
par(mar=c(5, 14, 1, 1))
plot(1, type="n", ylim=c(0,length(clinda_otus)*2), xlim=c(0,3), 
     ylab="", xlab="Normalized Abundance", xaxt="n", yaxt="n") # make blank plot
index <- 1
for(i in clinda_otus){
  stripchart(at=index-0.35, jitter(clinda_mock[,i], amount=1e-5), 
             pch=21, bg="royalblue1", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(clinda_infected[,i], amount=1e-5), 
             pch=21, bg="red", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  segments(median(clinda_mock[,i]), index-0.7, median(clinda_mock[,i]), index, lwd=3) #adds line for median
  segments(median(clinda_infected[,i]), index+0.7, median(clinda_infected[,i]), index, lwd=3)
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=clinda_otus, las=1, line=-0.5, tick=F, cex.axis=0.8) 
axis(1, at=c(0, 1, 2, 3), label=c('0','10', '100', "1000"))
legend('topright', legend=c("630 Infected", "Mock Infected"), 
       pch=c(21, 21), pt.bg=c("red","royalblue1"), pt.cex=2, cex=1.2)

#-----------------#

# Streptomycin plot
par(mar=c(5, 15, 1, 1))
plot(1, type="n", ylim=c(0,length(strep_otus)*2), xlim=c(0,4), 
     ylab="", xlab="Normalized Abundance", xaxt="n", yaxt="n") # make blank plot
index <- 1
for(i in strep_otus){
  stripchart(at=index-0.35, jitter(strep_mock[,i], amount=1e-5), 
             pch=21, bg="royalblue1", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(strep_infected[,i], amount=1e-5), 
             pch=21, bg="red", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  segments(median(strep_mock[,i]), index-0.7, median(strep_mock[,i]), index, lwd=3) #adds line for median
  segments(median(strep_infected[,i]), index+0.7, median(strep_infected[,i]), index, lwd=3)
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=strep_otus, las=1, line=-0.5, tick=F, cex.axis=0.8) 
axis(1, at=c(0, 1, 2, 3, 4), label=c('0','10', '100', "1000", "10000"))
legend('topright', legend=c("630 Infected", "Mock Infected"), 
       pch=c(21, 21), pt.bg=c("red","royalblue1"), pt.cex=2, cex=1.2)




mtext('C', side=2, line=2, las=2, adj=1.7, padj=-10.5, cex=1.3)


dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
#rm(plot_file, metadata_axes, mock_axes, infected_axes, metadata_summary, 
#   strep_div, cef_div, clinda_div, conv_div)
#for (dep in deps){
#  pkg <- paste('package:', dep, sep='')
#  detach(pkg, character.only = TRUE)
#}
#rm(dep, deps, pkg)
#gc()



