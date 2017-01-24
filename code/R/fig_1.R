
# Load dependencies
deps <- c('wesanderson', 'randomForest', 'vegan', 'shape')
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
cfu <- read.delim(cfu_file, sep='\t', header=T)
toxin <- read.delim(toxin_file, sep='\t', header=T)

rm(shared_otu_file, taxonomy_otu_file, shared_family_file, taxonomy_family_file, metadata_file, cfu_file, toxin_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format data (previously subsampled and filtered in mothur)

# OTU shared file
metadata_shared_otu <- clean_merge(metadata, shared_otu)
#res_shared_otu <- metadata_shared_otu
#res_shared_otu$infection <- NULL
#res_shared_otu$abx <- NULL
#abx_shared_otu <- subset(metadata_shared_otu, abx != 'none')
#abx_shared_otu$infection <- NULL
#abx_shared_otu$susceptibility <- NULL
cef_shared_otu <- subset(metadata_shared_otu, abx == 'cefoperazone')
cef_shared_otu$abx <- NULL
cef_shared_otu$susceptibility <- NULL
strep_shared_otu <- subset(metadata_shared_otu, abx == 'streptomycin')
strep_shared_otu$abx <- NULL
strep_shared_otu$susceptibility <- NULL
clinda_shared_otu <- subset(metadata_shared_otu, abx == 'clindamycin')
clinda_shared_otu$abx <- NULL
clinda_shared_otu$susceptibility <- NULL
rm(metadata_shared_otu, shared_otu)

# Run random forest
cef_features <- featureselect_RF(cef_shared_otu, 'infection')
strep_features <- featureselect_RF(strep_shared_otu, 'infection')
clinda_features <- featureselect_RF(clinda_shared_otu, 'infection')

# Filter for significant OTUs with random forest from shared
cef_feat_shared <- cef_shared_otu[,rownames(cef_features)]
cef_feat_tax <- clean_merge(t(cef_feat_shared), taxonomy_otu)
cef_feat_tax$Taxonomy <- gsub('_', ' ', cef_feat_tax$Taxonomy)
cef_feat_tax$Taxonomy <- gsub(',', ', ', cef_feat_tax$Taxonomy)
rownames(cef_feat_tax) <- cef_feat_tax$Taxonomy
cef_feat_tax$Taxonomy <- NULL
cef_feat_tax <- as.data.frame(t(cef_feat_tax))
cef_feat_tax <- log10(cef_feat_tax + 1)
cef_feat_tax$infection <- cef_shared_otu$infection
clinda_feat_shared <- clinda_shared_otu[,rownames(clinda_features)]
clinda_feat_tax <- clean_merge(t(clinda_feat_shared), taxonomy_otu)
clinda_feat_tax$Taxonomy <- gsub('_', ' ', clinda_feat_tax$Taxonomy)
clinda_feat_tax$Taxonomy <- gsub(',', ', ', clinda_feat_tax$Taxonomy)
rownames(clinda_feat_tax) <- clinda_feat_tax$Taxonomy
clinda_feat_tax$Taxonomy <- NULL
clinda_feat_tax <- as.data.frame(t(clinda_feat_tax))
clinda_feat_tax <- log10(clinda_feat_tax + 1)
clinda_feat_tax$infection <- clinda_shared_otu$infection
strep_feat_shared <- strep_shared_otu[,rownames(strep_features)]
strep_feat_tax <- clean_merge(t(strep_feat_shared), taxonomy_otu)
strep_feat_tax$Taxonomy <- gsub('_', ' ', strep_feat_tax$Taxonomy)
strep_feat_tax$Taxonomy <- gsub(',', ', ', strep_feat_tax$Taxonomy)
rownames(strep_feat_tax) <- strep_feat_tax$Taxonomy
strep_feat_tax$Taxonomy <- NULL
strep_feat_tax <- as.data.frame(t(strep_feat_tax))
strep_feat_tax <- log10(strep_feat_tax + 1)
strep_feat_tax$infection <- strep_shared_otu$infection
rm(cef_features, strep_features, clinda_features, 
   cef_feat_shared, clinda_feat_shared, strep_feat_shared,
   cef_shared_otu, clinda_shared_otu, strep_shared_otu)

# Subset to experimental groups
cef_infected_otu <- subset(cef_feat_tax, infection == '630')
cef_infected_otu$infection <- NULL
cef_mock_otu <- subset(cef_feat_tax, infection == 'mock')
cef_mock_otu$infection <- NULL
rm(cef_feat_tax)
clinda_infected_otu <- subset(clinda_feat_tax, infection == '630')
clinda_infected_otu$infection <- NULL
clinda_mock_otu <- subset(clinda_feat_tax, infection == 'mock')
clinda_mock_otu$infection <- NULL
rm(clinda_feat_tax)
strep_infected_otu <- subset(strep_feat_tax, infection == '630')
strep_infected_otu$infection <- NULL
strep_mock_otu <- subset(strep_feat_tax, infection == 'mock')
strep_mock_otu$infection <- NULL
rm(strep_feat_tax)
cef_otus <- colnames(cef_infected_otu)
clinda_otus <- colnames(clinda_infected_otu)
strep_otus <- colnames(strep_infected_otu)

# Phylotype family-level shared file
conv_shared_family <- clean_merge(metadata, shared_family)
conv_shared_family$type <- NULL
conv_shared_family$infection <- NULL
conv_shared_family$abx <- NULL
conv_shared_family$susceptibility <- NULL
rm(shared_family)

# Convert to relative abundance
relabund_shared <- (conv_shared_family / rowSums(conv_shared_family)) * 100
rm(conv_shared_family)

# Bin lowly abundant OTUs into an 'Other' category
relabund_shared[relabund_shared < 1] <- 0
relabund_shared <- relabund_shared[, colSums(relabund_shared != 0) > 0]
relabund_shared$Other <- 100 - rowSums(relabund_shared)

# Format CFU data and collect summary statistics
cfu[cfu == 0] <- 100
cfu$cfu_vegetative <- log10(cfu$cfu_vegetative)
cfu$cfu_spore <- log10(cfu$cfu_spore)
cfu$mouse <- NULL
cfu <- subset(cfu, cage < 4 ) # Remove uninfected controls
cfu$cage <- NULL
cfu$treatment <- factor(cfu$treatment, levels=c('streptomycin', 'cefoperazone', 'clindamycin', 'germfree', 'conventional'))
vegetative_cfu <- cfu
vegetative_cfu$cfu_spore <- NULL
spore_cfu <- cfu
spore_cfu$cfu_vegetative <- NULL
cef <- as.numeric(median(vegetative_cfu[vegetative_cfu$treatment == 'cefoperazone', 2]))
strep <- as.numeric(median(vegetative_cfu[vegetative_cfu$treatment == 'streptomycin', 2]))
clinda <- as.numeric(median(vegetative_cfu[vegetative_cfu$treatment == 'clindamycin', 2]))
gf <- as.numeric(median(vegetative_cfu[vegetative_cfu$treatment == 'germfree', 2]))
conv <- as.numeric(median(vegetative_cfu[vegetative_cfu$treatment == 'conventional', 2]))
vege_medians <- c(strep, cef, clinda, gf, conv)
vege_medians[vege_medians == 2.0] <- 1.6
cef <- as.numeric(median(spore_cfu[spore_cfu$treatment == 'cefoperazone', 2]))
strep <- as.numeric(median(spore_cfu[spore_cfu$treatment == 'streptomycin', 2]))
clinda <- as.numeric(median(spore_cfu[spore_cfu$treatment == 'clindamycin', 2]))
gf <- as.numeric(median(spore_cfu[spore_cfu$treatment == 'germfree', 2]))
conv <- as.numeric(median(spore_cfu[spore_cfu$treatment == 'conventional', 2]))
spore_medians <- c(strep, cef, clinda, gf, conv)
spore_medians[spore_medians == 2.0] <- 1.6
rm(cfu, cef, strep, clinda, gf, conv)
vegetative_cfu$color <- ifelse(vegetative_cfu$cfu_vegetative == 2.0, 'gray50', 'black')
vegetative_cfu$cfu_vegetative[vegetative_cfu$cfu_vegetative == 2.0] <- 1.6
spore_cfu$cfu_spore[spore_cfu$cfu_spore == 2.0] <- 1.6

# Format toxin data and find summary statistics
toxin$mouse <- NULL
toxin$cage <- NULL
toxin$treatment <- factor(toxin$treatment, levels=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'Germfree', 'Conventional'))
cef <- as.numeric(median(toxin[toxin$treatment == 'Cefoperazone', 2]))
strep <- as.numeric(median(toxin[toxin$treatment == 'Streptomycin', 2]))
clinda <- as.numeric(median(toxin[toxin$treatment == 'Clindamycin', 2]))
gf <- as.numeric(median(toxin[toxin$treatment == 'Germfree', 2]))
conv <- as.numeric(median(toxin[toxin$treatment == 'Conventional', 2]))
toxin_medians <- c(strep, cef, clinda, gf, conv)
toxin_medians[toxin_medians <= 2.0] <- 1.9
toxin$titer[toxin$titer <= 2.0] <- 1.9
rm(cef, strep, clinda, gf, conv)

#----------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/figures/figure_1.pdf'
pdf(file=plot_file, width=12, height=7)
layout(matrix(c(1,2,2,
                3,4,5),
              nrow=2, ncol=3, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Create an empty plot
par(mar=c(0,0,0,0))
plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-4.75,4), ylim=c(-2,5))

# Abx in drinking water timeline
rect(xleft=-4, ybottom=2.8, xright=0, ytop=3.2, col='red', border='black')
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
points(x=c(-4,-3,-2.25), y=c(1,1,1), pch=c(25,25,25), bg=c('blue','white','black'), col='black', cex=2.5)
text(x=c(-4,-3,-2.25), y=c(-0.8,-0.8,-0.8), c('Day -1', 'Day 0', '18 hrs'), cex=0.9)
text(x=-4.5, y=0, 'Clindamycin', cex=0.8)

# Legend
legend(x=0, y=0.7, legend=expression('Antibiotic in Drinking Water', 'IP Injection of Antibiotic',paste(italic('C. difficile'), ' Spore Gavage'), 'Sacrifice & Necropsy'), 
       pt.bg=c('red','blue','white','black'), pch=c(22,25,25,25), cex=1.2, pt.cex=c(3.2,2.2,2.2,2.2), bty='n')

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
legend(5.025, 75, legend=taxonomy_family, pt.bg=final_colors, pch=22, pt.cex=1.3, cex=0.7)

mtext('B', side=2, line=2, las=2, adj=1.7, padj=-10.5, cex=1.3)

#-------------------------------------------------------------------#

# Random Forest results - Significant OTUs and relative abundance changes

# Cefoperazone plot
par(mar=c(5, 15, 1, 1))

plot(1, type="n", ylim=c(0,length(cef_otus)*2), xlim=c(0,4), 
     ylab="", xlab="Normalized Abundance", xaxt="n", yaxt="n") # make blank plot
index <- 1
for(i in cef_otus){
  stripchart(at=index-0.35, jitter(cef_mock_otu[,i], amount=1e-5), 
             pch=21, bg="royalblue1", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(cef_infected_otu[,i], amount=1e-5), 
             pch=21, bg="red", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  segments(median(cef_mock_otu[,i]), index-0.7, median(cef_mock_otu[,i]), index, lwd=3) #adds line for median
  segments(median(cef_infected_otu[,i]), index+0.7, median(cef_infected_otu[,i]), index, lwd=3)
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=cef_otus, las=1, line=-0.5, tick=F, cex.axis=0.8) 
axis(1, at=c(0, 1, 2, 3, 4), label=c('0','10', '100', "1000", '10000'))
legend('topright', legend=c("630 Infected", "Mock Infected"), 
       pch=c(21, 21), pt.bg=c("red","royalblue1"), pt.cex=2, cex=1.2)

mtext('C', side=2, line=2, las=2, adj=1.7, padj=-10.5, cex=1.3)

#-----------------#

# Clindamycin plot
par(mar=c(5, 14, 1, 1))
plot(1, type="n", ylim=c(0,length(clinda_otus)*2), xlim=c(0,3), 
     ylab="", xlab="Normalized Abundance", xaxt="n", yaxt="n") # make blank plot
index <- 1
for(i in clinda_otus){
  stripchart(at=index-0.35, jitter(clinda_mock_otu[,i], amount=1e-5), 
             pch=21, bg="royalblue1", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(clinda_infected_otu[,i], amount=1e-5), 
             pch=21, bg="red", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  segments(median(clinda_mock_otu[,i]), index-0.7, median(clinda_mock_otu[,i]), index, lwd=3) #adds line for median
  segments(median(clinda_infected_otu[,i]), index+0.7, median(clinda_infected_otu[,i]), index, lwd=3)
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
  stripchart(at=index-0.35, jitter(strep_mock_otu[,i], amount=1e-5), 
             pch=21, bg="royalblue1", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(strep_infected_otu[,i], amount=1e-5), 
             pch=21, bg="red", method="jitter", jitter=0.2, add=T, cex=1, lwd=0.5)
  segments(median(strep_mock_otu[,i]), index-0.7, median(strep_mock_otu[,i]), index, lwd=3) #adds line for median
  segments(median(strep_infected_otu[,i]), index+0.7, median(strep_infected_otu[,i]), index, lwd=3)
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=strep_otus, las=1, line=-0.5, tick=F, cex.axis=0.8) 
axis(1, at=c(0, 1, 2, 3, 4), label=c('0','10', '100', "1000", "10000"))
legend('topright', legend=c("630 Infected", "Mock Infected"), 
       pch=c(21, 21), pt.bg=c("red","royalblue1"), pt.cex=2, cex=1.2)

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



