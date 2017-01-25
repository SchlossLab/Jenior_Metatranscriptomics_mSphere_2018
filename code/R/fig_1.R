
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
shared_otu_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'
taxonomy_otu_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.cons.genus.format.taxonomy'
shared_family_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.family.subsample.shared'
taxonomy_family_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.family.cons.family.format.taxonomy'
wetlab_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/wetlab_assays.tsv'
metadata_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metadata.tsv'
otu_blast_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/'

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
taxonomy_family <- read.delim(taxonomy_family_file, sep='\t', header=T)
taxonomy_family$Size <- NULL
shared_family <- read.delim(shared_family_file, sep='\t', header=T, row.names=2)
shared_family <- shared_family[!rownames(shared_family) %in% c('CefC5M2'), ]  # Remove contaminated sample
shared_family$numOtus <- NULL
shared_family$label <- NULL
wetlab <- read.delim(wetlab_file, sep='\t', header=T, row.names=1)

rm(shared_otu_file, taxonomy_otu_file, shared_family_file, taxonomy_family_file, metadata_file, wetlab_file)

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
   cef_shared_otu, clinda_shared_otu, strep_shared_otu, taxonomy_otu)

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

#--------------------------#

# Phylotype family-level shared file

# Convert to relative abundance
relabund_shared <- (shared_family / rowSums(shared_family)) * 100
rm(shared_family)

# Bin lowly abundant OTUs into an 'Other' category
relabund_shared[relabund_shared < 2] <- 0
relabund_shared <- relabund_shared[, colSums(relabund_shared != 0) > 0]
top_otus <- colnames(relabund_shared)

# Subset family-level taxonomy
rownames(taxonomy_family) <- taxonomy_family$OTU
taxonomy_family <- subset(taxonomy_family, rownames(taxonomy_family) %in% top_otus)
taxonomy_family$Taxonomy <- gsub('_', ' ', taxonomy_family$Taxonomy)
taxonomy_family$Taxonomy <- gsub(',', ', ', taxonomy_family$Taxonomy)
taxonomy_family$Taxonomy <- gsub('\\(OTU\\d+\\)', '', taxonomy_family$Taxonomy)
taxonomy_family <- taxonomy_family[order(taxonomy_family$Taxonomy),]
rm(top_otus)

# Sort shared based on OTUs and add 'Other' category
relabund_shared <- as.data.frame(t(relabund_shared))
relabund_shared$otu <- rownames(relabund_shared)
relabund_shared <- relabund_shared[match(as.vector(taxonomy_family$OTU), relabund_shared$otu),]
relabund_shared$otu <- NULL
relabund_shared <- as.data.frame(t(relabund_shared))
relabund_shared$Other <- 100 - rowSums(relabund_shared)
taxonomy_family <- as.vector(taxonomy_family$Taxonomy)
taxonomy_family <- append(taxonomy_family, 'Other')

# Add empty columns for plotting and sort table
relabund_shared <- clean_merge(metadata, relabund_shared)
relabund_shared$type <- NULL
relabund_shared$infection <- NULL
relabund_shared$susceptibility <- NULL
empty_columns <- as.data.frame(rbind(c('col_1', rep(0, ncol(relabund_shared)-1)),
                                     c('col_2', rep(0, ncol(relabund_shared)-1)),
                                     c('col_3', rep(0, ncol(relabund_shared)-1))))
colnames(empty_columns) <- colnames(relabund_shared)
relabund_shared <- rbind(relabund_shared, empty_columns)
relabund_shared <- relabund_shared[order(match(relabund_shared$abx, c('none', 'col_1', 'streptomycin', 'col_2', 'cefoperazone', 'col_3', 'clindamycin'))),] # sort by treatment
relabund_shared$abx <- NULL
rm(empty_columns)


# Format taxonomy legend
family_colors <- c('forestgreen','chartreuse2',
                   'darkred','firebrick3','red2','tomato2','salmon',
                   'navy','mediumblue','royalblue3','dodgerblue2','lightslateblue','deepskyblue','powderblue',
                   'gold1',
                   'chocolate1',
                   'gray40')
taxonomy_color <- as.data.frame(cbind(taxonomy_family, family_colors))
taxonomy_color <- taxonomy_color[order(taxonomy_color$taxonomy_family),]
other <- taxonomy_color[which(taxonomy_color$taxonomy_family == 'Other'),]
taxonomy_color <- subset(taxonomy_color, taxonomy_family != 'Other')
taxonomy_color <- rbind(other, taxonomy_color)
colnames(taxonomy_color) <- c('taxonomy', 'color')
rm(other)

#--------------------------#

# Format wetlab assay data
wetlab <- subset(wetlab, infection == '630') # Remove uninfected controls
wetlab$infection <- NULL
wetlab <- subset(wetlab, treatment != 'germfree')
wetlab$cfu_vegetative <- as.numeric(wetlab$cfu_vegetative)
wetlab$cfu_vegetative[wetlab$cfu_vegetative == 0] <- 100
wetlab$cfu_vegetative <- log10(wetlab$cfu_vegetative)
wetlab$cfu_spore <- as.numeric(wetlab$cfu_spore)
wetlab$cfu_spore[wetlab$cfu_spore == 0] <- 100
wetlab$cfu_spore <- log10(wetlab$cfu_spore)
wetlab$treatment <- factor(wetlab$treatment, levels=c('conventional', 'streptomycin', 'cefoperazone', 'clindamycin'))

wetlab$cfu_vegetative[wetlab$cfu_vegetative <= 2.0] <- 1.7 # Undetectable points below LOD
wetlab$cfu_spore[wetlab$cfu_spore <= 2.0] <- 1.7 # Undetectable points below LOD
wetlab$toxin_titer[wetlab$toxin_titer <= 2.0] <- 1.94 # Undetectable points below LOD

#----------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/figures/figure_1.pdf'
pdf(file=plot_file, width=12, height=10)
layout(matrix(c(1,2,2,
                3,3,4,
                5,6,7),
              nrow=3, ncol=3, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Create an empty plot
par(mar=c(0,0,0,0))
plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-4.75,4), ylim=c(-2,5))

# Abx in drinking water timeline
rect(xleft=-4, ybottom=2.8, xright=0, ytop=3.2, col='cadetblue3', border='black')
Arrows(x0=-4, y0=3, x1=3.5, y1=3, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,0,2,2.75), y0=c(3.4,3.4,3.4,3.4), x1=c(-4,0,2,2.75), y1=c(2.6,2.6,2.6,2.6), lwd=4, col=c('black','black','black','gray40'))
segments(x0=c(-4,-3,-2,-1,1), y0=c(3.25,3.25,3.25,3.25,3.25), x1=c(-4,-3,-2,-1,1), y1=c(2.75,2.75,2.75,2.75,2.75), lwd=2)
points(x=c(2,2.75), y=c(3.8,3.8), pch=25, bg=c('white','black'), col='black', cex=2.5)
text(x=c(-4,0,2), y=c(2.1,2.1,2.1), c('Day -7', 'Day -2', 'Day 0'), cex=1.1)
text(x=-3.3, y=3.7, 'Streptomycin', cex=1.1)
text(x=-2.1, y=3.72, 'or', cex=1.2, font=2)
text(x=-0.8, y=3.7, 'Cefoperazone', cex=1.1)
text(-4.4, 3, '1', font=2, cex=1.5)

# IP injection abx timeline
Arrows(x0=-4, y0=0, x1=-1.5, y1=0, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,-3,-2.25), y0=c(-0.4,-0.4,-0.4), x1=c(-4,-3,-2.25), y1=c(0.4,0.4,0.4), lwd=4, col=c('black','black','gray40'))
points(x=c(-4,-3,-2.25), y=c(0.8,0.8,0.8), pch=c(25,25,25), bg=c('coral1','white','black'), col='black', cex=2.5)
text(x=c(-4,-3), y=c(-0.8,-0.8), c('Day -1', 'Day 0'), cex=1.1)
text(-4.4, 0, '2', font=2, cex=1.5)

# Legend
legend(x=-0.6, y=1, legend=expression('Antibiotic in Drinking Water', 'Clindamycin IP Injection',paste(italic('C. difficile'), ' Spore Gavage'), 'Necropsy (18 hours)'), 
       pt.bg=c('cadetblue3','coral1','white','black'), pch=c(22,25,25,25), cex=1.2, pt.cex=c(3,2,2,2), bty='n')

# Plot label
text(-4.7, 4.88, 'A', cex=2)

#-------------------------------------------------------------------#

# CFU and toxin data

# Vegetative cells
par(mar=c(3,4,1,4), mgp=c(2.5, 1, 0))
stripchart(cfu_vegetative~treatment, data=wetlab, col='black', bg='firebrick2', xlim=c(0,22), ylim=c(0,9), pch=21,
           vertical=TRUE, at=c(0.5, 6.5, 12.5, 18.5), xaxt='n', yaxt='n', ylab='CFU/g Cecal Content', cex=1.7, method='jitter', jitter=0.2)
abline(h=2, lwd=1.5, col='gray25') # LOD
abline(v=c(5,11,17), lty=2, lwd=1.5) # dividers
axis(side=2, at=seq(0,9,1), labels=c(0, parse(text=paste(rep(10,9), '^', seq(1,9,1), sep=''))), las=1)
axis(side=1, at=c(2,8,14,20), cex.axis=1.1, tick=FALSE,
     labels=c('No Antibiotics','Streptomycin-treated','Cefoperazone-treated','Clindamycin-treated'))

# Median lines
segments(x0=c(0.5, 6.5, 12.5, 18.5)-0.6, y0=c(
  as.numeric(median(wetlab[wetlab$treatment == 'conventional', 2])),
  as.numeric(median(wetlab[wetlab$treatment == 'streptomycin', 2])),
  as.numeric(median(wetlab[wetlab$treatment == 'cefoperazone', 2])),
  as.numeric(median(wetlab[wetlab$treatment == 'clindamycin', 2]))),
x1=c(0.5, 6.5, 12.5, 18.5)+0.6, y1=c(
  as.numeric(median(wetlab[wetlab$treatment == 'conventional', 2])),
  as.numeric(median(wetlab[wetlab$treatment == 'streptomycin', 2])),
  as.numeric(median(wetlab[wetlab$treatment == 'cefoperazone', 2])),
  as.numeric(median(wetlab[wetlab$treatment == 'clindamycin', 2]))), lwd=3)

# Spores
stripchart(cfu_spore~treatment, data=wetlab, col='black', bg='blue2', xlim=c(0,22), ylim=c(0,9), pch=21,
           vertical=TRUE, at=c(2, 8, 14, 20), xaxt='n', yaxt='n', ylab='', cex=1.7, method='jitter', jitter=0.2, add=TRUE)
# Median lines
segments(x0=c(2, 8, 14, 20)-0.6, y0=c(
  as.numeric(median(wetlab[wetlab$treatment == 'conventional', 3])),
  as.numeric(median(wetlab[wetlab$treatment == 'streptomycin', 3])),
  as.numeric(median(wetlab[wetlab$treatment == 'cefoperazone', 3])),
  as.numeric(median(wetlab[wetlab$treatment == 'clindamycin', 3]))),
  x1=c(2, 8, 14, 20)+0.6, y1=c(
    as.numeric(median(wetlab[wetlab$treatment == 'conventional', 3])),
    as.numeric(median(wetlab[wetlab$treatment == 'streptomycin', 3])),
    as.numeric(median(wetlab[wetlab$treatment == 'cefoperazone', 3])),
    as.numeric(median(wetlab[wetlab$treatment == 'clindamycin', 3]))), lwd=3)

# Toxin
par(mar=c(3,4,1,4), new=TRUE, xpd=TRUE)
stripchart(toxin_titer~treatment, data=wetlab, col='black', bg='green3', xlim=c(0,22), ylim=c(1.6,3.4), pch=21,
           vertical=TRUE, at=c(3.5, 9.5, 15.5, 21.5), xaxt='n', yaxt='n', ylab='', cex=1.7, method='jitter', jitter=0.2)
# Median lines
segments(x0=c(3.5, 9.5, 15.5, 21.5)-0.6, y0=c(
  as.numeric(median(wetlab[wetlab$treatment == 'conventional', 4])),
  as.numeric(median(wetlab[wetlab$treatment == 'streptomycin', 4])),
  as.numeric(median(wetlab[wetlab$treatment == 'cefoperazone', 4])),
  as.numeric(median(wetlab[wetlab$treatment == 'clindamycin', 4]))),
  x1=c(3.5, 9.5, 15.5, 21.5)+0.6, y1=c(
    as.numeric(median(wetlab[wetlab$treatment == 'conventional', 4])),
    as.numeric(median(wetlab[wetlab$treatment == 'streptomycin', 4])),
    as.numeric(median(wetlab[wetlab$treatment == 'cefoperazone', 4])),
    as.numeric(median(wetlab[wetlab$treatment == 'clindamycin', 4]))), lwd=3)
axis(side=4, at=seq(1.6,3.4,0.2), las=1,
     labels=c('1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2','3.4'))
mtext(expression(paste('Toxin Titer/g Cecal Content (',log[10],')')), side=4, line=3, cex=0.7)

legend('topleft', legend=c('Vegetative cells (CFU)','Spores (CFU)'), bty='n',
       pch=21, col='black', pt.bg=c('firebrick2','blue2'), cex=1.1, pt.cex=2)
legend('topright', legend='Toxin titer', bty='n',
       pch=21, col='black', pt.bg='green3', cex=1.1, pt.cex=2)

mtext('B', side=2, line=2, las=2, adj=1.7, padj=-8, cex=1.3)

#-------------------------------------------------------------------#

# Family-level phylotype bar chart

par(mar=c(5,5,1,1), new=FALSE, xpd=FALSE)
barplot(t(relabund_shared), col=family_colors, yaxt='n', xaxt='n', 
        ylim=c(0,100), ylab='% Relative Abundance', cex.names=1.2)
box()
axis(side=2, at=seq(0,100,20), tick=TRUE, las=1)
abline(h=c(20,40,60,80), lty=2)

# Lines under plot
arrows(0.5, -2, 6.9, -2, angle=0, length=0, lwd=2, xpd=TRUE)
arrows(9.1, -2, 18.9, -2, angle=0, length=0, lwd=2, xpd=TRUE)
arrows(19.9, -2, 29.7, -2, angle=0, length=0, lwd=2, xpd=TRUE)
arrows(32, -2, 40.3, -2, angle=0, length=0, lwd=2, xpd=TRUE)
arrows(41.3, -2, 51.3, -2, angle=0, length=0, lwd=2, xpd=TRUE)
arrows(53.4, -2, 63.1, -2, angle=0, length=0, lwd=2, xpd=TRUE)
arrows(64.4, -2, 74, -2, angle=0, length=0, lwd=2, xpd=TRUE)

mtext(rep('630 infected',4), side=1, at=c(3.7,14,36,58.3), adj=0.5, padj=1, cex=0.7)
mtext(rep('Mock infected',3), side=1, at=c(24.7,46.5,69.1), adj=0.5, padj=1, cex=0.7)
mtext(c('No Antibiotics','Streptomycin-treated','Cefoperazone-treated','Clindamycin-treated'), side=1, at=c(3.7,19.3,41.5,63.7), adj=0.5, padj=2.5, cex=0.75)
mtext('C', side=2, line=2, las=2, adj=1.7, padj=-8, cex=1.3)

# Create a figure legend in empty plot
par(mar=c(4,0,0,1))
plot(0, type='n', ylim=c(-5,5), xlim=c(5,5), ylab='', xlab='', xaxt='n', yaxt='n', axes=FALSE)
legend('center', legend=as.vector(taxonomy_color$taxonomy), pt.bg=as.vector(taxonomy_color$color), pch=22, pt.cex=2.5, cex=1.1)

#-------------------------------------------------------------------#

# Random Forest results - Significant OTUs and relative abundance changes

# Streptomycin plot
par(mar=c(4, 15, 2, 1), mgp=c(2.3, 1, 0), xpd=FALSE)
plot(1, type="n", ylim=c(0,length(strep_otus)*2), xlim=c(0,4), 
     ylab="", xlab="Normalized Abundance", xaxt="n", yaxt="n") # make blank plot
title('Streptomycin-treated', line=0.5, cex.main=1.1, font.main=1)
index <- 1
for(i in strep_otus){
  stripchart(at=index-0.35, jitter(strep_mock_otu[,i], amount=1e-5), 
             pch=21, bg="skyblue1", method="jitter", jitter=0.12, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(strep_infected_otu[,i], amount=1e-5), 
             pch=21, bg="tomato", method="jitter", jitter=0.12, add=T, cex=1, lwd=0.5)
  if (i != strep_otus[length(strep_otus)]){
    abline(h=index+1, lty=2)
  }
  segments(median(strep_mock_otu[,i]), index-0.7, median(strep_mock_otu[,i]), index, lwd=2) #adds line for median
  segments(median(strep_infected_otu[,i]), index+0.7, median(strep_infected_otu[,i]), index, lwd=2)
  index_p <- wilcox.test(strep_mock_otu[,i], strep_infected_otu[,i], exact=FALSE)$p.value
  if (index_p <= 0.05){
    text(ceiling(max(c(as.vector(strep_mock_otu[,i]), as.vector(strep_infected_otu[,i])))), index, labels='*', cex=2)
  }
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=strep_otus, las=1, line=-0.5, tick=F, cex.axis=0.8) 
axis(1, at=c(0, 1, 2, 3, 4), label=c('0','10', '100', "1000", "10000"))
legend('topright', legend=c("630 infected", "Mock infected"), 
       pch=c(21, 21), pt.bg=c("tomato","skyblue1"), bg='white', pt.cex=1.5, cex=0.9)

mtext('D', side=2, line=2, las=2, adj=9.8, padj=-8, cex=1.3)

#-----------------#

# Cefoperazone plot
par(mar=c(4, 13, 2, 1), mgp=c(2.3, 1, 0))
plot(1, type='n', ylim=c(0,length(cef_otus)*2), xlim=c(0,4), 
     ylab='', xlab='Normalized Abundance', xaxt='n', yaxt='n') # make blank plot
title('Cefoperazone-treated', line=0.5, cex.main=1.1, font.main=1)
index <- 1
for(i in cef_otus){
  stripchart(at=index-0.35, jitter(cef_mock_otu[,i], amount=1e-5), 
             pch=21, bg="skyblue1", method="jitter", jitter=0.12, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(cef_infected_otu[,i], amount=1e-5), 
             pch=21, bg="tomato", method="jitter", jitter=0.12, add=T, cex=1, lwd=0.5)
  if (i != cef_otus[length(cef_otus)]){
  abline(h=index+1, lty=2)
  }
  segments(median(cef_mock_otu[,i]), index-0.7, median(cef_mock_otu[,i]), index, lwd=2) #adds line for median
  segments(median(cef_infected_otu[,i]), index+0.7, median(cef_infected_otu[,i]), index, lwd=2)
  index_p <- wilcox.test(cef_mock_otu[,i], cef_infected_otu[,i], exact=FALSE)$p.value
  if (index_p <= 0.05){
    text(ceiling(max(c(as.vector(cef_mock_otu[,i]), as.vector(cef_infected_otu[,i])))), index, labels='*', cex=2)
  }
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=cef_otus, las=1, line=-0.5, tick=F, cex.axis=0.8) 
axis(1, at=c(0, 1, 2, 3, 4), label=c('0','10', '100', "1000", '10000'))
legend('topright', legend=c("630 infected", "Mock infected"), 
       pch=c(21, 21), pt.bg=c("tomato","skyblue1"), bg='white', pt.cex=1.5, cex=0.9)

#-----------------#

# Clindamycin plot
par(mar=c(4, 12, 2, 1), mgp=c(2.3, 1, 0))
plot(1, type="n", ylim=c(0,length(clinda_otus)*2), xlim=c(0,3), 
     ylab="", xlab="Normalized Abundance", xaxt="n", yaxt="n") # make blank plot
title('Clindamycin-treated', line=0.5, cex.main=1.1, font.main=1)
index <- 1
for(i in clinda_otus){
  stripchart(at=index-0.35, jitter(clinda_mock_otu[,i], amount=1e-5), 
             pch=21, bg="skyblue1", method="jitter", jitter=0.12, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(clinda_infected_otu[,i], amount=1e-5), 
             pch=21, bg="tomato", method="jitter", jitter=0.12, add=T, cex=1, lwd=0.5)
  if (i != clinda_otus[length(clinda_otus)]){
    abline(h=index+1, lty=2)
  }
  segments(median(clinda_mock_otu[,i]), index-0.7, median(clinda_mock_otu[,i]), index, lwd=2) #adds line for median
  segments(median(clinda_infected_otu[,i]), index+0.7, median(clinda_infected_otu[,i]), index, lwd=2)
  index_p <- wilcox.test(clinda_mock_otu[,i], clinda_infected_otu[,i], exact=FALSE)$p.value
  if (index_p <= 0.05){
    text(ceiling(max(c(as.vector(clinda_mock_otu[,i]), as.vector(clinda_infected_otu[,i])))), index, labels='*', cex=2)
  }
  index <- index + 2
}
axis(2, at=seq(1,index-2,2), labels=clinda_otus, las=1, line=-0.5, tick=F, cex.axis=0.8) 
axis(1, at=c(0, 1, 2, 3), label=c('0','10', '100', "1000"))
legend('topright', legend=c("630 infected", "Mock infected"), 
       pch=c(21, 21), pt.bg=c("tomato","skyblue1"), bg='white', pt.cex=1.5, cex=0.9)

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



