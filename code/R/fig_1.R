
# Load dependencies
deps <- c('vegan', 'shape', 'Matrix')
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Load in functions
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define files
shared_otu_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'
shared_family_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.family.subsample.shared'
taxonomy_family_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.family.cons.family.format.taxonomy'
wetlab_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/wetlab_assays.tsv'
metadata_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metadata.tsv'
otu_tax_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/formatted.all_treatments.0.03.cons.taxonomy'

#-------------------------------------------------------------------------------------------------------------------------------------#

# Load in data
metadata <- read.delim(metadata_file, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2'), ] # Remove contaminated sample
metadata <- subset(metadata, type == 'conventional') # remove germfree
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
shared_otu <- read.delim(shared_otu_file, sep='\t', header=T, row.names=2)
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('CefC5M2'), ]  # Remove contaminated sample
shared_otu <- shared_otu[ , !(names(shared_otu) == 'Otu0004')] # Remove residual C. difficile OTU
shared_otu$numOtus <- NULL
shared_otu$label <- NULL
taxonomy_family <- read.delim(taxonomy_family_file, sep='\t', header=T)
taxonomy_family$Size <- NULL
shared_family <- read.delim(shared_family_file, sep='\t', header=T, row.names=2)
shared_family <- shared_family[!rownames(shared_family) %in% c('CefC5M2'), ]  # Remove contaminated sample
shared_family <- shared_family[ , !(names(shared_family) == 'Otu008')] # Remove residual C. difficile OTU
shared_family$numOtus <- NULL
shared_family$label <- NULL
wetlab <- read.delim(wetlab_file, sep='\t', header=T, row.names=1)
otu_tax <- read.delim(otu_tax_file, sep='\t', header=T, row.names=1)
rm(shared_otu_file, shared_family_file, taxonomy_family_file, metadata_file, wetlab_file, otu_tax_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format data (previously subsampled and filtered in mothur)

# OTU shared file
metadata_shared_otu <- clean_merge(metadata, shared_otu)
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

# Filter OTUs and transform
cef_shared_otu <- filter_table(cef_shared_otu)
cef_shared_otu[,2:ncol(cef_shared_otu)] <- log10(cef_shared_otu[,2:ncol(cef_shared_otu)] + 1)
strep_shared_otu <- filter_table(strep_shared_otu)
strep_shared_otu[,2:ncol(strep_shared_otu)] <- log10(strep_shared_otu[,2:ncol(strep_shared_otu)] + 1)
clinda_shared_otu <- filter_table(clinda_shared_otu)
clinda_shared_otu[,2:ncol(clinda_shared_otu)] <- log10(clinda_shared_otu[,2:ncol(clinda_shared_otu)] + 1)




# Filter for significant OTUs by Wilcoxon
cef_infected_otu <- subset(cef_shared_otu, infection == '630')
cef_infected_otu$infection <- NULL
cef_mock_otu <- subset(cef_shared_otu, infection == 'mock')
cef_mock_otu$infection <- NULL
cef_pvalues <- c()
index2 <- 1
for (index in colnames(cef_infected_otu)){
  if (wilcox.test(cef_mock_otu[,index], cef_infected_otu[,index], exact=FALSE)$p.value > 0.05){
    cef_infected_otu <- cef_infected_otu[, !(names(cef_infected_otu) == index)]
    cef_mock_otu <- cef_mock_otu[, !(names(cef_mock_otu) == index)]
    next
  }
  else if (median(cef_mock_otu[,index]) == median(cef_infected_otu[,index])){
    cef_infected_otu <- cef_infected_otu[, !(names(cef_infected_otu) == index)]
    cef_mock_otu <- cef_mock_otu[, !(names(cef_mock_otu) == index)]
    next
  }
  else {
    cef_pvalues[index2] <- wilcox.test(cef_mock_otu[,index], cef_infected_otu[,index], exact=FALSE)$p.value
    index2 <- index2 + 1
  }
}
cef_pvalues <- p.adjust(cef_pvalues, method='BH')
cef_final_p <- c()
for (index in 1:length(cef_pvalues)){
  if (cef_pvalues[index] > 0.05){
    cef_infected_otu[,index] <- NULL
    cef_mock_otu[,index] <- NULL
  }
  else {
    cef_final_p <- c(cef_final_p, cef_pvalues[index])
  }
}
cef_final_p <- as.character(round(cef_final_p, 3))
rm(cef_feat_shared, cef_features)

clinda_infected_otu <- subset(clinda_shared_otu, infection == '630')
clinda_infected_otu$infection <- NULL
clinda_mock_otu <- subset(clinda_shared_otu, infection == 'mock')
clinda_mock_otu$infection <- NULL
clinda_pvalues <- c()
index2 <- 1
for (index in colnames(clinda_infected_otu)){
  if (wilcox.test(clinda_mock_otu[,index], clinda_infected_otu[,index], exact=FALSE)$p.value > 0.05){
    clinda_infected_otu <- clinda_infected_otu[, !(names(clinda_infected_otu) == index)]
    clinda_mock_otu <- clinda_mock_otu[, !(names(clinda_mock_otu) == index)]
    next
  }
  else if (median(clinda_mock_otu[,index]) == median(clinda_infected_otu[,index])){
    clinda_infected_otu <- clinda_infected_otu[, !(names(clinda_infected_otu) == index)]
    clinda_mock_otu <- clinda_mock_otu[, !(names(clinda_mock_otu) == index)]
    next
  }
  else {
    clinda_pvalues[index2] <- wilcox.test(clinda_mock_otu[,index], clinda_infected_otu[,index], exact=FALSE)$p.value
    index2 <- index2 + 1
  }
}
clinda_pvalues <- p.adjust(clinda_pvalues, method='BH')
clinda_final_p <- c()
for (index in 1:length(clinda_pvalues)){
  if (cef_pvalues[index] > 0.05){
    clinda_infected_otu[,index] <- NULL
    clinda_mock_otu[,index] <- NULL
  }
  else {
    clinda_final_p <- c(clinda_final_p, clinda_pvalues[index])
  }
}
clinda_final_p <- as.character(round(clinda_final_p, 3))
rm(clinda_feat_shared, clinda_features)

strep_infected_otu <- subset(strep_shared_otu, infection == '630')
strep_infected_otu$infection <- NULL
strep_mock_otu <- subset(strep_shared_otu, infection == 'mock')
strep_mock_otu$infection <- NULL
strep_pvalues <- c()
index2 <- 1
for (index in colnames(strep_infected_otu)){
  if (wilcox.test(strep_mock_otu[,index], strep_infected_otu[,index], exact=FALSE)$p.value > 0.05){
    strep_infected_otu <- strep_infected_otu[, !(names(strep_infected_otu) == index)]
    strep_mock_otu <- strep_mock_otu[, !(names(strep_mock_otu) == index)]
    next
  }
  else if (median(strep_mock_otu[,index]) == median(strep_infected_otu[,index])){
    strep_infected_otu <- strep_infected_otu[, !(names(strep_infected_otu) == index)]
    strep_mock_otu <- strep_mock_otu[, !(names(strep_mock_otu) == index)]
    next
  }
  else {
    strep_pvalues[index2] <- wilcox.test(strep_mock_otu[,index], strep_infected_otu[,index], exact=FALSE)$p.value
    index2 <- index2 + 1
  }
}
strep_pvalues <- p.adjust(strep_pvalues, method='BH')
strep_final_p <- c()
for (index in 1:length(strep_pvalues)){
  if (strep_pvalues[index] > 0.05){
    strep_infected_otu[,index] <- NULL
    strep_mock_otu[,index] <- NULL
  }
  else {
    strep_final_p <- c(strep_final_p, strep_pvalues[index])
  }
}
strep_final_p <- as.character(round(strep_final_p, 3))
rm(strep_feat_shared, strep_features)

# Rename OTUs with species-level identifier
cef_infected_otu <- clean_merge(t(cef_infected_otu), otu_tax)
cef_genera <- as.vector(cef_infected_otu$genus)
cef_phyla <- as.vector(cef_infected_otu$phylum)
cef_otus <- as.vector(cef_infected_otu$OTU_short)
cef_infected_otu$genus <- NULL
cef_infected_otu$phylum <- NULL
cef_infected_otu$OTU_short <- NULL
cef_infected_otu <- as.data.frame(t(cef_infected_otu))
cef_mock_otu <- clean_merge(t(cef_mock_otu), otu_tax)
cef_mock_otu$genus <- NULL
cef_mock_otu$phylum <- NULL
cef_mock_otu$OTU_short <- NULL
cef_mock_otu <- as.data.frame(t(cef_mock_otu))
clinda_infected_otu <- clean_merge(t(clinda_infected_otu), otu_tax)
clinda_genera <- as.vector(clinda_infected_otu$genus)
clinda_phyla <- as.vector(clinda_infected_otu$phylum)
clinda_otus <- as.vector(clinda_infected_otu$OTU_short)
clinda_infected_otu$genus <- NULL
clinda_infected_otu$phylum <- NULL
clinda_infected_otu$OTU_short <- NULL
clinda_infected_otu <- as.data.frame(t(clinda_infected_otu))
clinda_mock_otu <- clean_merge(t(clinda_mock_otu), otu_tax)
clinda_mock_otu$genus <- NULL
clinda_mock_otu$phylum <- NULL
clinda_mock_otu$OTU_short <- NULL
clinda_mock_otu <- as.data.frame(t(clinda_mock_otu))
strep_infected_otu <- clean_merge(t(strep_infected_otu), otu_tax)
strep_genera <- as.vector(strep_infected_otu$genus)
strep_phyla <- as.vector(strep_infected_otu$phylum)
strep_otus <- as.vector(strep_infected_otu$OTU_short)
strep_infected_otu$genus <- NULL
strep_infected_otu$phylum <- NULL
strep_infected_otu$OTU_short <- NULL
strep_infected_otu <- as.data.frame(t(strep_infected_otu))
strep_mock_otu <- clean_merge(t(strep_mock_otu), otu_tax)
strep_mock_otu$genus <- NULL
strep_mock_otu$phylum <- NULL
strep_mock_otu$OTU_short <- NULL
strep_mock_otu <- as.data.frame(t(strep_mock_otu))
rm(otu_tax)

#--------------------------#

# Phylotype family-level shared file

# Convert to relative abundance
relabund_shared <- (shared_family / rowSums(shared_family)) * 100
rm(shared_family)

# Bin lowly abundant OTUs into an 'Other' category
relabund_shared[relabund_shared < 1] <- 0
relabund_shared <- relabund_shared[, colSums(relabund_shared != 0) > 0]
top_otus <- colnames(relabund_shared)

# Subset family-level taxonomy
rownames(taxonomy_family) <- taxonomy_family$OTU
taxonomy_family <- subset(taxonomy_family, rownames(taxonomy_family) %in% top_otus)
taxonomy_family$Taxonomy <- gsub('_', ' ', taxonomy_family$Taxonomy)
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
taxonomy_family <- append(taxonomy_family, 'Other (<1% each)')

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
family_colors <- c('chartreuse3',
                   'mediumorchid3',
                   'darkred','firebrick3','salmon','tomato2','red2',
                   'navy','mediumblue','royalblue3','lightslateblue','dodgerblue2','deepskyblue','powderblue',
                   'chocolate2','darkgoldenrod1',
                   'magenta2',
                   'gray50')
taxonomy_color <- as.data.frame(cbind(taxonomy_family, family_colors))
taxonomy_color <- taxonomy_color[order(taxonomy_color$taxonomy_family),]
bacteria <- taxonomy_color[which(taxonomy_color$taxonomy_family == 'Bacteria unclassified'),]
taxonomy_color <- subset(taxonomy_color, taxonomy_family != 'Bacteria unclassified')
other <- taxonomy_color[which(taxonomy_color$taxonomy_family == 'Other (<1% each)'),]
taxonomy_color <- subset(taxonomy_color, taxonomy_family != 'Other (<1% each)')
taxonomy_color <- rbind(bacteria, taxonomy_color)
taxonomy_color <- rbind(other, taxonomy_color)
colnames(taxonomy_color) <- c('taxonomy', 'color')
taxonomy_color$taxonomy <- gsub('.*,', '', taxonomy_color$taxonomy)
rm(metadata, bacteria, other)

#--------------------------#

# Format wetlab assay data
wetlab <- subset(wetlab, infection == '630') # Remove uninfected controls
wetlab$infection <- NULL
wetlab$cfu_vegetative <- as.numeric(wetlab$cfu_vegetative)
wetlab$cfu_vegetative[wetlab$cfu_vegetative == 0] <- 100
wetlab$cfu_vegetative <- log10(wetlab$cfu_vegetative)
wetlab$cfu_spore <- as.numeric(wetlab$cfu_spore)
wetlab$cfu_spore[wetlab$cfu_spore == 0] <- 100
wetlab$cfu_spore <- log10(wetlab$cfu_spore)
wetlab$treatment <- factor(wetlab$treatment, levels=c('germfree','conventional', 'streptomycin', 'cefoperazone', 'clindamycin'))

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
segments(x0=c(-4,0,2,2.75), y0=c(3.4,3.4,3.4,3.3), x1=c(-4,0,2,2.75), y1=c(2.6,2.6,2.6,2.7), 
         lwd=4, col=c('black','black','black','black'))
segments(x0=c(-4,-3,-2,-1,1), y0=c(3.25,3.25,3.25,3.25,3.25), x1=c(-4,-3,-2,-1,1), y1=c(2.75,2.75,2.75,2.75,2.75), lwd=2)
points(x=c(2,2.75), y=c(3.8,3.8), pch=25, bg=c('white','black'), col='black', cex=2.5)
text(x=c(-4,0,2), y=c(2.1,2.1,2.1), c('Day -7', 'Day -2', 'Day 0'), cex=1.1)
text(x=-3.3, y=3.7, 'Streptomycin', cex=1.1)
text(x=-2.1, y=3.72, 'or', cex=1.2, font=2)
text(x=-0.8, y=3.7, 'Cefoperazone', cex=1.1)
text(-4.4, 3, '1', font=2, cex=1.5)

# IP injection abx timeline
Arrows(x0=-4, y0=0, x1=-1.5, y1=0, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,-3,-2.25), y0=c(-0.4,-0.4,-0.3), x1=c(-4,-3,-2.25), y1=c(0.4,0.4,0.3), lwd=4, col=c('black','black','black'))
points(x=c(-4,-3,-2.25), y=c(0.8,0.8,0.8), pch=c(25,25,25), bg=c('coral1','white','black'), col='black', cex=2.5)
text(x=c(-4,-3), y=c(-0.8,-0.8), c('Day -1', 'Day 0'), cex=1.1)
text(-4.4, 0, '2', font=2, cex=1.5)

# Legend
legend(x=-0.6, y=1, legend=expression('Antibiotic in Drinking Water', 'Clindamycin IP Injection',paste(italic('C. difficile'), ' Spore Gavage'), 'Necropsy (18 hours-post)'), 
       pt.bg=c('cadetblue3','coral1','white','black'), pch=c(22,25,25,25), cex=1.2, pt.cex=c(3,2,2,2), bty='n')

# Plot label
text(-4.7, 4.88, 'A', cex=2)

#-------------------------------------------------------------------#

# CFU and toxin data

# Vegetative cells
par(mar=c(3,4,1,4), mgp=c(2.5, 1, 0))
stripchart(cfu_vegetative~treatment, data=wetlab, col='black', bg='firebrick2', xlim=c(0,28), ylim=c(0,9), pch=21,
           vertical=TRUE, at=c(0.5, 6.5, 12.5, 18.5, 24.5), xaxt='n', yaxt='n', ylab='CFU/g Cecal Content', 
           cex=1.7, method='jitter', jitter=0.2)
abline(h=2, lwd=1.5, col='gray25') # LOD
abline(v=c(5,11,17,23), lty=2, lwd=1.5) # dividers
axis(side=2, at=seq(0,9,1), labels=c(0, parse(text=paste(rep(10,9), '^', seq(1,9,1), sep=''))), las=1)
axis(side=1, at=c(2,8,14,20,26), cex.axis=1.1, tick=FALSE,
     labels=c('Gnotobiotic','No Antibiotics (SPF)','Streptomycin (SPF)','Cefoperazone (SPF)','Clindamycin (SPF)'))

# Median lines
segments(x0=c(0.5, 6.5, 12.5, 18.5, 24.5)-0.6, y0=c(
  as.numeric(median(wetlab[wetlab$treatment == 'germfree', 2])),
  as.numeric(median(wetlab[wetlab$treatment == 'conventional', 2])),
  as.numeric(median(wetlab[wetlab$treatment == 'streptomycin', 2])),
  as.numeric(median(wetlab[wetlab$treatment == 'cefoperazone', 2])),
  as.numeric(median(wetlab[wetlab$treatment == 'clindamycin', 2]))),
x1=c(0.5, 6.5, 12.5, 18.5, 24.5)+0.6, y1=c(
  as.numeric(median(wetlab[wetlab$treatment == 'germfree', 2])),
  as.numeric(median(wetlab[wetlab$treatment == 'conventional', 2])),
  as.numeric(median(wetlab[wetlab$treatment == 'streptomycin', 2])),
  as.numeric(median(wetlab[wetlab$treatment == 'cefoperazone', 2])),
  as.numeric(median(wetlab[wetlab$treatment == 'clindamycin', 2]))), lwd=3)

# Spores
stripchart(cfu_spore~treatment, data=wetlab, col='black', bg='blue2', xlim=c(0,28), ylim=c(0,9), pch=21,
           vertical=TRUE, at=c(2, 8, 14, 20, 26), xaxt='n', yaxt='n', ylab='', cex=1.7, method='jitter', jitter=0.2, add=TRUE)
# Median lines
segments(x0=c(2, 8, 14, 20, 26)-0.6, y0=c(
  as.numeric(median(wetlab[wetlab$treatment == 'germfree', 3])),
  as.numeric(median(wetlab[wetlab$treatment == 'conventional', 3])),
  as.numeric(median(wetlab[wetlab$treatment == 'streptomycin', 3])),
  as.numeric(median(wetlab[wetlab$treatment == 'cefoperazone', 3])),
  as.numeric(median(wetlab[wetlab$treatment == 'clindamycin', 3]))),
  x1=c(2, 8, 14, 20, 26)+0.6, y1=c(
    as.numeric(median(wetlab[wetlab$treatment == 'germfree', 3])),
    as.numeric(median(wetlab[wetlab$treatment == 'conventional', 3])),
    as.numeric(median(wetlab[wetlab$treatment == 'streptomycin', 3])),
    as.numeric(median(wetlab[wetlab$treatment == 'cefoperazone', 3])),
    as.numeric(median(wetlab[wetlab$treatment == 'clindamycin', 3]))), lwd=3)

# Toxin
par(mar=c(3,4,1,4), new=TRUE, xpd=TRUE)
stripchart(toxin_titer~treatment, data=wetlab, col='black', bg='green3', xlim=c(0,28), ylim=c(1.6,3.4), pch=23,
           vertical=TRUE, at=c(3.5, 9.5, 15.5, 21.5, 27.5), xaxt='n', yaxt='n', ylab='', cex=1.7, method='jitter', jitter=0.2)
# Median lines
segments(x0=c(3.5, 9.5, 15.5, 21.5, 27.5)-0.6, y0=c(
  as.numeric(median(wetlab[wetlab$treatment == 'germfree', 4])),
  as.numeric(median(wetlab[wetlab$treatment == 'conventional', 4])),
  as.numeric(median(wetlab[wetlab$treatment == 'streptomycin', 4])),
  as.numeric(median(wetlab[wetlab$treatment == 'cefoperazone', 4])),
  as.numeric(median(wetlab[wetlab$treatment == 'clindamycin', 4]))),
  x1=c(3.5, 9.5, 15.5, 21.5, 27.5)+0.6, y1=c(
    as.numeric(median(wetlab[wetlab$treatment == 'germfree', 4])),
    as.numeric(median(wetlab[wetlab$treatment == 'conventional', 4])),
    as.numeric(median(wetlab[wetlab$treatment == 'streptomycin', 4])),
    as.numeric(median(wetlab[wetlab$treatment == 'cefoperazone', 4])),
    as.numeric(median(wetlab[wetlab$treatment == 'clindamycin', 4]))), lwd=3)
axis(side=4, at=seq(1.6,3.4,0.2), las=1,
     labels=c('1.6','1.8','2.0','2.2','2.4','2.6','2.8','3.0','3.2','3.4'))
mtext(expression(paste('Toxin Titer/g Cecal Content (',log[10],')')), side=4, line=3, cex=0.7)

legend('topleft', legend=c('Vegetative cells (CFU)','Spores (CFU)'), ncol=2, box.lwd=0, box.col='white',
       pch=21, col='black', pt.bg=c('firebrick2','blue2'), pt.cex=1.9, bg='white')
legend('topright', legend='Toxin titer', bty='n',
       pch=23, col='black', pt.bg='green3', pt.cex=1.9)
box()
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
mtext(c('No Antibiotics','Streptomycin-pretreated','Cefoperazone-pretreated','Clindamycin-pretreated'), 
      side=1, at=c(3.7,19.3,41.5,63.7), adj=0.5, padj=2.5, cex=0.75)
mtext('C', side=2, line=2, las=2, adj=1.7, padj=-8, cex=1.3)

# Create a figure legend in empty plot
par(mar=c(4,0,0.3,5))
plot(0, type='n', ylim=c(-5,5), xlim=c(5,5), ylab='', xlab='', xaxt='n', yaxt='n', axes=FALSE)
legend('right', legend=as.vector(taxonomy_color$taxonomy), pt.bg=as.vector(taxonomy_color$color), 
       pch=22, pt.cex=2.5, cex=1.1, bty='n')

# Add in phylum classifications
segments(x0=c(4.8,4.8,4.8,4.8,4.8), x1=c(4.8,4.8,4.8,4.8,4.8), 
         y0=c(3.88,3.1,0.3, -3.6,-4.64), 
         y1=c(3.5, 0.9,-3.1,-4.3,-5), lwd=3)
segments(x0=c(4.4,4.4,4.28,4.43,4.48), x1=c(4.8,4.8,4.8,4.8,4.8), 
         y0=c(3.69,2,-1.4,-3.95,-4.82), 
         y1=c(3.69,2,-1.4,-3.95,-4.82), lwd=2)

text(x=c(3.8,3.8,3.8,3.8,3.8), y=c(3.69,2,-1.4,-3.95,-4.82), cex=1.2,
     labels=c('Actinobacteria', 'Bacteroidetes', 'Firmicutes', 'Proteobacteria', 'Verrucomicrobia'))

#-------------------------------------------------------------------#

# Feature selection results - Significant OTUs and relative abundance changes

# Streptomycin plot
par(mar=c(4, 16, 2, 1), mgp=c(2.3, 1, 0), xpd=FALSE)
plot(1, type="n", ylim=c(0,length(strep_otus)*2), xlim=c(0,4), 
     ylab="", xlab="Normalized Abundance", xaxt="n", yaxt="n") # make blank plot
title('Streptomycin-pretreated', line=0.5, cex.main=1.1, font.main=1)
index <- 1
p_values <- c()
maxes <- c()
for(i in colnames(strep_mock_otu)){
  stripchart(at=index-0.35, jitter(strep_mock_otu[,i], amount=1e-5), 
             pch=21, bg="mediumseagreen", method="jitter", jitter=0.12, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(strep_infected_otu[,i], amount=1e-5), 
             pch=21, bg="mediumorchid3", method="jitter", jitter=0.12, add=T, cex=1, lwd=0.5)
  if (i != colnames(strep_mock_otu)[length(colnames(strep_mock_otu))]){
    abline(h=index+1, lty=2)
  }
  segments(median(strep_mock_otu[,i]), index-0.5, median(strep_mock_otu[,i]), index, lwd=2) #adds line for median
  segments(median(strep_infected_otu[,i]), index+0.5, median(strep_infected_otu[,i]), index, lwd=2)
  maxes <- append(maxes, ceiling(max(c(as.vector(strep_mock_otu[,i]), as.vector(strep_infected_otu[,i])))))
  index <- index + 2
}
axis(1, at=c(0, 1, 2, 3, 4), label=c('0','10', '100', "1000", "10000"))
legend('topright', legend=c("630 infected", "Mock infected"), cex=0.8,
       pch=c(21, 21), pt.bg=c("mediumorchid3","mediumseagreen"), bg='white', pt.cex=1.5)
formatted <- lapply(1:length(strep_otus), function(i) bquote(atop(paste(.(strep_phyla[i]), '; ', italic(.(strep_genera[i])), sep=''),
                                                                paste(.(strep_otus[i]), italic(' p'),' = ',.(strep_final_p[i]), sep=''))))
axis(2, at=seq(1,index-2,2), labels=do.call(expression, formatted), las=1, line=-0.5, tick=F, cex.axis=0.9, font=3) 

mtext('D', side=2, line=2, las=2, adj=9.8, padj=-8, cex=1.3)

#-----------------#

# Cefoperazone plot
par(mar=c(4, 15, 2, 1), mgp=c(2.3, 1, 0))
plot(1, type='n', ylim=c(0,length(cef_otus)*2), xlim=c(0,4), 
     ylab='', xlab='Normalized Abundance', xaxt='n', yaxt='n') # make blank plot
title('Cefoperazone-pretreated', line=0.5, cex.main=1.1, font.main=1)
index <- 1
p_values <- c()
maxes <- c()
for(i in colnames(cef_mock_otu)){
  stripchart(at=index-0.35, jitter(cef_mock_otu[,i], amount=1e-5), 
             pch=21, bg="mediumseagreen", method="jitter", jitter=0.12, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(cef_infected_otu[,i], amount=1e-5), 
             pch=21, bg="mediumorchid3", method="jitter", jitter=0.12, add=T, cex=1, lwd=0.5)
  if (i != colnames(cef_mock_otu)[length(colnames(cef_mock_otu))]){
  abline(h=index+1, lty=2)
  }
  segments(median(cef_mock_otu[,i]), index-0.5, median(cef_mock_otu[,i]), index, lwd=2) #adds line for median
  segments(median(cef_infected_otu[,i]), index+0.5, median(cef_infected_otu[,i]), index, lwd=2)
  maxes <- append(maxes, ceiling(max(c(as.vector(cef_mock_otu[,i]), as.vector(cef_infected_otu[,i])))))
  index <- index + 2
}
axis(1, at=c(0, 1, 2, 3, 4), label=c('0','10', '100', "1000", '10000'))
legend('topright', legend=c("630 infected", "Mock infected"), cex=0.8,
       pch=c(21, 21), pt.bg=c("mediumorchid3","mediumseagreen"), bg='white', pt.cex=1.5)
formatted <- lapply(1:length(cef_otus), function(i) bquote(atop(paste(.(cef_phyla[i]), '; ', italic(.(cef_genera[i])), sep=''),
                                                                   paste(.(cef_otus[i]), italic(' p'),' = ',.(cef_final_p[i]), sep=''))))
axis(2, at=seq(1,index-2,2), labels=do.call(expression, formatted), las=1, line=-0.5, tick=F, cex.axis=0.9, font=3) 

#-----------------#

# Clindamycin plot
par(mar=c(4, 13, 2, 1), mgp=c(2.3, 1, 0))
plot(1, type="n", ylim=c(0,length(clinda_otus)*2), xlim=c(0,3), 
     ylab="", xlab="Normalized Abundance", xaxt="n", yaxt="n") # make blank plot
title('Clindamycin-pretreated', line=0.5, cex.main=1.1, font.main=1)
index <- 1
p_values <- c()
maxes <- c()
for(i in colnames(clinda_mock_otu)){
  stripchart(at=index-0.35, jitter(clinda_mock_otu[,i], amount=1e-5), 
             pch=21, bg="mediumseagreen", method="jitter", jitter=0.12, add=T, cex=1, lwd=0.5)
  stripchart(at=index+0.35, jitter(clinda_infected_otu[,i], amount=1e-5), 
             pch=21, bg="mediumorchid3", method="jitter", jitter=0.12, add=T, cex=1, lwd=0.5)
  if (i != colnames(clinda_mock_otu)[length(colnames(clinda_mock_otu))]){
    abline(h=index+1, lty=2)
  }
  segments(median(clinda_mock_otu[,i]), index-0.5, median(clinda_mock_otu[,i]), index, lwd=2) #adds line for median
  segments(median(clinda_infected_otu[,i]), index+0.5, median(clinda_infected_otu[,i]), index, lwd=2)
  maxes <- append(maxes, ceiling(max(c(as.vector(clinda_mock_otu[,i]), as.vector(clinda_infected_otu[,i])))))
  index <- index + 2
}
axis(1, at=c(0, 1, 2, 3), label=c('0','10', '100', "1000"))
legend('topright', legend=c("630 infected", "Mock infected"), cex=0.8,
       pch=c(21, 21), pt.bg=c("mediumorchid3","mediumseagreen"), bg='white', pt.cex=1.5)
formatted <- lapply(1:length(clinda_otus), function(i) bquote(atop(paste(.(clinda_phyla[i]), '; ', italic(.(clinda_genera[i])), sep=''),
                                                              paste(.(clinda_otus[i]), italic(' p'),' = ',.(clinda_final_p[i]), sep=''))))
axis(2, at=seq(1,index-2,2), labels=do.call(expression, formatted), las=1, line=-0.5, tick=F, cex.axis=0.9, font=3) 

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up

for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
#rm(list=ls())
gc()


