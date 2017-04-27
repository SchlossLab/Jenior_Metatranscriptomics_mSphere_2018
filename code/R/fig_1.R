
# Start with clear environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define input files
shared_otu_file <- 'data/16S_analysis/all_treatments.0.03.unique_list.conventional.shared'
otu_tax_file <- 'data/16S_analysis/formatted.all_treatments.0.03.cons.taxonomy'
shared_family_file <- 'data/16S_analysis/all_treatments.family.subsample.shared'
taxonomy_family_file <- 'data/16S_analysis/all_treatments.family.cons.format.taxonomy'
wetlab_file <- 'data/wetlab_assays.tsv'
metadata_file <- 'data/metadata.tsv'

# Define output files
plot_file <- 'results/figures/figure_1.pdf'

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
shared_otu <- shared_otu[,!(names(shared_otu) %in% c('Otu0004','Otu0308'))] # Remove residual C. difficile OTUs
shared_otu$numOtus <- NULL
shared_otu$label <- NULL
taxonomy_family <- read.delim(taxonomy_family_file, sep='\t', header=T)
shared_family <- read.delim(shared_family_file, sep='\t', header=T, row.names=2)
shared_family <- shared_family[!rownames(shared_family) %in% c('CefC5M2'), ]  # Remove contaminated sample
shared_family <- shared_family[ , !(names(shared_family) == 'Otu008')] # Remove residual C. difficile OTU
shared_family$numOtus <- NULL
shared_family$label <- NULL
wetlab <- read.delim(wetlab_file, sep='\t', header=T, row.names=1)
otu_tax <- read.delim(otu_tax_file, sep='\t', header=T, row.names=1)
rm(shared_family_file, taxonomy_family_file, metadata_file, wetlab_file, otu_tax_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format data

# OTU shared file
shared_otu <- rrarefy(shared_otu, ceiling(min(rowSums(shared_otu)) * 0.9))
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
strep_shared_otu <- filter_table(strep_shared_otu)
clinda_shared_otu <- filter_table(clinda_shared_otu)
cef_shared_otu[,2:ncol(cef_shared_otu)] <- log10(cef_shared_otu[,2:ncol(cef_shared_otu)] + 1)
strep_shared_otu[,2:ncol(strep_shared_otu)] <- log10(strep_shared_otu[,2:ncol(strep_shared_otu)] + 1)
clinda_shared_otu[,2:ncol(clinda_shared_otu)] <- log10(clinda_shared_otu[,2:ncol(clinda_shared_otu)] + 1)

# Filter for significant OTUs by Wilcoxon
# Cef
cef_infected_otu <- subset(cef_shared_otu, infection == '630')
cef_infected_otu$infection <- NULL
cef_mock_otu <- subset(cef_shared_otu, infection == 'mock')
cef_mock_otu$infection <- NULL
rm(cef_shared_otu)
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
cef_pvalues <- round(p.adjust(cef_pvalues, method='BH'), 3)
for (index in 1:length(cef_pvalues)){
  if (cef_pvalues[index] > 0.05){
    cef_infected_otu[,index] <- NULL
    cef_mock_otu[,index] <- NULL
  }
  else {
    if (cef_pvalues[index] <= 0.001) {
      cef_pvalues[index] <- paste('= ', as.character(cef_pvalues[index]), '***', sep='')
    }
    else if (cef_pvalues[index] <= 0.01) {
      cef_pvalues[index] <- paste('= ', as.character(cef_pvalues[index]), '**', sep='')
    }
    else {
      cef_pvalues[index] <- paste('= ', as.character(cef_pvalues[index]), '*', sep='')
    }
  }
}
# Clinda
clinda_infected_otu <- subset(clinda_shared_otu, infection == '630')
clinda_infected_otu$infection <- NULL
clinda_mock_otu <- subset(clinda_shared_otu, infection == 'mock')
clinda_mock_otu$infection <- NULL
rm(clinda_shared_otu)
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
clinda_pvalues <- round(p.adjust(clinda_pvalues, method='BH'), 3)
for (index in 1:length(clinda_pvalues)){
  if (cef_pvalues[index] > 0.05){
    clinda_infected_otu[,index] <- NULL
    clinda_mock_otu[,index] <- NULL
  }
  else {
    if (clinda_pvalues[index] <= 0.001) {
      clinda_pvalues[index] <- paste('= ', as.character(clinda_pvalues[index]), '***', sep='')
    }
    else if (clinda_pvalues[index] <= 0.01) {
      clinda_pvalues[index] <- paste('= ', as.character(clinda_pvalues[index]), '**', sep='')
    }
    else {
      clinda_pvalues[index] <- paste('= ', as.character(clinda_pvalues[index]), '*', sep='')
    }
  }
}
# Strep
strep_infected_otu <- subset(strep_shared_otu, infection == '630')
strep_infected_otu$infection <- NULL
strep_mock_otu <- subset(strep_shared_otu, infection == 'mock')
strep_mock_otu$infection <- NULL
rm(strep_shared_otu)
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
strep_pvalues <- round(p.adjust(strep_pvalues, method='BH'), 3)
for (index in 1:length(strep_pvalues)){
  if (strep_pvalues[index] > 0.05){
    strep_infected_otu[,index] <- NULL
    strep_mock_otu[,index] <- NULL
  }
  else {
    if (strep_pvalues[index] <= 0.001) {
      strep_pvalues[index] <- paste('= ', as.character(clinda_pvalues[index]), '***', sep='')
    }
    else if (strep_pvalues[index] <= 0.01) {
      strep_pvalues[index] <- paste('= ', as.character(strep_pvalues[index]), '**', sep='')
    }
    else {
      strep_pvalues[index] <- paste('= ', as.character(strep_pvalues[index]), '*', sep='')
    }
  }
}

# Rename OTUs with species-level identifier
cef_infected_otu <- clean_merge(otu_tax, t(cef_infected_otu))
cef_genera <- as.vector(cef_infected_otu$genus)
cef_phyla <- as.vector(cef_infected_otu$phylum)
cef_otus <- as.vector(cef_infected_otu$OTU_short)
cef_infected_otu$genus <- NULL
cef_infected_otu$phylum <- NULL
cef_infected_otu$OTU_short <- NULL
cef_infected_otu <- as.data.frame(t(cef_infected_otu))
cef_mock_otu <- clean_merge(otu_tax, t(cef_mock_otu))
cef_mock_otu$genus <- NULL
cef_mock_otu$phylum <- NULL
cef_mock_otu$OTU_short <- NULL
cef_mock_otu <- as.data.frame(t(cef_mock_otu))
clinda_infected_otu <- clean_merge(otu_tax, t(clinda_infected_otu))
clinda_genera <- as.vector(clinda_infected_otu$genus)
clinda_phyla <- as.vector(clinda_infected_otu$phylum)
clinda_otus <- as.vector(clinda_infected_otu$OTU_short)
clinda_infected_otu$genus <- NULL
clinda_infected_otu$phylum <- NULL
clinda_infected_otu$OTU_short <- NULL
clinda_infected_otu <- as.data.frame(t(clinda_infected_otu))
clinda_mock_otu <- clean_merge(otu_tax, t(clinda_mock_otu))
clinda_mock_otu$genus <- NULL
clinda_mock_otu$phylum <- NULL
clinda_mock_otu$OTU_short <- NULL
clinda_mock_otu <- as.data.frame(t(clinda_mock_otu))
strep_infected_otu <- clean_merge(otu_tax, t(strep_infected_otu))
strep_genera <- as.vector(strep_infected_otu$genus)
strep_phyla <- as.vector(strep_infected_otu$phylum)
strep_otus <- as.vector(strep_infected_otu$OTU_short)
strep_infected_otu$genus <- NULL
strep_infected_otu$phylum <- NULL
strep_infected_otu$OTU_short <- NULL
strep_infected_otu <- as.data.frame(t(strep_infected_otu))
strep_mock_otu <- clean_merge(otu_tax, t(strep_mock_otu))
strep_mock_otu$genus <- NULL
strep_mock_otu$phylum <- NULL
strep_mock_otu$OTU_short <- NULL
strep_mock_otu <- as.data.frame(t(strep_mock_otu))
rm(otu_tax)

#--------------------------#

# Phylotype family-level shared file

# Convert to relative abundance
relabund_family <- (shared_family / rowSums(shared_family)) * 100
rm(shared_family)

# Bin lowly abundant OTUs into an 'Other' category
relabund_family[relabund_family < 1] <- 0
relabund_family <- relabund_family[, colSums(relabund_family != 0) > 0]
top_otus <- colnames(relabund_family)

# Subset family-level taxonomy
taxonomy_family <- subset(taxonomy_family, taxonomy_family$otu %in% top_otus)
taxonomy_family$phylum <- gsub('_', ' ', taxonomy_family$phylum)
taxonomy_family$family <- gsub('_', ' ', taxonomy_family$family)
taxonomy_family <- taxonomy_family[order(taxonomy_family$family),]
rm(top_otus)

# Calculate left over abundances
relabund_family$Other <- 100 - rowSums(relabund_family)

# Define group colors
taxonomy_family[] <- lapply(taxonomy_family, as.character)
taxonomy_family <- rbind(taxonomy_family, c('Other', 'Other', 'Other (<1% each)'))
taxonomy_family <- taxonomy_family[order(taxonomy_family$phylum), ]
family_colors <- c('chartreuse3',
                   'navajowhite4',
                   'mediumblue','royalblue3','dodgerblue2','deepskyblue','powderblue',
                   'darkred','firebrick3','red2','brown3','tomato2','coral3','salmon',
                   'white',
                   'darkgoldenrod1','#CCCC00',
                   'darkmagenta')
taxonomy_family$color <- family_colors
bacteria <- taxonomy_family[which(taxonomy_family$family == 'Bacteria unclassified'),]
taxonomy_family <- subset(taxonomy_family, family != 'Bacteria unclassified')
other <- taxonomy_family[which(taxonomy_family$family == 'Other (<1% each)'),]
taxonomy_family <- subset(taxonomy_family, family != 'Other (<1% each)')
taxonomy_family <- rbind(taxonomy_family, bacteria)
taxonomy_family <- rbind(taxonomy_family, other)
rm(bacteria, other)

# Add empty columns for plotting and sort table
relabund_family <- relabund_family[, taxonomy_family$otu] # reorder shared according to taxonomy
relabund_family <- clean_merge(metadata, relabund_family)
relabund_family$type <- NULL
relabund_family$infection <- NULL
relabund_family$susceptibility <- NULL
empty_columns <- as.data.frame(rbind(c('col_1', rep(0, ncol(relabund_family)-1)),
                                     c('col_2', rep(0, ncol(relabund_family)-1)),
                                     c('col_3', rep(0, ncol(relabund_family)-1))))
colnames(empty_columns) <- colnames(relabund_family)
relabund_family <- rbind(relabund_family, empty_columns)
relabund_family <- relabund_family[order(match(relabund_family$abx, c('none', 'col_1', 'streptomycin', 'col_2', 'cefoperazone', 'col_3', 'clindamycin'))),] # sort by treatment
relabund_family$abx <- NULL
rm(metadata, empty_columns)

#--------------------------#

# Format wetlab assay data
wetlab <- subset(wetlab, infection == '630') # Remove uninfected controls
wetlab <- subset(wetlab, treatment != 'germfree') # Remove germfree
wetlab$infection <- NULL
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
rect(xleft=-4, ybottom=3.8, xright=0, ytop=4.2, col='gray70', border='black')
Arrows(x0=-4, y0=4, x1=3.5, y1=4, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,0,2,2.75), y0=c(4.4,4.4,4.4,4.3), x1=c(-4,0,2,2.75), y1=c(3.6,3.6,3.6,3.7), 
         lwd=4, col=c('black','black','black','black'))
segments(x0=c(-4,-3,-2,-1,1), y0=c(4.25,4.25,4.25,4.25,4.25), x1=c(-4,-3,-2,-1,1), y1=c(3.75,3.75,3.75,3.75,3.75), lwd=2)
points(x=c(2,2.75), y=c(4.8,4.8), pch=25, bg=c('white','black'), col='black', cex=2.5)
text(x=c(-4,0,2), y=c(3.1,3.1,3.1), c('Day -7', 'Day -2', 'Day 0'), cex=1.1)
text(x=-4.4, y=4, '1', font=2, cex=1.5)

# IP injection abx timeline
Arrows(x0=-4, y0=1, x1=-1.5, y1=1, lwd=4, arr.type='triangle', arr.length=0.6, arr.width=0.2)
segments(x0=c(-4,-3,-2.25), y0=c(0.6,0.6,0.7), x1=c(-4,-3,-2.25), y1=c(1.4,1.4,1.3), lwd=4, col=c('black','black','black'))
points(x=c(-4,-3,-2.25), y=c(1.8,1.8,1.8), pch=c(25,25,25), bg=c('gray70','white','black'), col='black', cex=2.5)
text(x=c(-4,-3), y=c(0.2,0.2), c('Day -1', 'Day 0'), cex=1.1)
text(-4.4, 1, '2', font=2, cex=1.5)

# Legend
legend(x=-0.6, y=2, legend=expression('Antibiotic in Drinking Water', 'Antibiotic IP Injection',paste(italic('C. difficile'), ' Spore Gavage'), 'Necropsy (18 hpi)'), 
       pt.bg=c('gray70','gray70','white','black'), pch=c(22,25,25,25), cex=1.2, pt.cex=c(3,2,2,2), bty='n')

# Route of administration
text(x=c(-2, 2), y=c(-0.6,-0.6), c('In Drinking Water:', 'IP Injected:'), cex=1.2, font=2)
text(x=-2, y=c(-1,-1.4), c('Streptomycin','Cefoperazone'), cex=1.2, col=c(cef_col, strep_col), font=2)
text(x=2, y=-1, 'Clindamycin', cex=1.2, col=clinda_col, font=2)

# Plot label
text(-4.7, 4.88, 'a', cex=1.5, font=2)

#-------------------------------------------------------------------#

# CFU and toxin data

# Vegetative cells
par(mar=c(3,4,1,4), mgp=c(2.5, 0.75, 0))
stripchart(cfu_vegetative~treatment, data=wetlab, col='black', bg='firebrick2', xlim=c(0,22), ylim=c(0,9), pch=21,
           vertical=TRUE, at=c(0.5, 6.5, 12.5, 18.5), xaxt='n', yaxt='n', ylab='CFU/g Cecal Content', 
           cex=1.7, method='jitter', jitter=0.2)
abline(h=2, lwd=1.5, col='gray30', lty=5) # LOD
abline(v=c(5,11,17), lwd=1.5) # dividers
axis(side=2, at=seq(0,9,1), labels=c(0, parse(text=paste(rep(10,9), '^', seq(1,9,1), sep=''))), las=1)

mtext(c('No Antibiotics','Streptomycin-pretreated','Cefoperazone-pretreated','Clindamycin-pretreated'), col=c('black',strep_col, cef_col, clinda_col), 
      at=c(2,8,14,20), side=1, font=2, cex=0.75, padj=0.75)

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
par(new=TRUE, xpd=TRUE)
stripchart(toxin_titer~treatment, data=wetlab, col='black', bg='green3', xlim=c(0,22), ylim=c(1.6,3.4), pch=23,
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

legend('topleft', legend=c('Vegetative cells (CFU)','Spores (CFU)'), ncol=1, bty='n',
       pch=21, col='black', pt.bg=c('firebrick2','blue2'), pt.cex=1.9)
legend('topright', legend='Toxin titer', bty='n',
       pch=23, col='black', pt.bg='green3', pt.cex=1.9)
box()

# Add significance
text(x=c(6.5,12.5,18.5,8,14,20,15.5,21.5), y=c(3.4,3.4,3.4,3,3,3,3.2,2.5), labels=rep('*',8), col=noabx_col, font=2, cex=2.2)

mtext('b', side=2, line=2, las=2, adj=3, padj=-11, cex=1.0, font=2)

#-------------------------------------------------------------------#

# Family-level phylotype bar chart
par(mar=c(5,5,1,1), mgp=c(2.5, 0.25, 0), new=FALSE, xpd=FALSE)
barplot(t(rev(relabund_family)), col=rev(taxonomy_family$color), yaxt='n', xaxt='n', 
        ylim=c(0,100), ylab='Relative Abundance', cex.names=1.2, space=0)
box()
axis(side=2, at=seq(0,100,20), labels=c('0%','20%','40%','60%','80%','100%'), tick=FALSE, las=1)
abline(h=c(20,40,60,80), lty=2)

# Lines under plot
arrows(0.3, -2, 7.7, -2, angle=0, length=0, lwd=2, xpd=TRUE) # no abx
arrows(9.3, -2, 17.7, -2, angle=0, length=0, lwd=2, xpd=TRUE) # strep - cdi
arrows(18.3, -2, 26.7, -2, angle=0, length=0, lwd=2, xpd=TRUE) # strep - mock
arrows(28.3, -2, 35.7, -2, angle=0, length=0, lwd=2, xpd=TRUE) # cef - cdi
arrows(36.3, -2, 44.7, -2, angle=0, length=0, lwd=2, xpd=TRUE) # cef - mock
arrows(46.3, -2, 54.7, -2, angle=0, length=0, lwd=2, xpd=TRUE) # clinda - cdi
arrows(55.3, -2, 63.7, -2, angle=0, length=0, lwd=2, xpd=TRUE) # clinda - mock

mtext(rep('CDI',3), side=1, at=c(13.5,32,50.5), adj=0.5, padj=1, cex=0.7)
mtext(rep('Mock',4), side=1, at=c(4,22.5,40.5,59.5), adj=0.5, padj=1, cex=0.7)
mtext(c('No Antibiotics','Streptomycin-pretreated','Cefoperazone-pretreated','Clindamycin-pretreated'), 
      side=1, at=c(4,18,36,55), adj=0.5, padj=2.5, cex=0.75, font=2, col=c('black',strep_col, cef_col, clinda_col))
mtext('c', side=2, line=2, las=2, adj=3, padj=-11, cex=1.0, font=2)

# Create a figure legend in empty plot
par(mar=c(4,0,0.3,5))
plot(0, type='n', ylim=c(-5,5), xlim=c(5,5), ylab='', xlab='', xaxt='n', yaxt='n', axes=FALSE)
legend('right', legend=taxonomy_family$family, pt.bg=taxonomy_family$color, 
       pch=22, pt.cex=2.5, cex=1.1, bty='n')

# Add in phylum classifications
segments(x0=c(4.8,4.8,4.8,4.8,4.8), x1=c(4.8,4.8,4.8,4.8,4.8), 
         y0=c(5,4.4,1.55,-2.4,-3.5), 
         y1=c(4.6,1.85,-2.1,-3.25,-3.9), lwd=3) # vertical
segments(x0=c(4.4,4.4,4.28,4.43,4.48), x1=c(4.8,4.8,4.8,4.8,4.8), 
         y0=c(4.8,3.125,-0.275,-2.8,-3.7), 
         y1=c(4.8,3.125,-0.275,-2.8,-3.7), lwd=2) # horizontal
text(x=c(3.75,3.75,3.75,3.75,3.75), y=c(4.8,3.125,-0.275,-2.8,-3.7), cex=1.2,
     labels=c('Actinobacteria', 'Bacteroidetes', 'Firmicutes', 'Proteobacteria', 'Verrucomicrobia'))

#-------------------------------------------------------------------#

# Significant OTUs and relative abundance changes

# Streptomycin plot
par(mar=c(4, 14, 2, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE, yaxs='i')
plot(1, type='n', ylim=c(0,length(clinda_otus)*2), xlim=c(0,4), 
     ylab='', xlab=expression(paste('Normalized Abundance (',log[10],')')), xaxt='n', yaxt='n')
title('Streptomycin-pretreated', line=0.5, cex.main=1.1, font.main=2, col.main=strep_col)
index <- 1
for(i in colnames(strep_mock_otu)){
  stripchart(at=index-0.35, jitter(strep_mock_otu[,i], amount=1e-5), 
             pch=21, bg="mediumseagreen", method="jitter", jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  stripchart(at=index+0.35, jitter(strep_infected_otu[,i], amount=1e-5), 
             pch=21, bg="mediumorchid3", method="jitter", jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  if (i != colnames(strep_mock_otu)[length(colnames(strep_mock_otu))]){
    abline(h=index+1, lty=2)
  }
  segments(median(strep_mock_otu[,i]), index-0.4, median(strep_mock_otu[,i]), index, lwd=2) #adds line for median
  segments(median(strep_infected_otu[,i]), index+0.4, median(strep_infected_otu[,i]), index, lwd=2)
  index <- index + 2
}
axis(1, at=seq(0,4,1), label=c('0','10','100','1000','10000'))
minor.ticks.axis(1, 10, mn=0, mx=4)
legend('topright', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))), 'Mock-infected'),
       pch=c(21, 21), pt.bg=c("mediumorchid3","mediumseagreen"), bg='white', pt.cex=1.4, cex=0.9)
formatted <- lapply(1:length(strep_otus), function(i) bquote(paste(italic(.(strep_genera[i])), .(strep_otus[i]), sep=' ')))
axis(2, at=seq(1,index-2,2)+0.4, labels=do.call(expression, formatted), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 
axis(2, at=seq(1,index-2,2), labels=strep_phyla, las=1, line=-0.5, tick=F, cex.axis=1.1) 
formatted <- lapply(1:length(strep_pvalues), function(i) bquote(paste(italic('p'), .(strep_pvalues[i]), sep=' ')))
axis(2, at=seq(1,index-2,2)-0.5, labels=do.call(expression, formatted), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 

mtext('d', side=2, line=2, las=2, adj=14, padj=-11, cex=1.0, font=2)

#-----------------#

# Cefoperazone plot
plot(1, type='n', ylim=c(0,length(clinda_otus)*2), xlim=c(0,4), 
     ylab='', xlab=expression(paste('Normalized Abundance (',log[10],')')), xaxt='n', yaxt='n')
title('Cefoperazone-pretreated', line=0.5, cex.main=1.1, font.main=2, col.main=cef_col)
index <- 1
for(i in colnames(cef_mock_otu)){
  stripchart(at=index-0.35, jitter(cef_mock_otu[,i], amount=1e-5), 
             pch=21, bg="mediumseagreen", method="jitter", jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  stripchart(at=index+0.35, jitter(cef_infected_otu[,i], amount=1e-5), 
             pch=21, bg="mediumorchid3", method="jitter", jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  if (i != colnames(cef_mock_otu)[length(colnames(cef_mock_otu))]){
  abline(h=index+1, lty=2)
  }
  segments(median(cef_mock_otu[,i]), index-0.4, median(cef_mock_otu[,i]), index, lwd=2) #adds line for median
  segments(median(cef_infected_otu[,i]), index+0.4, median(cef_infected_otu[,i]), index, lwd=2)
  index <- index + 2
}
axis(1, at=seq(0,4,1), label=c('0','10','100','1000','10000'))
minor.ticks.axis(1, 10, mn=0, mx=4)
legend('topright', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))), 'Mock-infected'),
       pch=c(21, 21), pt.bg=c("mediumorchid3","mediumseagreen"), bg='white', pt.cex=1.4, cex=0.9)
formatted <- lapply(1:length(cef_otus), function(i) bquote(paste(italic(.(cef_genera[i])), .(cef_otus[i]), sep=' ')))
axis(2, at=seq(1,index-2,2)+0.4, labels=do.call(expression, formatted), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 
axis(2, at=seq(1,index-2,2), labels=cef_phyla, las=1, line=-0.5, tick=F, cex.axis=1.1) 
formatted <- lapply(1:length(cef_pvalues), function(i) bquote(paste(italic('p'), .(cef_pvalues[i]), sep=' ')))
axis(2, at=seq(1,index-2,2)-0.5, labels=do.call(expression, formatted), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 

mtext('e', side=2, line=2, las=2, adj=16, padj=-11, cex=1.0, font=2)
  
#-----------------#

# Clindamycin plot
plot(1, type='n', ylim=c(0,length(clinda_otus)*2), xlim=c(0,4), 
     ylab='', xlab=expression(paste('Normalized Abundance (',log[10],')')), xaxt='n', yaxt='n')
title('Clindamycin-pretreated', line=0.5, cex.main=1.1, font.main=2, col.main=clinda_col)
index <- 1
for(i in colnames(clinda_mock_otu)){
  stripchart(at=index-0.35, jitter(clinda_mock_otu[,i], amount=1e-5), 
             pch=21, bg="mediumseagreen", method="jitter", jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  stripchart(at=index+0.35, jitter(clinda_infected_otu[,i], amount=1e-5), 
             pch=21, bg="mediumorchid3", method="jitter", jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  if (i != colnames(clinda_mock_otu)[length(colnames(clinda_mock_otu))]){
    abline(h=index+1, lty=2)
  }
  segments(median(clinda_mock_otu[,i]), index-0.4, median(clinda_mock_otu[,i]), index, lwd=2) #adds line for median
  segments(median(clinda_infected_otu[,i]), index+0.4, median(clinda_infected_otu[,i]), index, lwd=2)
  index <- index + 2
}
axis(1, at=seq(0,4,1), label=c('0','10','100','1000','10000'))
minor.ticks.axis(1, 10, mn=0, mx=4)
legend('topright', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))), 'Mock-infected'),
       pch=c(21, 21), pt.bg=c("mediumorchid3","mediumseagreen"), bg='white', pt.cex=1.4, cex=0.9)
formatted <- lapply(1:length(clinda_otus), function(i) bquote(paste(italic(.(clinda_genera[i])), .(clinda_otus[i]), sep=' ')))
axis(2, at=seq(1,index-2,2)+0.4, labels=do.call(expression, formatted), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 
axis(2, at=seq(1,index-2,2), labels=clinda_phyla, las=1, line=-0.5, tick=F, cex.axis=1.1) 
formatted <- lapply(1:length(clinda_pvalues), function(i) bquote(paste(italic('p'), .(clinda_pvalues[i]), sep=' ')))
axis(2, at=seq(1,index-2,2)-0.5, labels=do.call(expression, formatted), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 

mtext('f', side=2, line=2, las=2, adj=27, padj=-11, cex=1.0, font=2)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()
setwd(starting_dir)

