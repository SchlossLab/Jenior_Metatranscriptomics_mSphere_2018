
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
metadata_file <- 'data/metadata.tsv'

# LEfSe results
cef_lefse <- 'data/16S_analysis/lefse/cef.0.03.filter.0.03.lefse_summary'
clinda_lefse <- 'data/16S_analysis/lefse/clinda.0.03.filter.0.03.lefse_summary'
strep_lefse <- 'data/16S_analysis/lefse/strep.0.03.filter.0.03.lefse_summary'

# Define output files
plot_file <- 'results/figures/figure_2.pdf'

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
shared_family <- shared_family[ ,!(names(shared_family) == 'Otu008')] # Remove residual C. difficile OTU
shared_family$numOtus <- NULL
shared_family$label <- NULL
otu_tax <- read.delim(otu_tax_file, sep='\t', header=T, row.names=1)
cef_lefse <- read.delim(cef_lefse, sep='\t', header=T, row.names=1)
clinda_lefse <- read.delim(clinda_lefse, sep='\t', header=T, row.names=1)
strep_lefse <- read.delim(strep_lefse, sep='\t', header=T, row.names=1)
rm(shared_family_file, taxonomy_family_file, metadata_file, otu_tax_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format data

# OTU shared files
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

# OTU names
otu_tax$genus <- gsub('_', ' ', otu_tax$genus)
otu_tax$genus <- gsub('Ruminococcus2', 'Ruminococcus', otu_tax$genus)

# Normalize reads
cef_size <- ceiling(min(rowSums(cef_shared_otu[,2:ncol(cef_shared_otu)])) * 0.9)
strep_size <- ceiling(min(rowSums(strep_shared_otu[,2:ncol(strep_shared_otu)])) * 0.9)
clinda_size <- ceiling(min(rowSums(clinda_shared_otu[,2:ncol(clinda_shared_otu)])) * 0.9)
cef_shared_otu[,2:ncol(cef_shared_otu)] <- rrarefy(cef_shared_otu[,2:ncol(cef_shared_otu)], cef_size)
cef_shared_otu[,2:ncol(cef_shared_otu)] <- filter_table(cef_shared_otu[,2:ncol(cef_shared_otu)])
strep_shared_otu[,2:ncol(strep_shared_otu)] <- rrarefy(strep_shared_otu[,2:ncol(strep_shared_otu)], strep_size)
strep_shared_otu[,2:ncol(strep_shared_otu)] <- filter_table(strep_shared_otu[,2:ncol(strep_shared_otu)])
clinda_shared_otu[,2:ncol(clinda_shared_otu)] <- rrarefy(clinda_shared_otu[,2:ncol(clinda_shared_otu)], clinda_size)
clinda_shared_otu[,2:ncol(clinda_shared_otu)] <- filter_table(clinda_shared_otu[,2:ncol(clinda_shared_otu)])

# Retrieve significant lefse results
cef_lefse <- clean_merge(cef_lefse, otu_tax)
cef_lefse <- cef_lefse[complete.cases(cef_lefse),]
cef_lefse$LogMaxMean <- NULL
cef_lefse$pValue <- round(cef_lefse$pValue, 3)
cef_lefse <- cef_lefse[!(rownames(cef_lefse) %in% c('Otu0004','Otu0308','Otu0110')),]
cef_lefse <- cef_lefse[order(-cef_lefse$LDA),]
clinda_lefse <- clean_merge(clinda_lefse, otu_tax)
clinda_lefse <- clinda_lefse[complete.cases(clinda_lefse),]
clinda_lefse$LogMaxMean <- NULL
clinda_lefse$pValue <- round(clinda_lefse$pValue, 3)
clinda_lefse <- clinda_lefse[!(rownames(clinda_lefse) %in% c('Otu0004','Otu0308','Otu0110')),]
clinda_lefse <- clinda_lefse[order(-clinda_lefse$LDA),]
strep_lefse <- clean_merge(strep_lefse, otu_tax)
strep_lefse <- strep_lefse[complete.cases(strep_lefse),]
strep_lefse$LogMaxMean <- NULL
strep_lefse$pValue <- round(strep_lefse$pValue, 3)
strep_lefse <- strep_lefse[!(rownames(strep_lefse) %in% c('Otu0004','Otu0308','Otu0110')),]
strep_lefse <- strep_lefse[order(-strep_lefse$LDA),]

# Subset OTU abundances
cef_shared_otu <- cef_shared_otu[,c('infection', rownames(cef_lefse))]
cef_infected_otu <- subset(cef_shared_otu, infection == '630')
cef_infected_otu$infection <- NULL
cef_infected_otu <- cef_infected_otu[, rownames(cef_lefse)]
cef_mock_otu <- subset(cef_shared_otu, infection == 'mock')
cef_mock_otu$infection <- NULL
cef_mock_otu <- cef_mock_otu[, rownames(cef_lefse)]
rm(cef_shared_otu)
clinda_shared_otu <- clinda_shared_otu[,c('infection',rownames(clinda_lefse))]
clinda_infected_otu <- subset(clinda_shared_otu, infection == '630')
clinda_infected_otu$infection <- NULL
clinda_infected_otu <- clinda_infected_otu[, rownames(clinda_lefse)]
clinda_mock_otu <- subset(clinda_shared_otu, infection == 'mock')
clinda_mock_otu$infection <- NULL
clinda_mock_otu <- clinda_mock_otu[, rownames(clinda_lefse)]
rm(clinda_shared_otu)
strep_shared_otu <- strep_shared_otu[,c('infection',rownames(strep_lefse))]
strep_infected_otu <- subset(strep_shared_otu, infection == '630')
strep_infected_otu$infection <- NULL
strep_infected_otu <- strep_infected_otu[, rownames(strep_lefse)]
strep_mock_otu <- subset(strep_shared_otu, infection == 'mock')
strep_mock_otu <- strep_mock_otu[, rownames(strep_lefse)]
strep_mock_otu$infection <- NULL
rm(strep_shared_otu)

# Transform abundances
cef_infected_otu <- log10(cef_infected_otu + 1)
cef_mock_otu <- log10(cef_mock_otu + 1)
clinda_infected_otu <- log10(clinda_infected_otu + 1)
clinda_mock_otu <- log10(clinda_mock_otu + 1)
strep_infected_otu <- log10(strep_infected_otu + 1)
strep_mock_otu <- log10(strep_mock_otu + 1)

#--------------------------#

# Family-level shared file

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

#----------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
pdf(file=plot_file, width=11, height=7)
layout(matrix(c(1,1,2,
                3,4,5),
              nrow=2, ncol=3, byrow = TRUE))
minor_ticks <- c(0.18,0.34,0.48,0.6,0.7,0.78,0.84,0.9,0.94,0.98)
short_ticks <- c(log10(strep_size/1000)-0.33,log10(strep_size/1000)-0.2,log10(strep_size/1000)-0.1,log10(strep_size/1000)-0.05)

#-------------------------------------------------------------------#

# Family-level phylotype bar chart
par(mar=c(5,5,1,1), mgp=c(2.5, 0.25, 0), new=FALSE, xpd=FALSE)
barplot(t(rev(relabund_family)), col=rev(taxonomy_family$color), yaxt='n', xaxt='n', cex.lab=1.3,
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
      side=1, at=c(4,18,36,55), adj=0.5, padj=2.5, cex=0.85, col=c('black',strep_col, cef_col, clinda_col))
mtext('A', side=2, line=2, las=2, adj=2.5, padj=-8, cex=1.3, font=2)

# Create a figure legend in empty plot
par(mar=c(4,0,0.3,3))
plot(0, type='n', ylim=c(-5,5), xlim=c(5,5), ylab='', xlab='', xaxt='n', yaxt='n', axes=FALSE)
legend('topright', legend=taxonomy_family$family, pt.bg=taxonomy_family$color, 
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
par(mar=c(4, 13, 2, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE, yaxs='i')
plot(1, type='n', ylim=c(0,nrow(strep_lefse)*2), xlim=c(0,log10(strep_size)), 
     ylab='', xlab='Relative Abundance %', xaxt='n', yaxt='n')
title('Streptomycin-pretreated', line=0.5, cex.main=1.3, col.main=strep_col, font.main=1)
index <- 1
for(i in c(1:ncol(strep_mock_otu))){
  stripchart(at=index-0.35, jitter(strep_mock_otu[,i], amount=1e-5), 
             pch=21, bg='chartreuse2', method='jitter', jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  stripchart(at=index+0.35, jitter(strep_infected_otu[,i], amount=1e-5), 
             pch=21, bg='mediumorchid4', method='jitter', jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  if (i != ncol(strep_mock_otu)){
    abline(h=index+1, lty=2)
  }
  segments(median(strep_mock_otu[,i]), index-0.6, median(strep_mock_otu[,i]), index, lwd=2) #adds line for median
  segments(median(strep_infected_otu[,i]), index+0.6, median(strep_infected_otu[,i]), index, lwd=2)
  #if (wilcox.test(strep_mock_otu[,i], strep_infected_otu[,i], exact=FALSE)$p.value <= 0.05){
  #  text(x=log10(strep_size)+0.25, y=index, labels='*', font=2, cex=1.8, xpd=TRUE)
  #}
  index <- index + 2
}
axis(1, at=c(0,(log10(strep_size/1000)),(log10(strep_size/100)),(log10(strep_size/10)),log10(strep_size)), labels=c('0','0.1','1','10','100')) 
legend('topright', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))), 'Mock-infected'),
       pch=c(21, 21), pt.bg=c('mediumorchid4','chartreuse2'), bg='white', pt.cex=1.4, cex=0.9)
formatted_names <- lapply(1:nrow(strep_lefse), function(i) bquote(paste(italic(.(strep_lefse$genus[i])), ' ', .(as.vector(strep_lefse$OTU)[i]), sep='')))
axis(2, at=seq(1,index-2,2)+0.4, labels=do.call(expression, formatted_names), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 
formatted_p <- lapply(1:nrow(strep_lefse), function(i) bquote(paste(.(as.vector(strep_lefse$phylum)[i]), '; ', italic('p'), ' = ', .(strep_lefse$pValue[i]), sep='')))
axis(2, at=seq(1,index-2,2)-0.4, labels=do.call(expression, formatted_p), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 
axis(side=1, at=short_ticks+0.03, label=rep('',length(short_ticks)), tck=-0.01)
axis(side=1, at=minor_ticks+0.42, label=rep('',length(minor_ticks)), tck=-0.01)
axis(side=1, at=minor_ticks+1.42, label=rep('',length(minor_ticks)), tck=-0.01)
axis(side=1, at=minor_ticks+2.42, label=rep('',length(minor_ticks)), tck=-0.01)

mtext('B', side=2, line=2, las=2, adj=9, padj=-9, cex=1.3, font=2)

#-----------------#

# Cefoperazone plot
plot(1, type='n', ylim=c(0,nrow(cef_lefse)*2), xlim=c(0,log10(cef_size)), 
     ylab='', xlab='Relative Abundance %', xaxt='n', yaxt='n')
title('Cefoperazone-pretreated', line=0.5, cex.main=1.3, col.main=cef_col, font.main=1)
index <- 1
for(i in c(1:ncol(cef_mock_otu))){
  stripchart(at=index-0.35, jitter(cef_mock_otu[,i], amount=1e-5), 
             pch=21, bg='chartreuse2', method='jitter', jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  stripchart(at=index+0.35, jitter(cef_infected_otu[,i], amount=1e-5), 
             pch=21, bg='mediumorchid4', method='jitter', jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  if (i != ncol(cef_mock_otu)){
    abline(h=index+1, lty=2)
  }
  segments(median(cef_mock_otu[,i]), index-0.6, median(cef_mock_otu[,i]), index, lwd=2) #adds line for median
  segments(median(cef_infected_otu[,i]), index+0.6, median(cef_infected_otu[,i]), index, lwd=2)
  #if (wilcox.test(cef_mock_otu[,i], cef_infected_otu[,i], exact=FALSE)$p.value <= 0.05){
  #  text(x=log10(cef_size)+0.25, y=index, labels='*', font=2, cex=1.8, xpd=TRUE)
  #}
  index <- index + 2
}
axis(1, at=c(0,(log10(cef_size/1000)),(log10(cef_size/100)),(log10(cef_size/10)),log10(cef_size)), labels=c('0','0.1','1','10','100')) 
legend('topright', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))), 'Mock-infected'),
       pch=c(21, 21), pt.bg=c('mediumorchid4','chartreuse2'), bg='white', pt.cex=1.4, cex=0.9)
formatted_names <- lapply(1:nrow(cef_lefse), function(i) bquote(paste(italic(.(cef_lefse$genus[i])), ' ', .(as.vector(cef_lefse$OTU)[i]), sep='')))
axis(2, at=seq(1,index-2,2)+0.4, labels=do.call(expression, formatted_names), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 
formatted_p <- lapply(1:nrow(cef_lefse), function(i) bquote(paste(.(as.vector(cef_lefse$phylum)[i]), '; ', italic('p'), ' = ', .(cef_lefse$pValue[i]), sep='')))
axis(2, at=seq(1,index-2,2)-0.4, labels=do.call(expression, formatted_p), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 
cef_ticks <- c(log10(cef_size/1000)-0.55,log10(cef_size/1000)-0.35,log10(cef_size/1000)-0.2,log10(cef_size/1000)-0.1,log10(cef_size/1000)-0.05)
axis(side=1, at=cef_ticks+0.02, label=rep('',length(cef_ticks)), tck=-0.01)
axis(side=1, at=minor_ticks+0.74, label=rep('',length(minor_ticks)), tck=-0.01)
axis(side=1, at=minor_ticks+1.74, label=rep('',length(minor_ticks)), tck=-0.01)
axis(side=1, at=minor_ticks+2.74, label=rep('',length(minor_ticks)), tck=-0.01)

mtext('C', side=2, line=2, las=2, adj=9, padj=-9, cex=1.3, font=2)

#-----------------#

# Clindamycin plot
plot(1, type='n', ylim=c(0,nrow(clinda_lefse)*2), xlim=c(0,log10(clinda_size)), 
     ylab='', xlab='Relative Abundance %', xaxt='n', yaxt='n')
title('Clindamycin-pretreated', line=0.5, cex.main=1.3, col.main=clinda_col, font.main=1)
index <- 1
for(i in c(1:ncol(clinda_mock_otu))){
  stripchart(at=index-0.25, jitter(clinda_mock_otu[,i], amount=1e-5), 
             pch=21, bg='chartreuse2', method='jitter', jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  stripchart(at=index+0.25, jitter(clinda_infected_otu[,i], amount=1e-5), 
             pch=21, bg='mediumorchid4', method='jitter', jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  if (i != ncol(clinda_mock_otu)){
    abline(h=index+1, lty=2)
  }
  segments(median(clinda_mock_otu[,i]), index-0.6, median(clinda_mock_otu[,i]), index, lwd=2) #adds line for median
  segments(median(clinda_infected_otu[,i]), index+0.6, median(clinda_infected_otu[,i]), index, lwd=2)
  #if (wilcox.test(clinda_mock_otu[,i], clinda_infected_otu[,i], exact=FALSE)$p.value <= 0.05){
  #  text(x=log10(clinda_size)+0.25, y=index, labels='*', font=2, cex=1.8, xpd=TRUE)
  #}
  index <- index + 2
}
axis(1, at=c(0,(log10(clinda_size/1000)),(log10(clinda_size/100)),(log10(clinda_size/10)),log10(clinda_size)), labels=c('0','0.1','1','10','100')) 
legend('topright', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))), 'Mock-infected'),
       pch=c(21, 21), pt.bg=c('mediumorchid4','chartreuse2'), bg='white', pt.cex=1.4, cex=0.9)
formatted_names <- lapply(1:nrow(clinda_lefse), function(i) bquote(paste(italic(.(clinda_lefse$genus[i])), ' ', .(as.vector(clinda_lefse$OTU)[i]), sep='')))
axis(2, at=seq(1,index-2,2)+0.25, labels=do.call(expression, formatted_names), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 
formatted_p <- lapply(1:nrow(clinda_lefse), function(i) bquote(paste(.(as.vector(clinda_lefse$phylum)[i]), '; ', italic('p'), ' = ', .(clinda_lefse$pValue[i]), sep='')))
axis(2, at=seq(1,index-2,2)-0.25, labels=do.call(expression, formatted_p), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 
axis(side=1, at=short_ticks, label=rep('',length(short_ticks)), tck=-0.01)
axis(side=1, at=minor_ticks+0.4, label=rep('',length(minor_ticks)), tck=-0.01)
axis(side=1, at=minor_ticks+1.4, label=rep('',length(minor_ticks)), tck=-0.01)
axis(side=1, at=minor_ticks+2.4, label=rep('',length(minor_ticks)), tck=-0.01)

mtext('D', side=2, line=2, las=2, adj=9, padj=-9, cex=1.3, font=2)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()


