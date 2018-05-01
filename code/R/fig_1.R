
# Start with clear environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')

# Define input files
wetlab_file <- 'data/wetlab_assays.tsv'
metadata_file <- 'data/metadata.tsv'
cfu <- 'data/cfu_time.tsv'
shared_family_file <- 'data/16S_analysis/all_treatments.family.subsample.shared'
taxonomy_family_file <- 'data/16S_analysis/all_treatments.family.cons.format.taxonomy'

# Define output files
plot_file <- 'results/figures/figure_1.pdf'

#-------------------------------------------------------------------------------------------------------------------------------------#

# Load in data
metadata <- read.delim(metadata_file, sep='\t', header=T, row.names=1)
metadata <- subset(metadata, type == 'conventional') # remove germfree
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
wetlab <- read.delim(wetlab_file, sep='\t', header=T, row.names=1)
cfu <- read.delim(cfu, sep='\t', header=TRUE)
taxonomy_family <- read.delim(taxonomy_family_file, sep='\t', header=T)
shared_family <- read.delim(shared_family_file, sep='\t', header=T, row.names=2)
shared_family <- shared_family[!rownames(shared_family) %in% c('CefC5M2'), ]  # Remove contaminated sample
shared_family <- shared_family[ ,!(names(shared_family) == 'Otu008')] # Remove residual C. difficile OTU
shared_family$numOtus <- NULL
shared_family$label <- NULL
rm(metadata_file, wetlab_file, taxonomy_family_file, shared_family_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format wetlab assay data
wetlab <- subset(wetlab, infection == '630') # Remove uninfected controls
wetlab <- subset(wetlab, treatment != 'germfree') # Remove germfree
wetlab$infection <- NULL
wetlab$cfu_spore <- NULL
wetlab$toxin_titer <- NULL
wetlab$cfu_vegetative <- as.numeric(wetlab$cfu_vegetative)
wetlab$cfu_vegetative[wetlab$cfu_vegetative == 0] <- 100
wetlab$cfu_vegetative <- log10(wetlab$cfu_vegetative)
wetlab$treatment <- droplevels(wetlab$treatment)
wetlab$treatment <- factor(wetlab$treatment, levels=c('conventional', 'streptomycin', 'cefoperazone', 'clindamycin'))
wetlab$cfu_vegetative[wetlab$cfu_vegetative <= 2.0] <- 1.7 # Undetectable points below LOD
noabx_veg <- as.numeric(subset(wetlab, treatment == 'conventional')$cfu_vegetative)
strep_veg <- as.numeric(subset(wetlab, treatment == 'streptomycin')$cfu_vegetative)
cef_veg <- as.numeric(subset(wetlab, treatment == 'cefoperazone')$cfu_vegetative)
clinda_veg <- as.numeric(subset(wetlab, treatment == 'clindamycin')$cfu_vegetative)
rm(wetlab)

# Test differences in time series
cfu_dtw <- read.delim('data/cfu_time_test.tsv', sep='\t', header=TRUE)
cfu_dtw_dist <- dist(cfu_dtw[,2:ncol(cfu_dtw)], method='DTW')
cfu_dtw_pval <- adonis(cfu_dtw_dist ~ cfu_dtw$abx, cfu_dtw[,2:ncol(cfu_dtw)], perm=999)$aov.tab[[6]][1]
cfu_dtw_pval <- as.character(round(cfu_dtw_pval, 3))
rm(cfu_dtw_dist, cfu_dtw)
# Format CFU over time
cfu[,2:ncol(cfu)] <- log10(cfu[,2:ncol(cfu)] + 1)
cfu_median <- aggregate(cfu[,2:ncol(cfu)], by=list(cfu$abx), FUN=quantile, probs=0.5)
rownames(cfu_median) <- cfu_median$Group.1
cfu_median$Group.1 <- NULL
cfu_median <- as.data.frame(t(cfu_median))
cfu_q25 <- aggregate(cfu[,2:ncol(cfu)], by=list(cfu$abx), FUN=quantile, probs=0.25)
rownames(cfu_q25) <- cfu_q25$Group.1
cfu_q25$Group.1 <- NULL
cfu_q25 <- as.data.frame(t(cfu_q25))
cfu_q75 <- aggregate(cfu[,2:ncol(cfu)], by=list(cfu$abx), FUN=quantile, probs=0.75)
rownames(cfu_q75) <- cfu_q75$Group.1
cfu_q75$Group.1 <- NULL
cfu_q75 <- as.data.frame(t(cfu_q75))
rm(cfu)

#--------------------------#

# Family-level shared file

# Convert to relative abundance
relabund_family <- (shared_family / rowSums(shared_family)) * 100
rm(shared_family)

# Bin lowly abundant OTUs into an 'Other' category
relabund_family[relabund_family < 1] <- 0
relabund_family <- relabund_family[, which(colSums(relabund_family) > 0)]
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
pdf(file=plot_file, width=12, height=4)
layout(matrix(c(1,2,2,3,4),
              nrow=1, ncol=5, byrow = TRUE))
minor_ticks <- c(0.18,0.34,0.48,0.6,0.7,0.78,0.84,0.9,0.94,0.98)
short_ticks <- c(log10(strep_size/1000)-0.33,log10(strep_size/1000)-0.2,log10(strep_size/1000)-0.1,log10(strep_size/1000)-0.05)

#----------------------------------------------------------------------------------------------------------------------#

# CFU over time
par(mar=c(4,4,1,1), las=1, mgp=c(2.5, 0.75, 0))
plot(0, type='n', xlab='Days Post-Infection', ylab='Total cfu/g Feces', xaxt='n', yaxt='n', xlim=c(0,11), ylim=c(0,10))
lines(cfu_median$streptomycin, lwd=2.5, col=strep_col, type='b', pch=19) 
segments(x0=c(1:11), y0=cfu_q25$streptomycin, x1=c(1:11), y1=cfu_q75$streptomycin, col=strep_col, lwd=2.5)
lines(cfu_median$cefoperazone, lwd=2.5, col=cef_col, type='b', pch=19) 
segments(x0=c(1:11), y0=cfu_q25$cefoperazone, x1=c(1:11), y1=cfu_q75$cefoperazone, col=cef_col, lwd=2.5)
lines(cfu_median$clindamycin, lwd=2.5, col=clinda_col, type='b', pch=19) 
segments(x0=c(1:11), y0=cfu_q25$clindamycin, x1=c(1:11), y1=cfu_q75$clindamycin, col=clinda_col, lwd=2.5)
lines(cfu_median$none, lwd=2.5, col=noabx_col, type='b', pch=19) 
segments(x0=c(1:11), y0=cfu_q25$none, x1=c(1:11), y1=cfu_q75$none, col=noabx_col, lwd=2.5)
abline(h=2, lwd=2, lty=2)
axis(side=1, at=c(0:11), labels=c(-1:10))
axis(side=2, at=seq(0,10,1), labels=c(0, parse(text=paste(rep(10,10), '^', seq(1,10,1), sep=''))), las=1)
legend(x=1.7, y=4, legend=c('Streptomycin', 'Cefoperazone', 'Clindamycin', 'No Antibiotics'),
       pch=16, col=c(strep_col, cef_col, clinda_col, noabx_col), cex=0.9, pt.cex=1.5)
text(x=9, y=4, '*', cex=2, font=2)
legend('topright', legend=as.expression(bquote(paste(italic('p'),' = ', .(cfu_dtw_pval)))), 
       pch=1, cex=1.2, pt.cex=0, bty='n')
mtext('A', side=2, line=2, las=2, adj=1.4, padj=-9.2, cex=1.4, font=2)
box(lwd=1.5)

#-----------------------------#

# Family-level phylotype bar chart
par(mar=c(3,4,1,1), mgp=c(2.5, 0.25, 0), new=FALSE, xpd=FALSE)
barplot(t(rev(relabund_family)), col=rev(taxonomy_family$color), yaxt='n', xaxt='n', cex.lab=1.3,
        ylim=c(0,100), ylab='Relative Abundance', cex.names=1.2, space=0)
box(lwd=1.5)
axis(side=2, at=seq(0,100,20), labels=c('0%','20%','40%','60%','80%','100%'), tick=FALSE, las=1)
abline(h=c(20,40,60,80), lty=2)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin'), font=2,
      side=1, at=c(4,18,36,55), adj=0.5, padj=1, cex=0.7, col='black')
mtext('B', side=2, line=2, las=2, adj=1.4, padj=-9.8, cex=1.4, font=2)

# Create a figure legend in empty plot
par(mar=c(0,0,0,1))
plot(0, type='n', ylim=c(-10,10), xlim=c(5,5), ylab='', xlab='', xaxt='n', yaxt='n', axes=FALSE)
legend('right', legend=taxonomy_family$family, pt.bg=taxonomy_family$color, 
       pch=22, pt.cex=2, cex=0.9, bty='n')
# Add in phylum classifications
segments(x0=rep(4.45,5), x1=rep(4.45,5), 
         y0=c(5.65,4.9,1.7, -2.75,-4), y1=c(5.25,2.2,-2.25, -3.55,-4.4), 
         lwd=3) # vertical
text(x=rep(3.7,5), y=c(5.45,3.5,-0.3,-3.15,-4.15), cex=0.9,
     labels=c('Actinobacteria', 'Bacteroidetes', 'Firmicutes', 'Proteobacteria', 'Verrucomicrobia'))
text(x=c(3.75,5.5), y=-7, labels=c('Phylum','Family'), cex=1.1, font=2, xpd=TRUE)

#-----------------------------#

# Vegetative cells
par(mar=c(5,4,1,1), mgp=c(2.2, 0.75, 0))
stripchart(noabx_veg, bg=noabx_col, xlim=c(0.25,2.25), ylim=c(0,9), pch=21,
           vertical=TRUE, at=0.5, xaxt='n', yaxt='n', ylab='Vegetative cfu/g Cecal Content', lwd=2,
           cex=1.7, method='jitter', jitter=0.12)
stripchart(strep_veg, bg=strep_col, xlim=c(0.25,2.25), ylim=c(0,9), pch=21,
           vertical=TRUE, at=1, xaxt='n', yaxt='n', ylab='Vegetative cfu/g Cecal Content', lwd=2,
           cex=1.7, method='jitter', jitter=0.12, add=TRUE)
stripchart(cef_veg, bg=cef_col, xlim=c(0.25,2.25), ylim=c(0,9), pch=21,
           vertical=TRUE, at=1.5, xaxt='n', yaxt='n', ylab='Vegetative cfu/g Cecal Content', lwd=2,
           cex=1.7, method='jitter', jitter=0.12, add=TRUE)
stripchart(clinda_veg, bg=clinda_col, xlim=c(0.25,2.25), ylim=c(0,9), pch=21,
           vertical=TRUE, at=2, xaxt='n', yaxt='n', ylab='Vegetative cfu/g Cecal Content', lwd=2,
           cex=1.7, method='jitter', jitter=0.12, add=TRUE)
abline(h=2, lwd=1.5, col='gray30', lty=5) # LOD
axis(side=2, at=seq(0,9,1), labels=c(0, parse(text=paste(rep(10,9), '^', seq(1,9,1), sep=''))), las=1)
text(x=c(1,1.5,2), y=rep(9,3), labels=rep('*',3), font=2, cex=2)
legend('bottomleft', legend='18 hours post-infection', bty='n', pt.cex=0)
box(lwd=1.5)
text(c(0.1,0.6,1.1,1.6), -2.1, adj = 0, srt=45, xpd = TRUE, font=2,
     labels=c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin'),
     col='black')

# Median lines
segments(x0=c(0.3, 0.8, 1.3, 1.8), 
         y0=c(as.numeric(median(noabx_veg)),
              as.numeric(median(strep_veg)),
              as.numeric(median(cef_veg)),
              as.numeric(median(clinda_veg))),
         x1=c(0.7, 1.2, 1.7, 2.2), 
         y1=c(as.numeric(median(noabx_veg)),
              as.numeric(median(strep_veg)),
              as.numeric(median(cef_veg)),
              as.numeric(median(clinda_veg))), 
         lwd=3)

mtext('C', side=2, line=2, las=2, adj=1, padj=-9, cex=1.4, font=2)


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


