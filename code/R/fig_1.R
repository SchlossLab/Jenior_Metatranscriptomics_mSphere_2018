
# Start with clear environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define input files
wetlab_file <- 'data/wetlab_assays.tsv'
metadata_file <- 'data/metadata.tsv'
cfu <- 'data/cfu_time.tsv'

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
wetlab <- read.delim(wetlab_file, sep='\t', header=T, row.names=1)
cfu <- read.delim(cfu, sep='\t', header=TRUE)
rm(metadata_file, wetlab_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

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

# Format CFU oveer time
# Aggregate median by abx
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

#----------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
pdf(file=plot_file, width=8, height=7)
layout(matrix(c(1,2,
                3,3),
              nrow=2, ncol=2, byrow = TRUE))
minor_ticks <- c(0.18,0.34,0.48,0.6,0.7,0.78,0.84,0.9,0.94,0.98)
short_ticks <- c(log10(strep_size/1000)-0.33,log10(strep_size/1000)-0.2,log10(strep_size/1000)-0.1,log10(strep_size/1000)-0.05)

#----------------------------------------------------------------------------------------------------------------------#

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
       pt.bg=c('gray70','gray70','white','black'), pch=c(22,25,25,25), cex=1.1, pt.cex=c(3,2,2,2), bty='n')

# Route of administration
text(x=c(-2.5, 2.3), y=c(-0.6,-0.6), c('In Drinking Water:', 'IP Injected:'), cex=1.2, font=2)
text(x=-2.5, y=c(-1,-1.4), c('Streptomycin (5.0 mg/ml)','Cefoperazone (0.5 mg/ml)'), cex=1.2, col=c(cef_col, strep_col))
text(x=2.3, y=-1, 'Clindamycin (10 mg/kg)', cex=1.2, col=clinda_col)

# Plot label
text(-4.7, 4.88, 'A', cex=1.5)

#-------------------------------------------------------------------#

# CFU over time
par(mar=c(4,4,1,1), las=1, mgp=c(2.5, 0.75, 0))
plot(0, type='n', xlab='Days Post-Infection', ylab='Total CFU/g Cecal Content', xaxt='n', yaxt='n', xlim=c(0,11), ylim=c(0,10))
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
legend(x=6.5, y=8, legend=c('Streptomycin (5 mg/ml)', 'Cefoperazone (0.5 mg/ml)', 'Clindamycin (10 mg/kg)', 'No Antibiotics'),
       pch=16, col=c(strep_col, cef_col, clinda_col, noabx_col), cex=0.7, pt.cex=1.2)
mtext('B', side=2, line=2, las=2, adj=2.3, padj=-9, cex=1.2)

#-------------------------------------------------------------------#

# CFU and toxin data

# Vegetative cells
par(mar=c(3,4,1,4), mgp=c(2.5, 0.75, 0))
stripchart(cfu_vegetative~treatment, data=wetlab, col='black', bg='chocolate2', xlim=c(0,22), ylim=c(0,9), pch=21,
           vertical=TRUE, at=c(0.5, 6.5, 12.5, 18.5), xaxt='n', yaxt='n', ylab='CFU/g Cecal Content', cex.lab=1.2,
           cex=1.7, method='jitter', jitter=0.2)
abline(h=2, lwd=1.5, col='gray30', lty=5) # LOD
abline(v=c(5,11,17), lwd=1.5) # dividers
axis(side=2, at=seq(0,9,1), labels=c(0, parse(text=paste(rep(10,9), '^', seq(1,9,1), sep=''))), las=1)

mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin'), col=c('black',strep_col, cef_col, clinda_col), 
      at=c(2,8,14,20), side=1, cex=1, padj=0.75, font=2)

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
stripchart(cfu_spore~treatment, data=wetlab, col='black', bg='chartreuse3', xlim=c(0,22), ylim=c(0,9), pch=21,
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
par(new=TRUE, xpd=TRUE, las=0)
stripchart(toxin_titer~treatment, data=wetlab, col='black', bg='deeppink3', xlim=c(0,22), ylim=c(1.6,3.4), pch=23,
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
mtext(expression(paste('Toxin Titer/g Cecal Content (',log[10],')')), side=4, line=3, padj=-0.5, cex=0.9)
legend('topleft', legend=c('Vegetative cells (CFU)','Spores (CFU)'), ncol=1, bty='n',
       pch=21, col='black', pt.bg=c('chocolate2','chartreuse3'), pt.cex=1.9)
legend('topright', legend='Toxin titer', bty='n',
       pch=23, col='black', pt.bg='deeppink3', pt.cex=1.9)
box()

# Add significance
text(x=c(6.5,12.5,18.5,8,14,20,15.5,21.5), y=c(3.4,3.4,3.4,3,3,3,3.2,2.5), labels=rep('*',8), col=noabx_col, font=2, cex=2.2)

mtext('C', side=2, line=2, las=2, adj=2.3, padj=-9.5, cex=1.2)

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


