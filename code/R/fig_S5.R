
# Start with a blank slate
rm(list=ls())
gc()

# Load dependencies
deps <- c('wesanderson', 'plotrix');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Function for population variance of columns in a matrix
pop_var <- function(data) {
  vars <- c()
  for (x in 1:ncol(data)){
    vars[x] <- sum((data[,x] - mean(data[,x]))^2) / length(data[,x])
  }
  return(vars)
}

# Function for sample variance of columns in a matrix
samp_var <- function(data) {
  vars <- c()
  for (x in 1:ncol(data)){
    vars[x] <- var(data[,x])
  }
  return(vars)
}

# Define files
metadata <- '~/Desktop/Repositories/Jenior_Modeling_mSystems_2017/data/metadata.tsv'
metabolome <- '~/Desktop/Repositories/Jenior_Modeling_mSystems_2017/data/wetlab_assays/metabolomics.scaled_intensities.tsv'
shared <- '~/Desktop/Repositories/Jenior_Modeling_mSystems_2017/data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'
ko_var <- '~/Desktop/Repositories/Jenior_Modeling_mSystems_2017/data/mapping/variance_ko.tsv'
cfu <- '~/Desktop/Repositories/Jenior_Modeling_mSystems_2017/data/wetlab_assays/cfu.dat'

#----------------------------------------#

# Read data
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2','StrepC4M1'), ] # Remove possible contamination
metabolome <- read.delim(metabolome, sep='\t', header=T, row.names=1)
metabolome <- metabolome[, !colnames(metabolome) %in% c('CefC5M2','StrepC4M1')] # Remove possible contamination
shared <- read.delim(shared, sep='\t', header=T, row.names=2)
shared <- shared[!rownames(shared) %in% c('CefC5M2','StrepC4M1'), ] # Remove possible contamination
ko_var <- read.delim(ko_var, sep='\t', header=T, row.names=1)
cfu <- read.delim(cfu, sep='\t', header=T)

# Subset and format KO variation data
clpP <- as.numeric(ko_var['K01358',c(1:3)])
thrS <- as.numeric(ko_var['K01868',c(1:3)])
gyrA <- as.numeric(ko_var['K02469',c(1:3)])
kar_var <- as.data.frame(cbind(gyrA,thrS,clpP))
kar_var$abx <- c('Streptomycin','Cefoperazone','Clindamycin')
rm(gyrA,thrS,clpP,ko_var)

# Format data
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$type <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome <- as.data.frame(t(metabolome))
metabolites <- colnames(metabolome)
shared$label <- NULL
shared$numOtus <- NULL
shared <- log10(shared + 10)
otus <- colnames(shared)
cfu$mouse <- NULL
cfu$cfu_spore <- NULL
cfu <- subset(cfu, cage < 4 ) # Remove uninfected controls
cfu$cage <- NULL
cfu <- subset(cfu, cfu$treatment != 'conventional')
cfu <- subset(cfu, cfu$treatment != 'germfree')
cfu$treatment <- factor(cfu$treatment, levels=c('streptomycin', 'cefoperazone', 'clindamycin'))
cfu$cfu_vegetative <- log10(cfu$cfu_vegetative)

# Merge datasets
shared <- merge(metadata, shared, by='row.names')
rownames(shared) <- shared$Row.names
shared$Row.names <- NULL
metabolome <- merge(metadata, metabolome, by='row.names')
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
rm(metadata)

# Subset datasets
strep <- subset(metabolome, abx == 'streptomycin')
strep$abx <- NULL
strep_metabolome_mock <- subset(strep, infection == 'mock')
strep_metabolome_mock$infection <- NULL
strep_metabolome_630 <- subset(strep, infection == '630')
strep_metabolome_630$infection <- NULL
cef <- subset(metabolome, abx == 'cefoperazone')
cef$abx <- NULL
cef_metabolome_mock <- subset(cef, infection == 'mock')
cef_metabolome_mock$infection <- NULL
cef_metabolome_630 <- subset(cef, infection == '630')
cef_metabolome_630$infection <- NULL
clinda <- subset(metabolome, abx == 'clindamycin')
clinda$abx <- NULL
clinda_metabolome_mock <- subset(clinda, infection == 'mock')
clinda_metabolome_mock$infection <- NULL
clinda_metabolome_630 <- subset(clinda, infection == '630')
clinda_metabolome_630$infection <- NULL
conv <- subset(metabolome, abx == 'none')
conv$abx <- NULL
conv_metabolome_mock <- subset(conv, infection == 'mock')
conv_metabolome_mock$infection <- NULL
rm(metabolome)
strep <- subset(shared, abx == 'streptomycin')
strep$abx <- NULL
strep_shared_mock <- subset(strep, infection == 'mock')
strep_shared_mock$infection <- NULL
strep_shared_630 <- subset(strep, infection == '630')
strep_shared_630$infection <- NULL
cef <- subset(shared, abx == 'cefoperazone')
cef$abx <- NULL
cef_shared_mock <- subset(cef, infection == 'mock')
cef_shared_mock$infection <- NULL
cef_shared_630 <- subset(cef, infection == '630')
cef_shared_630$infection <- NULL
clinda <- subset(shared, abx == 'clindamycin')
clinda$abx <- NULL
clinda_shared_mock <- subset(clinda, infection == 'mock')
clinda_shared_mock$infection <- NULL
clinda_shared_630 <- subset(clinda, infection == '630')
clinda_shared_630$infection <- NULL
conv <- subset(shared, abx == 'none')
conv$abx <- NULL
conv_shared_mock <- subset(conv, infection == 'mock')
conv_shared_mock$infection <- NULL
rm(shared)
rm(strep, cef, clinda, conv)
cef_cfu <- as.numeric(cfu[cfu$treatment == 'cefoperazone', 2])
strep_cfu <- as.numeric(cfu[cfu$treatment == 'streptomycin', 2])
clinda_cfu <- as.numeric(cfu[cfu$treatment == 'clindamycin', 2])
rm(cfu)

#----------------------------------------#

# Calculate sample variance
strep_metabolome_mock <- samp_var(strep_metabolome_mock)
strep_metabolome_630 <- samp_var(strep_metabolome_630)
cef_metabolome_mock <- samp_var(cef_metabolome_mock)
cef_metabolome_630 <- samp_var(cef_metabolome_630)
clinda_metabolome_mock <- samp_var(clinda_metabolome_mock)
clinda_metabolome_630 <- samp_var(clinda_metabolome_630)
conv_metabolome_mock <- samp_var(conv_metabolome_mock)
strep_shared_mock <- samp_var(strep_shared_mock)
strep_shared_630 <- samp_var(strep_shared_630)
cef_shared_mock <- samp_var(cef_shared_mock)
cef_shared_630 <- samp_var(cef_shared_630)
clinda_shared_mock <- samp_var(clinda_shared_mock)
clinda_shared_630 <- samp_var(clinda_shared_630)
conv_shared_mock <- samp_var(conv_shared_mock)
cfu_var <- c(var(strep_cfu),var(cef_cfu),var(clinda_cfu)) 
rm(strep_cfu,cef_cfu,clinda_cfu) 

# Conserved colors across studies and figures
strep_col <- '#D37A1F'
cef_col <- '#3A9CBC'
clinda_col <- '#A40019'
noabx_col <- 'gray40'

#----------------------------------------#

# Generate plot
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/results/supplement/figures/figure_S5.pdf'
pdf(file=plot_file, width=9, height=7)
layout(matrix(c(1,
                2), nrow=2, ncol=1, byrow=TRUE))

# Housekeeping genes
#par(mar=c(3,5,1,1), las=1, mgp=c(3,0.7,0))
#plot(0, type='n', xlab='', xaxt='n', ylab='Normalized cDNA Abundance', xlim=c(0,14), ylim=c(0,100), cex.lab=1.4, cex.axis=1.4)
#legend('topright', legend=c('Streptomycin-pretreated','Cefoperazone-pretreated','Clindamycin-pretreated'),
#       pt.bg=c(strep_col, cef_col, clinda_col), pch=21, cex=1.6, pt.cex=2.6, col='black', bty='n')
#stripchart(at=2, kar_var[,1], vertical=T, pch=21, bg=c(strep_col,cef_col,clinda_col), 
#           method='jitter', jitter=0.5, cex=2.5, lwd=1, add=TRUE)
#stripchart(at=7, kar_var[,2], vertical=T, pch=21, bg=c(strep_col,cef_col,clinda_col), 
#           method='jitter', jitter=0.5, cex=2.5, lwd=1, add=TRUE)
#stripchart(at=12, kar_var[,3], vertical=T, pch=21, bg=c(strep_col,cef_col,clinda_col), 
#           method='jitter', jitter=0.5, cex=2.5, lwd=1, add=TRUE)
#segments(x0=c(1,6,11), y0=c(median(kar_var[,1]),median(kar_var[,2]),median(kar_var[,3])), 
#         x1=c(3,8,13), y1=c(median(kar_var[,1]),median(kar_var[,2]),median(kar_var[,3])), lwd=4)
#legend('topleft', legend='Housekeeping Genes', pt.cex=0, bty='n', cex=1.6)
#text(cex=1.8, x=c(2.8,7.6,12.3), y=-12, c('GyrA','ThrS','ClpP'), xpd=TRUE, pos=2)
#mtext('A', side=2, line=2, las=2, adj=2, padj=-7, cex=1.5)
#box(lwd=2)

# Vegetative C. difficile CFU
#par(las=1, mar=c(3,5,1,1), mgp=c(3,0.7,0))
#stripchart(at=1, cfu_var, vertical=T, pch=21, bg=c(strep_col,cef_col,clinda_col), 
#           method='jitter', jitter=0.75, cex=2.5, lwd=1, ylim=c(0,1), 
#           ylab='Sample Variance', cex.lab=1.5, cex.axis=1.4)
#segments(x0=0, y0=median(cfu_var), x1=2, y1=median(cfu_var), lwd=4)
#mtext('B', side=2, line=2, las=2, adj=2, padj=-7, cex=1.5)
#legend('topleft', legend='Vegetative CFU', pt.cex=0, bty='n', cex=1.6)
#box(lwd=2)

# 16S
par(las=1, mar=c(3,5,1,1), mgp=c(3,0.7,0))
plot(0, type='n', xlab='', xaxt='n', ylab='Sample Variance', 
     xlim=c(0,11), ylim=c(0,1.5), cex.lab=1.5, cex.axis=1.4)
abline(v=c(3,6,9), lty=2)
box(lwd=2)
stripchart(at=1, strep_shared_mock, vertical=T, pch=21, bg=strep_col, 
           method='jitter', jitter=0.3, cex=1.5, lwd=1, add=TRUE)
stripchart(at=2, strep_shared_630, vertical=T, pch=21, bg=strep_col, 
           method='jitter', jitter=0.3, cex=1.5, lwd=1, add=TRUE)
stripchart(at=4, cef_shared_mock, vertical=T, pch=21, bg=cef_col, 
           method='jitter', jitter=0.3, cex=1.5, lwd=1, add=TRUE)
stripchart(at=5, cef_shared_630, vertical=T, pch=21, bg=cef_col, 
           method='jitter', jitter=0.3, cex=1.5, lwd=1, add=TRUE)
stripchart(at=7, clinda_shared_mock, vertical=T, pch=21, bg=clinda_col, 
           method='jitter', jitter=0.3, cex=1.5, lwd=1, add=TRUE)
stripchart(at=8, clinda_shared_630, vertical=T, pch=21, bg=clinda_col, 
           method='jitter', jitter=0.3, cex=1.5, lwd=1, add=TRUE)
stripchart(at=10, conv_shared_mock, vertical=T, pch=21, bg=noabx_col, 
           method='jitter', jitter=0.3, cex=1.5, lwd=1, add=TRUE)
segments(x0=c(1,2,4,5,7,8,10)-0.45, y0=c(median(strep_shared_mock),median(strep_shared_630),median(cef_shared_mock),median(cef_shared_630),median(clinda_shared_mock),median(clinda_shared_630),median(conv_shared_mock)), 
         x1=c(1,2,4,5,7,8,10)+0.45, y1=c(median(strep_shared_mock),median(strep_shared_630),median(cef_shared_mock),median(cef_shared_630),median(clinda_shared_mock),median(clinda_shared_630),median(conv_shared_mock)), 
         lwd=4)
segments(x0=c(1,2,4,5,7,8,10)-0.3, y0=c(quantile(strep_shared_mock)[4],quantile(strep_shared_630)[4],quantile(cef_shared_mock)[4],quantile(cef_shared_630)[4],quantile(clinda_shared_mock)[4],quantile(clinda_shared_630)[4],quantile(conv_shared_mock)[4]), 
         x1=c(1,2,4,5,7,8,10)+0.3, y1=c(quantile(strep_shared_mock)[4],quantile(strep_shared_630)[4],quantile(cef_shared_mock)[4],quantile(cef_shared_630)[4],quantile(clinda_shared_mock)[4],quantile(clinda_shared_630)[4],quantile(conv_shared_mock)[4]), 
         lwd=4, col='gray')
mtext('CDI:', side=1, at=0, padj=0.8, cex=1)
mtext(c('+','-','+','-','+','-','-'), side=1, 
      at=c(1,2,4,5,7,8,10), padj=0.5, cex=1.5)
mtext(c('Streptomycin','Cefoperazone','Clindamycin','No Antibiotics'), side=1, 
      at=c(1.5,4.5,7.5,10), padj=2, cex=1.2)
mtext('A', side=2, line=2, las=2, adj=2.5, padj=-7, cex=1.5)
legend('topleft', legend='OTU Abundances', pt.cex=0, bty='n', cex=1.2)

# Metabolome
par(las=1, mar=c(4,5,1,1), mgp=c(3,0.7,0))
plot(0, type='n', xlab='', xaxt='n', ylab='Sample Variance', 
     xlim=c(0,11), ylim=c(0,500), cex.lab=1.5, cex.axis=1.4)
abline(v=c(3,6,9), lty=2)
box(lwd=2)
stripchart(at=1, strep_metabolome_mock, vertical=T, pch=21, bg=strep_col, 
           method='jitter', jitter=0.3, cex=1.5, lwd=1, add=TRUE)
stripchart(at=2, strep_metabolome_630, vertical=T, pch=21, bg=strep_col, 
           method='jitter', jitter=0.3, cex=1.5, lwd=1, add=TRUE)
stripchart(at=4, cef_metabolome_mock, vertical=T, pch=21, bg=cef_col, 
           method='jitter', jitter=0.3, cex=1.5, lwd=1, add=TRUE)
stripchart(at=5, cef_metabolome_630, vertical=T, pch=21, bg=cef_col, 
           method='jitter', jitter=0.3, cex=1.5, lwd=1, add=TRUE)
stripchart(at=7, clinda_metabolome_mock, vertical=T, pch=21, bg=clinda_col, 
           method='jitter', jitter=0.3, cex=1.5, lwd=1, add=TRUE)
stripchart(at=8, clinda_metabolome_630, vertical=T, pch=21, bg=clinda_col, 
           method='jitter', jitter=0.3, cex=1.5, lwd=1, add=TRUE)
stripchart(at=10, conv_metabolome_mock, vertical=T, pch=21, bg=noabx_col, 
           method='jitter', jitter=0.3, cex=1.5, lwd=1, add=TRUE)
segments(x0=c(1,2,4,5,7,8,10)-0.45, y0=c(median(strep_metabolome_mock),median(strep_metabolome_630),median(cef_metabolome_mock),median(cef_metabolome_630),median(clinda_metabolome_mock),median(clinda_metabolome_630),median(conv_metabolome_mock)), 
         x1=c(1,2,4,5,7,8,10)+0.45, y1=c(median(strep_metabolome_mock),median(strep_metabolome_630),median(cef_metabolome_mock),median(cef_metabolome_630),median(clinda_metabolome_mock),median(clinda_metabolome_630),median(conv_metabolome_mock)), 
         lwd=4)
segments(x0=c(1,2,4,5,7,8,10)-0.3, y0=c(quantile(strep_metabolome_mock)[4],quantile(strep_metabolome_630)[4],quantile(cef_metabolome_mock)[4],quantile(cef_metabolome_630)[4],quantile(clinda_metabolome_mock)[4],quantile(clinda_metabolome_630)[4],quantile(conv_metabolome_mock)[4]), 
         x1=c(1,2,4,5,7,8,10)+0.3, y1=c(quantile(strep_metabolome_mock)[4],quantile(strep_metabolome_630)[4],quantile(cef_metabolome_mock)[4],quantile(cef_metabolome_630)[4],quantile(clinda_metabolome_mock)[4],quantile(clinda_metabolome_630)[4],quantile(conv_metabolome_mock)[4]), 
         lwd=4, col='gray')
mtext('CDI:', side=1, at=0, padj=0.8, cex=1)
mtext(c('+','-','+','-','+','-','-'), side=1, 
      at=c(1,2,4,5,7,8,10), padj=0.5, cex=1.5)
mtext(c('Streptomycin','Cefoperazone','Clindamycin','No Antibiotics'), side=1, 
      at=c(1.5,4.5,7.5,10), padj=2, cex=1.2)
mtext('B', side=2, line=2, las=2, adj=2.5, padj=-7, cex=1.5)
legend('topleft', legend='Metabolites', pt.cex=0, bty='n', cex=1.2)


dev.off()

#----------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()
