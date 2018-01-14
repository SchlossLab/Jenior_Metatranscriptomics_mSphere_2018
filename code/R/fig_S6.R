
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
metadata <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.tsv'
metabolome <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/metabolomics.scaled_intensities.tsv'
shared <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'
ko_var <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/mapping/variance_ko.tsv'
cfu <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/wetlab_assays/cfu.dat'

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
kar_var <- c(gyrA,0,thrS,0,clpP)
rm(gyrA,thrS,clpP)

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

# Calculate summary stats for barplots
metabolome <- rbind(quantile(strep_metabolome_mock)[2:4],
                    quantile(strep_metabolome_630)[2:4],
                    quantile(cef_metabolome_mock)[2:4],
                    quantile(cef_metabolome_630)[2:4],
                    quantile(clinda_metabolome_mock)[2:4],
                    quantile(clinda_metabolome_630)[2:4],
                    quantile(conv_metabolome_mock)[2:4])
metabolome <- as.data.frame(metabolome)
colnames(metabolome) <- c('q25','median','q75')
rownames(metabolome) <- c('strep_mock','strep_630','cef_mock','cef_630',
                          'clinda_mock','clinda_630','conv_mock')
metabolome$q75[nrow(metabolome)] <- 0.86
shared <- rbind(quantile(strep_shared_mock)[2:4],
                quantile(strep_shared_630)[2:4],
                quantile(cef_shared_mock)[2:4],
                quantile(cef_shared_630)[2:4],
                quantile(clinda_shared_mock)[2:4],
                quantile(clinda_shared_630)[2:4],
                quantile(conv_shared_mock)[2:4])
shared <- as.data.frame(shared)
colnames(shared) <- c('q25','median','q75')
rownames(shared) <- c('strep_mock','strep_630','cef_mock','cef_630',
                          'clinda_mock','clinda_630','conv_mock')
shared$q75[nrow(shared)] <- 0.00095
cfu_var <- c(var(strep_cfu),var(cef_cfu),var(clinda_cfu))
rm(strep_cfu,cef_cfu,clinda_cfu)

rm(strep_metabolome_mock,strep_metabolome_630,cef_metabolome_mock,cef_metabolome_630,
   clinda_metabolome_mock,clinda_metabolome_630,conv_metabolome_mock,strep_shared_mock,
   strep_shared_630,cef_shared_mock,cef_shared_630,clinda_shared_mock,clinda_shared_630,conv_shared_mock)

#----------------------------------------#

# Generate plot
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/supplement/figures/figure_S6.pdf'
pdf(file=plot_file, width=12, height=10)
layout(matrix(c(1,2,
                3,3,
                4,4), nrow=3, ncol=2, byrow=TRUE))

# Conserved colors across studies and figures
strep_col <- '#D37A1F'
cef_col <- '#3A9CBC'
clinda_col <- '#A40019'
noabx_col <- 'gray40'
gf_col <- 'forestgreen'

# Housekeeping genes
par(mar=c(3,5,1,1), las=1, mgp=c(3,0.7,0))
plot(0, type='n', xlab='', xaxt='n', ylab='Normalized cDNA Abundance', xlim=c(0,17), ylim=c(0,100), yaxs='i')
legend('topleft', legend=c('Streptomycin-pretreated','Cefoperazone-pretreated','Clindamycin-pretreated'),
       pt.bg=c(strep_col, cef_col, clinda_col), pch=22, cex=1.6, pt.cex=2.6, col='black', bty='n')
# Add groups
barplot(kar_var, col=c(strep_col,cef_col,clinda_col,'white',
                       strep_col,cef_col,clinda_col,'white',
                       strep_col,cef_col,clinda_col,'white'), yaxt='n', add=TRUE, yaxs='i')
text(cex=1.6, x=c(3.5,9.5,15.5), y=-9, c('GyrA','ThrS','ClpP'), xpd=TRUE, pos=2)
mtext('A', side=2, line=2, las=2, adj=2, padj=-7, cex=1.5)

# Vegetative C. difficile CFU
par(las=1, mar=c(3,5,1,1), mgp=c(3,0.7,0), yaxs='i')
barplot(cfu_var, ylim=c(0,1), ylab='Sample Variance',
        col=c(strep_col,cef_col,clinda_col))
box()
mtext(c('Streptomycin','Cefoperazone','Clindamycin'), side=1, 
      at=c(0.7,1.9,3.1,4.3), padj=2, cex=1.2)
mtext('B', side=2, line=2, las=2, adj=2, padj=-7, cex=1.5)
legend('topleft', legend='Vegetative CFU (Log10)', pt.cex=0, bty='n', cex=1.6)

# 16S
par(las=1, mar=c(3,5,1,1), mgp=c(3,0.7,0), yaxs='i')
barplot(shared$median, xaxt='n', yaxt='n', ylim=c(0,0.001), ylab='Sample Variance',
        col=c(strep_col,strep_col,cef_col,cef_col,clinda_col,clinda_col,noabx_col))
segments(x0=c(0.7,1.9,3.1,4.3,5.5,6.7,7.9), y0=shared$q25, x1=c(0.7,1.9,3.1,4.3,5.5,6.7,7.9), y1=shared$q75)
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.9)
mtext(c('+','-','+','-','+','-','-'), side=1, 
      at=c(0.7,1.9,3.1,4.3,5.5,6.7,7.9), padj=0.5, cex=1.5)
mtext(c('Streptomycin','Cefoperazone','Clindamycin','No Antibiotics'), side=1, 
      at=c(1.3,3.7,6.1,7.9), padj=2, cex=0.9)
abline(v=c(2.5,4.9,7.3), lty=2)
axis(side=2, at=c(0,0.0002,0.0004,0.0006,0.001), labels=c('0.0','0.0002','0.0004','0.0006','0.004'), cex.axis=0.8)
axis.break(2, 0.0008, style='slash') 
rect(xleft=7.7, xright=8.1, ytop=0.00081, ybottom=0.00079, col='white', border='white')
segments(x0=c(-1,-1,8.72),y0=c(0,0.001,0),x1=c(10,10,8.72),y1=c(0,0.001,0.001), lwd=2)
mtext('C', side=2, line=2, las=2, adj=2, padj=-7, cex=1.5)
legend('topleft', legend='OTU Abundance', pt.cex=0, bty='n', cex=1.2)
box(lwd=2)

# Metabolome
par(las=1, mar=c(3,5,1,1), mgp=c(3,0.7,0), yaxs='i')
barplot(metabolome$median, xaxt='n', yaxt='n', ylim=c(0,0.9), ylab='Sample Variance',
        col=c(strep_col,strep_col,cef_col,cef_col,clinda_col,clinda_col,noabx_col))
segments(x0=c(0.7,1.9,3.1,4.3,5.5,6.7,7.9), y0=metabolome$q25, x1=c(0.7,1.9,3.1,4.3,5.5,6.7,7.9), y1=metabolome$q75, lwd=2)
mtext('CDI:', side=1, at=0, padj=0.5, cex=0.9)
mtext(c('+','-','+','-','+','-','-'), side=1, 
      at=c(0.7,1.9,3.1,4.3,5.5,6.7,7.9), padj=0.5, cex=1.5)
mtext(c('Streptomycin','Cefoperazone','Clindamycin','No Antibiotics'), side=1, 
      at=c(1.3,3.7,6.1,7.9), padj=2, cex=0.9)
abline(v=c(2.5,4.9,7.3), lty=2)
axis(side=2, at=c(0,0.2,0.4,0.6,0.9), labels=c('0.0','0.2','0.4','0.6','9.0'), cex=1.5)
axis.break(2, 0.8, style='slash') 
segments(x0=c(-1,-1,8.72),y0=c(0,0.9,0),x1=c(10,10,8.72),y1=c(0,0.9,0.9), lwd=2)
rect(xleft=7.7, xright=8.1, ytop=0.81, ybottom=0.79, col='white', border='white')
mtext('D', side=2, line=2, las=2, adj=2, padj=-7, cex=1.5)
legend('topleft', legend='Metabolome (Log10)', pt.cex=0, bty='n', cex=1.2)
bow(lwd=2)

dev.off()

#----------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()
