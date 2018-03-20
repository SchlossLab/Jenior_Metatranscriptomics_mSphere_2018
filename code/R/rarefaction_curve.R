
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')

# Define files
# MetaGs
cef_630_metaG <- 'data/read_mapping/metagenome/cefoperazone_630.Cefoperazone.metaG.final.pool.norm.txt'
cef_mock_metaG <- 'data/read_mapping/metagenome/cefoperazone_mock.Cefoperazone.metaG.final.pool.norm.txt'
clinda_630_metaG <- 'data/read_mapping/metagenome/clindamycin_630.Clindamycin.metaG.final.pool.norm.txt'
clinda_mock_metaG <- 'data/read_mapping/metagenome/clindamycin_mock.Clindamycin.metaG.final.pool.norm.txt'
strep_630_metaG <- 'data/read_mapping/metagenome/streptomycin_630.Streptomycin.metaG.final.pool.norm.txt'
strep_mock_metaG <- 'data/read_mapping/metagenome/streptomycin_mock.Streptomycin.metaG.final.pool.norm.txt'
noabx_mock_metaG <- 'data/read_mapping/metagenome/conventional_mock.Conventional.metaG.final.pool.norm.txt'
# MetaTs
cef_630_metaT <- 'data/read_mapping/metatranscriptome/cefoperazone_630.Cefoperazone.metaT.final.pool.norm.txt'
cef_mock_metaT <- 'data/read_mapping/metatranscriptome/cefoperazone_mock.Cefoperazone.metaT.final.pool.norm.txt'
clinda_630_metaT <- 'data/read_mapping/metatranscriptome/clindamycin_630.Clindamycin.metaT.final.pool.norm.txt'
clinda_mock_metaT <- 'data/read_mapping/metatranscriptome/clindamycin_mock.Clindamycin.metaT.final.pool.norm.txt'
strep_630_metaT <- 'data/read_mapping/metatranscriptome/streptomycin_630.Streptomycin.metaT.final.pool.norm.txt'
strep_mock_metaT <- 'data/read_mapping/metatranscriptome/streptomycin_mock.Streptomycin.metaT.final.pool.norm.txt'
noabx_mock_metaT <- 'data/read_mapping/metatranscriptome/conventional.Conventional.metaT.final.pool.norm.txt'

#--------------------------------------------------------------------------------------------------#

# Read in data
# MetaGs
cef_630_metaG <- read.delim(cef_630_metaG, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
cef_630_metaG <- cef_630_metaG[cef_630_metaG != 0]
cef_mock_metaG <- read.delim(cef_mock_metaG, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
cef_mock_metaG <- cef_mock_metaG[cef_mock_metaG != 0]
clinda_630_metaG <- read.delim(clinda_630_metaG, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
clinda_630_metaG <- clinda_630_metaG[clinda_630_metaG != 0]
clinda_mock_metaG <- read.delim(clinda_mock_metaG, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
clinda_mock_metaG <- clinda_mock_metaG[clinda_mock_metaG != 0]
strep_630_metaG <- read.delim(strep_630_metaG, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
strep_630_metaG <- strep_630_metaG[strep_630_metaG != 0]
strep_mock_metaG <- read.delim(strep_mock_metaG, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
strep_mock_metaG <- strep_mock_metaG[strep_mock_metaG != 0]
noabx_mock_metaG <- read.delim(noabx_mock_metaG, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
noabx_mock_metaG <- noabx_mock_metaG[noabx_mock_metaG != 0]

# MetaTs
cef_630_metaT <- read.delim(cef_630_metaT, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
cef_630_metaT <- cef_630_metaT[cef_630_metaT != 0]
cef_mock_metaT <- read.delim(cef_mock_metaT, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
cef_mock_metaT <- cef_mock_metaT[cef_mock_metaT != 0]
clinda_630_metaT <- read.delim(clinda_630_metaT, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
clinda_630_metaT <- clinda_630_metaT[clinda_630_metaT != 0]
clinda_mock_metaT <- read.delim(clinda_mock_metaT, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
clinda_mock_metaT <- clinda_mock_metaT[clinda_mock_metaT != 0]
strep_630_metaT <- read.delim(strep_630_metaT, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
strep_630_metaT <- strep_630_metaT[strep_630_metaT != 0]
strep_mock_metaT <- read.delim(strep_mock_metaT, sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))[,1]
strep_mock_metaT <- strep_mock_metaT[strep_mock_metaT != 0]
noabx_mock_metaT <- read.delim(noabx_mock_metaT, sep='\t', header=FALSE, row.names=1)[,1]
noabx_mock_metaT <- noabx_mock_metaT[noabx_mock_metaT != 0]

#--------------------------------------------------------------------------------------------------#

# Stepwise rarefaction analysis (order of magnitude)
stepRarefy <- function(abundVect){
  options(warn=-1)
  rareVect <- c()
  subVect <- c(1)
  sub_max <- 10^ceiling(log10(signif(as.numeric(sum(abundVect)), digits=1)))
  current <- 1
  while (current != sub_max)
  {
    current <- current * 10
    subVect <- c(subVect, current)
  }
  for (x in 1:length(subVect)) {
    rareVect[x] <- as.numeric(sum(as.vector(rrarefy(abundVect, sample=subVect[x])) != 0))
  }
  options(warn=0)
  return(rareVect)
}

# Perform rarefaction
cef_630_metaG <- stepRarefy(cef_630_metaG)
cef_mock_metaG <- stepRarefy(cef_mock_metaG)
strep_630_metaG <- stepRarefy(strep_630_metaG)
strep_mock_metaG <- stepRarefy(strep_mock_metaG)
clinda_630_metaG <- stepRarefy(clinda_630_metaG)
clinda_mock_metaG <- stepRarefy(clinda_mock_metaG)
noabx_mock_metaG <- stepRarefy(noabx_mock_metaG)
cef_630_metaT <- stepRarefy(cef_630_metaT)
cef_mock_metaT <- stepRarefy(cef_mock_metaT)
strep_630_metaT <- stepRarefy(strep_630_metaT)
strep_mock_metaT <- stepRarefy(strep_mock_metaT)
clinda_630_metaT <- stepRarefy(clinda_630_metaT)
clinda_mock_metaT <- stepRarefy(clinda_mock_metaT)
noabx_mock_metaT <- stepRarefy(noabx_mock_metaT)

# Log transform
cef_630_metaG <- log10(cef_630_metaG)
cef_mock_metaG <- log10(cef_mock_metaG)
strep_630_metaG <- log10(strep_630_metaG)
strep_mock_metaG <- log10(strep_mock_metaG)
clinda_630_metaG <- log10(clinda_630_metaG)
clinda_mock_metaG <- log10(clinda_mock_metaG)
noabx_mock_metaG <- log10(noabx_mock_metaG)
cef_630_metaT <- log10(cef_630_metaT)
cef_mock_metaT <- log10(cef_mock_metaT)
strep_630_metaT <- log10(strep_630_metaT)
strep_mock_metaT <- log10(strep_mock_metaT)
clinda_630_metaT <- log10(clinda_630_metaT)
clinda_mock_metaT <- log10(clinda_mock_metaT)
noabx_mock_metaT <- log10(noabx_mock_metaT)

#--------------------------------------------------------------------------------------------------#

# Plot collector's curves
pdf(file='results/supplement/figures/figure_S6EF.pdf', width=10, height=5)
layout(matrix(c(1,2),
              nrow=1, ncol=2, byrow=TRUE))

# Metagenomes
par(mar=c(4,4,1,2), mgp=c(2.5, 0.75, 0), las=1, xaxs='i', yaxs='i')
plot(1, type='n', xlim=c(1,10), ylim=c(0,6), xaxt='n', yaxt='n',
     xlab='Sampling Depth', ylab='Detected Genes')
axis(1, at=seq(1,12,2), label=c(1,100,10000,1000000,100000000,1000000000), cex.axis=0.9)
minor.ticks.axis(1, 10, mn=0, mx=10)
axis(2, at=seq(0,6,2), label=c(1,100,10000,1000000), cex.axis=0.9, las=1)
minor.ticks.axis(2, 10, mn=0, mx=6)
legend('topleft', 'Metagenomes', bty='n', cex=1.2)
mtext('E', side=2, line=2, las=2, adj=2, padj=-9, cex=1.7, font=2)
legend('right', legend=c(expression(italic('C. difficile')),'Mock'), pch=c(17,19), pt.cex=1.5)
legend('bottomright', legend=c('Streptomycin-pretreated','Cefoperazone-pretreated','Clindamycin-pretreated'), 
       col=c(strep_col,cef_col,clinda_col), pch=22, pt.cex=0, lwd=3)
box()
# Add lines
lines(cef_630_metaG, col=cef_col, type='p', pch=17, cex=1.5)
lines(cef_630_metaG, col=cef_col, lwd=2.5)
lines(cef_mock_metaG, col=cef_col, type='p', pch=19, cex=1.5)
lines(cef_mock_metaG, col=cef_col, lwd=2.5)
lines(strep_630_metaG, col=strep_col, type='p', pch=17, cex=1.5)
lines(strep_630_metaG, col=strep_col, lwd=2.5)
lines(strep_mock_metaG, col=strep_col, type='p', pch=19, cex=1.5)
lines(strep_mock_metaG, col=strep_col, lwd=2.5)
lines(clinda_630_metaG, col=clinda_col, type='p', pch=17, cex=1.5)
lines(clinda_630_metaG, col=clinda_col, lwd=2.5)
lines(clinda_mock_metaG, col=clinda_col, type='p', pch=19, cex=1.5)
lines(clinda_mock_metaG, col=clinda_col, lwd=2.5)

# Metatranscriptomes
par(mar=c(4,4,1,2), mTp=c(2.5, 0.75, 0), las=1, xaxs='i', yaxs='i')
plot(1, type='n', xlim=c(1,10), ylim=c(0,6), xaxt='n', yaxt='n',
     xlab='Sampling Depth', ylab='Detected Genes')
axis(1, at=seq(1,12,2), label=c(1,100,10000,1000000,100000000,1000000000), cex.axis=0.9)
minor.ticks.axis(1, 10, mn=0, mx=10)
axis(2, at=seq(0,6,2), label=c(1,100,10000,1000000), cex.axis=0.9, las=1)
minor.ticks.axis(2, 10, mn=0, mx=6)
legend('topleft', 'Metatranscriptomes', bty='n', cex=1.2)
mtext('F', side=2, line=2, las=2, adj=2, padj=-9, cex=1.7, font=2)
legend('right', legend=c(expression(italic('C. difficile')),'Mock'), pch=c(17,19), pt.cex=1.5)
legend('bottomright', legend=c('Streptomycin-pretreated','Cefoperazone-pretreated','Clindamycin-pretreated'), 
       col=c(strep_col,cef_col,clinda_col), pch=22, pt.cex=0, lwd=3)
box()
# Add lines
lines(cef_630_metaT, col=cef_col, type='p', pch=17, cex=1.5)
lines(cef_630_metaT, col=cef_col, lwd=2.5)
lines(cef_mock_metaT, col=cef_col, type='p', pch=19, cex=1.5)
lines(cef_mock_metaT, col=cef_col, lwd=2.5)
lines(strep_630_metaT, col=strep_col, type='p', pch=17, cex=1.5)
lines(strep_630_metaT, col=strep_col, lwd=2.5)
lines(strep_mock_metaT, col=strep_col, type='p', pch=19, cex=1.5)
lines(strep_mock_metaT, col=strep_col, lwd=2.5)
lines(clinda_630_metaT, col=clinda_col, type='p', pch=17, cex=1.5)
lines(clinda_630_metaT, col=clinda_col, lwd=2.5)
lines(clinda_mock_metaT, col=clinda_col, type='p', pch=19, cex=1.5)
lines(clinda_mock_metaT, col=clinda_col, lwd=2.5)

dev.off()

#--------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()

