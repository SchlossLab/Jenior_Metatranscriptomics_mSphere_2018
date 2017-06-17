# Start with clear environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# CFU over time
cfu <- 'data/cfu_time.tsv'

# Define output files
plot_file <- 'results/supplement/figures/figure_S1.pdf'

#---------------------------------------------------------------------------------------------------------------#

# Load in data
cfu <- read.delim(cfu, sep='\t', header=TRUE)

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

#---------------------------------------------------------------------------------------------------------------#

# Generate figure
pdf(file=plot_file, width=9, height=7)

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
legend('right', legend=c('Streptomycin (5 mg/ml)', 'Cefoperazone (0.5 mg/ml)', 'Clindamycin (10 mg/kg', 'No Antibiotics'),
       pch=16, col=c(strep_col, cef_col, clinda_col, noabx_col), cex=0.9, pt.cex=1.5)

dev.off()

#---------------------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()

