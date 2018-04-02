
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')

# Metabolomes
aminovalerate <- 'exploratory/aminovalerate.tsv'

# Output files
plot_S2D <- 'results/supplement/figures/figure_S2D.pdf'

#----------------#

# Read in data
aminovalerate <- read.delim(aminovalerate, sep='\t', header=T)

# Subset data
aminovalerate <- subset(aminovalerate, abx == 'germfree')
aminovalerate$infection <- as.factor(aminovalerate$infection)
aminovalerate$aminovalerate <- as.numeric(aminovalerate$aminovalerate)
aminovalerate$abx <- NULL

# Test differences
amv_pval <- round(wilcox.test(as.numeric(subset(aminovalerate, infection == 'infected')$aminovalerate), 
                              as.numeric(subset(aminovalerate, infection == 'mock')$aminovalerate), 
                              exact=FALSE)$p.value, 4)

# Calculate medians
infected_med <- median(as.numeric(subset(aminovalerate, infection == 'infected')$aminovalerate))
mock_med <- median(as.numeric(subset(aminovalerate, infection == 'mock')$aminovalerate))

#----------------#

# Plot data
pdf(file=plot_S2D, width=3, height=4)
par(mar=c(3,3,1,1), xpd=FALSE, las=1, mgp=c(2,0.7,0))
stripchart(aminovalerate~infection, data=aminovalerate, vertical=T, pch=21, 
           xaxt='n', bg=gf_col, ylim=c(0,6), lwd=2,
           cex=1.7, ylab='Scaled Intensity (Log10)', method='jitter', jitter=0.1, cex.lab=1.2)
box(lwd=2)
mtext('CDI:', side=1, at=0.5, padj=0.3, cex=1.1)
mtext(c('+','-'), side=1, at=c(1,2), padj=0.3, cex=1.6)
mtext('Germfree', side=1, at=1.5, padj=2, cex=1.5)
segments(x0=c(0.8,1.8), x1=c(1.2,2.2), # Medians
         y0=c(infected_med, mock_med),
         y1=c(infected_med, mock_med),
         lwd=3)
segments(x0=1, x1=2, y0=5.5, y1=5.5, lwd=3) # Significance
text(x=1.5, y=5.8, '*', font=2, cex=2.5) # Significance
mtext('D', side=2, adj=2.5, padj=-7, cex=1.8, font=2)
dev.off()

#----------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()
