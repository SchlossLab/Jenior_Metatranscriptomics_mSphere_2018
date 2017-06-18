
# Start with blank slate
rm(list=ls())
gc()

# Load dependency
if ('wesanderson' %in% installed.packages()[,"Package"] == FALSE){
  install.packages(as.character('wesanderson'), quiet=TRUE);
}
library('wesanderson', verbose=FALSE, character.only=TRUE)

# Select files
aminovalerate <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/exploratory/aminovalerate.tsv'
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/exploratory/aminovalerate.pdf'

# Read in data
aminovalerate <- read.delim(aminovalerate, sep='\t', header=T)

# Subset metabolites
aminovalerate_untreated <- subset(aminovalerate, abx == 'none')
aminovalerate_untreated$abx <- NULL
colnames(aminovalerate_untreated) <- c('infection', 'substrate')
aminovalerate_untreated$infection <- factor(aminovalerate_untreated$infection, levels=c('mock','infected'))
aminovalerate_cef <- subset(aminovalerate, abx == 'cefoperazone')
aminovalerate_cef$abx <- NULL
colnames(aminovalerate_cef) <- c('infection', 'substrate')
aminovalerate_cef$infection <- factor(aminovalerate_cef$infection, levels=c('mock','infected'))
aminovalerate_strep <- subset(aminovalerate, abx == 'streptomycin')
aminovalerate_strep$abx <- NULL
colnames(aminovalerate_strep) <- c('infection', 'substrate')
aminovalerate_strep$infection <- factor(aminovalerate_strep$infection, levels=c('mock','infected'))
aminovalerate_clinda <- subset(aminovalerate, abx == 'clindamycin')
aminovalerate_clinda$abx <- NULL
colnames(aminovalerate_clinda) <- c('infection', 'substrate')
aminovalerate_clinda$infection <- factor(aminovalerate_clinda$infection, levels=c('mock','infected'))
aminovalerate_gf <- subset(aminovalerate, abx == 'germfree')
aminovalerate_gf$abx <- NULL
colnames(aminovalerate_gf) <- c('infection', 'substrate')
aminovalerate_gf$infection <- factor(aminovalerate_gf$infection, levels=c('mock','infected'))
rm(aminovalerate)
rm(metabolome)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Calculate significant differences

# Untreated vs Everything
p.adjust(c(wilcox.test(subset(aminovalerate_untreated, infection=='mock')[,2], subset(aminovalerate_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(aminovalerate_untreated, infection=='mock')[,2], subset(aminovalerate_strep, infection=='infected')[,2], exact=F)$p.value,
           wilcox.test(subset(aminovalerate_untreated, infection=='mock')[,2], subset(aminovalerate_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(aminovalerate_untreated, infection=='mock')[,2], subset(aminovalerate_cef, infection=='infected')[,2], exact=F)$p.value,
           wilcox.test(subset(aminovalerate_untreated, infection=='mock')[,2], subset(aminovalerate_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(aminovalerate_untreated, infection=='mock')[,2], subset(aminovalerate_clinda, infection=='infected')[,2], exact=F)$p.value,
           wilcox.test(subset(aminovalerate_untreated, infection=='mock')[,2], subset(aminovalerate_gf, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(aminovalerate_untreated, infection=='mock')[,2], subset(aminovalerate_gf, infection=='infected')[,2], exact=F)$p.value), method='BH')
# 0.0006596717 0.0010563576 0.0006596717 0.0016896741 0.0006596717 0.2509982704 0.0006596717 0.0006596717

#------------------#

# Mock vs Infected
p.adjust(c(wilcox.test(subset(aminovalerate_strep, infection=='infected')[,2], subset(aminovalerate_strep, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(aminovalerate_cef, infection=='infected')[,2], subset(aminovalerate_cef, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(aminovalerate_clinda, infection=='infected')[,2], subset(aminovalerate_clinda, infection=='mock')[,2], exact=F)$p.value,
           wilcox.test(subset(aminovalerate_gf, infection=='infected')[,2], subset(aminovalerate_gf, infection=='mock')[,2], exact=F)$p.value), method='BH')
# 0.0005497264 0.0005497264 0.0104442711 0.0005497264

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
pdf(file=plot_file, width=8, height=4)
par(mar=c(3.5,5,1.5,1), xpd=FALSE, las=1, mgp=c(3,0.7,0))

# N-aminovalerate
stripchart(substrate~infection, data=aminovalerate_untreated, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,6), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intesity', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(substrate~infection, data=aminovalerate_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[1], ylim=c(0,6), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=aminovalerate_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[3], ylim=c(0,6), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=aminovalerate_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=wes_palette('FantasticFox')[5], ylim=c(0,6), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=aminovalerate_gf, vertical=T, pch=19, at=c(12,13),
           xaxt='n', yaxt='n', col='forestgreen', ylim=c(0,6), xlim=c(0.5,13.5),
           cex=1.5, ylab='Scaled Intensity', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
axis(side=2, at=c(0:6), labels=c('0.0','1.0','2.0','3.0', '4.0','5.0','6.0'), cex.axis=1.2)
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext(c('CDI:','Group:'), side=1, at=-0.7, padj=c(0.3,2.5), cex=0.7)
mtext(c('-','-','+','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10,12,13), padj=0.3, cex=1.1)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin','ex-Germfree'), side=1, 
      at=c(1,3.5,6.5,9.5,12.5), padj=2, cex=0.9)
legend('topright', legend='5-Aminovalerate', pt.cex=0, bty='n')
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6,11.6,12.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4,12.4,13.4),
         y0=c(median(aminovalerate_untreated[,2]),
              median(subset(aminovalerate_strep, infection=='mock')[,2]), median(subset(aminovalerate_strep, infection=='infected')[,2]),
              median(subset(aminovalerate_cef, infection=='mock')[,2]), median(subset(aminovalerate_cef, infection=='infected')[,2]),
              median(subset(aminovalerate_clinda, infection=='mock')[,2]), median(subset(aminovalerate_clinda, infection=='infected')[,2]),
              median(subset(aminovalerate_gf, infection=='mock')[,2]), median(subset(aminovalerate_gf, infection=='infected')[,2])), 
         y1=c(median(aminovalerate_untreated[,2]),
              median(subset(aminovalerate_strep, infection=='mock')[,2]), median(subset(aminovalerate_strep, infection=='infected')[,2]),
              median(subset(aminovalerate_cef, infection=='mock')[,2]), median(subset(aminovalerate_cef, infection=='infected')[,2]),
              median(subset(aminovalerate_clinda, infection=='mock')[,2]), median(subset(aminovalerate_clinda, infection=='infected')[,2]),
              median(subset(aminovalerate_gf, infection=='mock')[,2]), median(subset(aminovalerate_gf, infection=='infected')[,2])),
         lwd=3)
segments(x0=c(3,6,9,12), y0=5, x1=c(4,7,10,13), y1=5, lwd=2)
text(x=c(3.5,6.5,9.5,12.5), y=5.2, '*', font=2, cex=2)
mtext(rep('*',7), side=3, adj=c(0.21,0.28,
                                0.43,0.5,
                                0.645,
                                0.863,0.933), padj=0.4, font=2, cex=1.5, col='gray40') # Untreated vs Mock significance

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
detach('package:wesanderson', character.only = TRUE)
rm(list=ls())
gc()

