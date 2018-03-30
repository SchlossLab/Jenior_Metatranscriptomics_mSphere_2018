
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')

# Metabolomes
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'
aminovalerate <- 'exploratory/aminovalerate.tsv'

# Input Metadata
metadata <- 'data/metadata.tsv'

# Output files
plot_abc <- 'results/figures/figure_4abc.pdf'
plot_d <- 'results/figures/figure_4d.pdf'
plot_e <- 'results/figures/figure_4e.pdf'
plot_f <- 'results/figures/figure_4f.pdf'
plot_g <- 'results/figures/figure_4g.pdf'

#----------------#

# Read in data
# Metabolomes
metabolome <- read.delim(metabolome, sep='\t', header=TRUE)
aminovalerate <- read.delim(aminovalerate, sep='\t', header=T)

# Metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data
# Metadata
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL

# Metabolomes
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.data.frame(t(metabolome))

aminovalerate <- subset(aminovalerate, abx != 'germfree')

#-------------------------------------------------------------------------------------------------------------------------#

# Ordination
# Prep data
metabolome <- clean_merge(metadata, metabolome)
metabolome <- subset(metabolome, abx != 'germfree')
cef_metabolome <- subset(metabolome, abx == 'cefoperazone')
cef_metabolome$abx <- NULL
cef_metabolome$susceptibility <- NULL
clinda_metabolome <- subset(metabolome, abx == 'clindamycin')
clinda_metabolome$abx <- NULL
clinda_metabolome$susceptibility <- NULL
strep_metabolome <- subset(metabolome, abx == 'streptomycin')
strep_metabolome$abx <- NULL
strep_metabolome$susceptibility <- NULL

# Calculate significant differences
strep_p <- round(adonis(strep_metabolome[,2:ncol(strep_metabolome)]~factor(strep_metabolome$infection), data=strep_metabolome, permutations=10000, method='jaccard')$aov.tab[[6]][1], 3)
strep_metabolome$infection <- NULL
cef_p <- round(adonis(cef_metabolome[,2:ncol(cef_metabolome)]~factor(cef_metabolome$infection), data=cef_metabolome, permutations=10000, method='jaccard')$aov.tab[[6]][1], 3)
cef_metabolome$infection <- NULL
clinda_p <- round(adonis(clinda_metabolome[,2:ncol(clinda_metabolome)]~factor(clinda_metabolome$infection), data=clinda_metabolome, permutations=10000, method='jaccard')$aov.tab[[6]][1], 3)
clinda_metabolome$infection <- NULL

# Calculate axes and merge with metadata
cef_metabolome_nmds <- metaMDS(cef_metabolome, k=2, trymax=100)$points
cef_metabolome_nmds <- clean_merge(metadata, cef_metabolome_nmds)
clinda_metabolome_nmds <- metaMDS(clinda_metabolome, k=2, trymax=100)$points
clinda_metabolome_nmds <- clean_merge(metadata, clinda_metabolome_nmds)
strep_metabolome_nmds <- metaMDS(strep_metabolome, k=2, trymax=100)$points
strep_metabolome_nmds <- clean_merge(metadata, strep_metabolome_nmds)
rm(cef_metabolome, clinda_metabolome, strep_metabolome)

# Subset to points for plot
cef_metabolome_nmds$MDS2 <- cef_metabolome_nmds$MDS2 - 0.02
clinda_metabolome_nmds$MDS1 <- clinda_metabolome_nmds$MDS1 + 0.1
clinda_metabolome_nmds$MDS2 <- clinda_metabolome_nmds$MDS2 + 0.02
cef_metabolome_nmds_630 <- subset(cef_metabolome_nmds, infection == '630')
cef_metabolome_nmds_mock <- subset(cef_metabolome_nmds, infection == 'mock')
clinda_metabolome_nmds_630 <- subset(clinda_metabolome_nmds, infection == '630')
clinda_metabolome_nmds_mock <- subset(clinda_metabolome_nmds, infection == 'mock')
strep_metabolome_nmds_630 <- subset(strep_metabolome_nmds, infection == '630')
strep_metabolome_nmds_mock <- subset(strep_metabolome_nmds, infection == 'mock')

# Calculate centroids
cef_metabolome_centoids <- aggregate(cbind(cef_metabolome_nmds$MDS1,cef_metabolome_nmds$MDS2)~cef_metabolome_nmds$infection, data=cef_metabolome_nmds, mean)
clinda_metabolome_centoids <- aggregate(cbind(clinda_metabolome_nmds$MDS1,clinda_metabolome_nmds$MDS2)~clinda_metabolome_nmds$infection, data=clinda_metabolome_nmds, mean)
strep_metabolome_centoids <- aggregate(cbind(strep_metabolome_nmds$MDS1,strep_metabolome_nmds$MDS2)~strep_metabolome_nmds$infection, data=strep_metabolome_nmds, mean)

# Amnovalerate data
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

#-------------------------------------------------------------------------------------------------------------------------#

# Feature selection
# Separate groups
cef_metabolome <- subset(metabolome, abx == 'cefoperazone')
cef_metabolome$abx <- NULL
cef_metabolome$susceptibility <- NULL
cef_metabolome$infection <- factor(cef_metabolome$infection)
clinda_metabolome <- subset(metabolome, abx == 'clindamycin')
clinda_metabolome$abx <- NULL
clinda_metabolome$susceptibility <- NULL
clinda_metabolome$infection <- factor(clinda_metabolome$infection)
strep_metabolome <- subset(metabolome, abx == 'streptomycin')
strep_metabolome$abx <- NULL
strep_metabolome$susceptibility <- NULL
strep_metabolome$infection <- factor(strep_metabolome$infection)
rm(metabolome)

# AUCRF
colnames(cef_metabolome) <- make.names(colnames(cef_metabolome))
cef_rf <- aucrfInfection(cef_metabolome)
colnames(clinda_metabolome) <- make.names(colnames(clinda_metabolome))
clinda_rf <- aucrfInfection(clinda_metabolome)
colnames(strep_metabolome) <- make.names(colnames(strep_metabolome))
strep_rf <- aucrfInfection(strep_metabolome)

# Get OOB error
cef_oob <- cef_rf$RFopt
cef_oob <- cef_oob$err.rate
cef_oob <- as.character(round(median(cef_oob[,1]) * 100, 3))
strep_oob <- strep_rf$RFopt
strep_oob <- strep_oob$err.rate
strep_oob <- as.character(round(median(strep_oob[,1]) * 100, 3))

# Get optimal features
strep_rf <- c('infection', as.character(OptimalSet(strep_rf)$Name))
strep_metabolome <- strep_metabolome[,strep_rf]
colnames(strep_metabolome) <- c("infection","Prolyl-4-hydroxyproline","Arabitol / Xylitol","5-Aminovalerate","Arabonate / Xylonate","Mevalonate")
inf_strep_metabolome <- subset(strep_metabolome, infection == '630')
inf_strep_metabolome$infection<- NULL
mock_strep_metabolome <- subset(strep_metabolome, infection == 'mock')
mock_strep_metabolome$infection<- NULL
rm(strep_metabolome)
cef_rf <- c('infection', as.character(OptimalSet(cef_rf)$Name))
cef_metabolome <- cef_metabolome[,cef_rf]
colnames(cef_metabolome) <- c("infection","Ectoine","N-Acetylglucosaminylasparagine","N-Acetylthreonine","beta-Cryptoxanthin","Adenosine")
inf_cef_metabolome <- subset(cef_metabolome, infection == '630')
inf_cef_metabolome$infection<- NULL
mock_cef_metabolome <- subset(cef_metabolome, infection == 'mock')
mock_cef_metabolome$infection<- NULL
rm(cef_metabolome)

# Run AUCRF again on features from first run through
clinda_rf <- c('infection', as.character(OptimalSet(clinda_rf)$Name))
clinda_metabolome <- clinda_metabolome[,clinda_rf]
clinda_rf <- aucrfInfection(clinda_metabolome)
# Get OOB error - clinda
clinda_oob <- clinda_rf$RFopt
clinda_oob <- clinda_oob$err.rate
clinda_oob <- as.character(round(median(clinda_oob[,1]) * 100, 3))

# Generates a slightly overfit model with low OOD and 10 features
clinda_rf <- c('infection', as.character(OptimalSet(clinda_rf)$Name))
clinda_metabolome <- clinda_metabolome[,clinda_rf]
colnames(clinda_metabolome) <- c('infection',"Adenine","Nicotinate","N6-Acetyllysine","Gulonate","N-Palmitoylsphingosine","2-Deoxyadenosine","1-Oleoylglycerophosphocholine","N-Acetylleucine","5-Aminovalerate","Cysteinylglycine")
inf_clinda_metabolome <- subset(clinda_metabolome, infection == '630')
inf_clinda_metabolome$infection<- NULL
mock_clinda_metabolome <- subset(clinda_metabolome, infection == 'mock')
mock_clinda_metabolome$infection<- NULL
rm(metadata, clinda_metabolome)

# Find significant differences
cef_pvalues <- c()
for (i in 1:ncol(inf_cef_metabolome)){cef_pvalues[i] <- wilcox.test(inf_cef_metabolome[,i], mock_cef_metabolome[,i], exact=FALSE)$p.value}
cef_pvalues <- p.adjust(cef_pvalues, method='BH')
clinda_pvalues <- c()
for (i in 1:ncol(inf_clinda_metabolome)){clinda_pvalues[i] <- wilcox.test(inf_clinda_metabolome[,i], mock_clinda_metabolome[,i], exact=FALSE)$p.value}
pvalues <- p.adjust(clinda_pvalues, method='BH')
strep_pvalues <- c()
for (i in 1:ncol(inf_strep_metabolome)){strep_pvalues[i] <- wilcox.test(inf_strep_metabolome[,i], mock_strep_metabolome[,i], exact=FALSE)$p.value}
strep_pvalues <- p.adjust(strep_pvalues, method='BH')
rm(clinda_rf, cef_rf, strep_rf)

#-------------------------------------------------------------------------------------------------------------------------#

# Generate figure panels
# Ordinations
pdf(file=plot_abc, width=12, height=4)
layout(matrix(c(1,2,3),
              nrow=1, ncol=3, byrow=TRUE))
par(mar=c(4.5,4,1.5,1), las=1, mgp=c(2.8,0.75,0))
#Streptomycin - Fig. 4a
plot(x=strep_metabolome_nmds$MDS1, y=strep_metabolome_nmds$MDS2, xlim=c(-0.25,0.25), ylim=c(-0.15,0.15),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('A', side=2, line=2, las=2, adj=1.4, padj=-8.3, cex=1.6, font=2)
legend('topleft', legend='Streptomycin-pretreated', pch=1, cex=1.5, pt.cex=0, bty='n')
segments(x0=strep_metabolome_nmds_630$MDS1, y0=strep_metabolome_nmds_630$MDS2, x1=strep_metabolome_centoids[1,2], y1=strep_metabolome_centoids[1,3], col='gray30')
segments(x0=strep_metabolome_nmds_mock$MDS1, y0=strep_metabolome_nmds_mock$MDS2, x1=strep_metabolome_centoids[2,2], y1=strep_metabolome_centoids[2,3], col='gray30')
points(x=strep_metabolome_nmds_630$MDS1, y=strep_metabolome_nmds_630$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=strep_metabolome_nmds_mock$MDS1, y=strep_metabolome_nmds_mock$MDS2, bg=strep_col, pch=24, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Mock vs Infected', as.expression(bquote(paste(italic('p'),' = 0.039*')))), 
       pch=1, cex=1.4, pt.cex=0, bty='n')
legend('bottomright', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(19,17), cex=1.2, pt.cex=2)
box()
mtext('D', side=2, line=2, las=2, adj=1.4, padj=-6, cex=1.6, font=2)

# Cefoperazone - Fig. 4b
plot(x=cef_metabolome_nmds$MDS1, y=cef_metabolome_nmds$MDS2, xlim=c(-0.15,0.15), ylim=c(-0.15,0.15),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('B', side=2, line=2, las=2, adj=1.4, padj=-8.3, cex=1.6, font=2)
legend('topleft', legend='Cefoperazone-pretreated', pch=1, cex=1.5, pt.cex=0, bty='n')
segments(x0=cef_metabolome_nmds_630$MDS1, y0=cef_metabolome_nmds_630$MDS2, x1=cef_metabolome_centoids[1,2], y1=cef_metabolome_centoids[1,3], col='gray30')
segments(x0=cef_metabolome_nmds_mock$MDS1, y0=cef_metabolome_nmds_mock$MDS2, x1=cef_metabolome_centoids[2,2], y1=cef_metabolome_centoids[2,3], col='gray30')
points(x=cef_metabolome_nmds_630$MDS1, y=cef_metabolome_nmds_630$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=cef_metabolome_nmds_mock$MDS1, y=cef_metabolome_nmds_mock$MDS2, bg=cef_col, pch=24, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Mock vs Infected', as.expression(bquote(paste(italic('p'),' = 0.016*')))), 
       pch=1, cex=1.4, pt.cex=0, bty='n')
legend('bottomright', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(19,17), cex=1.2, pt.cex=2)
box()
mtext('E', side=2, line=2, las=2, adj=1.4, padj=-6, cex=1.6, font=2)

# Clindamycin - Fig. 4c
plot(x=clinda_metabolome_nmds$MDS1-0.075, y=clinda_metabolome_nmds$MDS2, xlim=c(-0.15,0.15), ylim=c(-0.1,0.1),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex.axis=1.2, cex.lab=1.2)
mtext('C', side=2, line=2, las=2, adj=1.4, padj=-8.3, cex=1.6, font=2)
legend('topleft', legend='Clindamycin-pretreated', pch=1, cex=1.5, pt.cex=0, bty='n')
segments(x0=clinda_metabolome_nmds_630$MDS1-0.075, y0=clinda_metabolome_nmds_630$MDS2, x1=clinda_metabolome_centoids[1,2]-0.075, y1=clinda_metabolome_centoids[1,3], col='gray30')
segments(x0=clinda_metabolome_nmds_mock$MDS1-0.0751, y0=clinda_metabolome_nmds_mock$MDS2, x1=clinda_metabolome_centoids[2,2]-0.075, y1=clinda_metabolome_centoids[2,3], col='gray30')
points(x=clinda_metabolome_nmds_630$MDS1-0.075, y=clinda_metabolome_nmds_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=clinda_metabolome_nmds_mock$MDS1-0.075, y=clinda_metabolome_nmds_mock$MDS2, bg=clinda_col, pch=24, cex=2, lwd=1.2)
legend('bottomleft', legend=c('Mock vs Infected',as.expression(bquote(paste(italic('p'),' = 0.127 n.s.')))), 
       pch=1, cex=1.4, pt.cex=0, bty='n')
legend('bottomright', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(19,17), cex=1.2, pt.cex=2)
box()
mtext('F', side=2, line=2, las=2, adj=1.4, padj=-6, cex=1.6, font=2)
mtext('G', side=2, line=2, las=2, adj=1.4, padj=-4.5, cex=1.6, font=2)
dev.off()

#----------------#

# Feature Selection
# Strep Infected
metabolite_stripchart(plot_d, inf_strep_metabolome, mock_strep_metabolome, strep_pvalues,  
                      strep_oob, 'Infected', 'Mock', 'mediumorchid4', 'white', 'Streptomycin-pretreated', strep_col)
# Cef Infected
metabolite_stripchart(plot_e, inf_cef_metabolome, mock_cef_metabolome, cef_pvalues,  
                      cef_oob, 'Infected', 'Mock', 'mediumorchid4', 'white', 'Cefoperazone-pretreated', cef_col)
# Clinda Infected
metabolite_stripchart(plot_f, inf_clinda_metabolome, mock_clinda_metabolome, clinda_pvalues,  
                      clinda_oob, 'Infected', 'Mock', 'mediumorchid4', 'white', 'Clindamycin-pretreated', clinda_col)

#----------------#

# Aminovalerate
pdf(file=plot_g, width=6, height=4)
par(mar=c(3.5,5,1.5,1), xpd=FALSE, las=1, mgp=c(3,0.7,0))
stripchart(substrate~infection, data=aminovalerate_untreated, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,6), xlim=c(0.5,10.5),
           cex=1.5, ylab='', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(substrate~infection, data=aminovalerate_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=strep_col, ylim=c(0,6), xlim=c(0.5,10.5),
           cex=1.5, ylab='', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=aminovalerate_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=cef_col, ylim=c(0,6), xlim=c(0.5,10.5),
           cex=1.5, ylab='', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=aminovalerate_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=clinda_col, ylim=c(0,6), xlim=c(0.5,10.5),
           cex=1.5, ylab='', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
axis(side=2, at=c(0:6), labels=c('0.0','1.0','2.0','3.0', '4.0','5.0','6.0'), cex.axis=1.2)
box()
mtext(text=expression(paste('Scaled Intensity (',log[10],')')), side=2, cex=1.2, las=0, padj=-2.5)
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext(c('CDI:','Group:'), side=1, at=-0.7, padj=c(0.3,2.5), cex=0.7)
mtext(c('-','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10), padj=0.3, cex=1.1)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin'), side=1, 
      at=c(1,3.5,6.5,9.5), padj=2, cex=0.9)
legend('topright', legend='5-aminovalerate', pt.cex=0, bty='n', cex=0.9)
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4),
         y0=c(median(aminovalerate_untreated[,2]),
              median(subset(aminovalerate_strep, infection=='mock')[,2]), median(subset(aminovalerate_strep, infection=='infected')[,2]),
              median(subset(aminovalerate_cef, infection=='mock')[,2]), median(subset(aminovalerate_cef, infection=='infected')[,2]),
              median(subset(aminovalerate_clinda, infection=='mock')[,2]), median(subset(aminovalerate_clinda, infection=='infected')[,2])),
         y1=c(median(aminovalerate_untreated[,2]),
              median(subset(aminovalerate_strep, infection=='mock')[,2]), median(subset(aminovalerate_strep, infection=='infected')[,2]),
              median(subset(aminovalerate_cef, infection=='mock')[,2]), median(subset(aminovalerate_cef, infection=='infected')[,2]),
              median(subset(aminovalerate_clinda, infection=='mock')[,2]), median(subset(aminovalerate_clinda, infection=='infected')[,2])),
         lwd=3)
segments(x0=c(3,6,9), y0=5, x1=c(4,7,10), y1=5, lwd=2)
text(x=c(3.5,6.5,9.5), y=5.2, '*', font=2, cex=2)
mtext(c(rep('*',5),'n.s.'), side=3, adj=c(0.25,0.35,
                                0.54,0.64,
                                0.83,0.95), padj=c(rep(0.4,5),-0.1), font=2, cex=c(rep(1.6,5),1.1), col='chartreuse3') # Untreated vs Mock significance
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

