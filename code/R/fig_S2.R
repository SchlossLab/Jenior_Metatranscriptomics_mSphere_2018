
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')

# Output plot name
plot_a <- 'results/supplement/figures/figure_S2a.pdf'
plot_b <- 'results/supplement/figures/figure_S2b.pdf'
plot_c <- 'results/supplement/figures/figure_S2c.pdf'

# Input Metabolomes
metabolome_file <- 'data/metabolome/scaled_intensities.log10.tsv'

# Input Metadata
metadata_file <- 'data/metadata.tsv'

#----------------#

# Read in data

# Metabolomes
metabolome <- read.delim(metabolome_file, sep='\t', header=TRUE)
rm(metabolome_file)

# Metadata
metadata <- read.delim(metadata_file, sep='\t', header=T, row.names=1)
rm(metadata_file)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data

# Metadata
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL

# Metabolomes
metabolome$BIOCHEMICAL <- gsub('_',' ', metabolome$BIOCHEMICAL)
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.data.frame(t(metabolome))

#-------------------------------------------------------------------------------------------------------------------------#

# Get bile acids
bile_acids <- c('beta-muricholate','6-beta-hydroxylithocholate','lithocholate','deoxycholate',
                'lithocholate [6-oxo or 7-keto]','gamma-muricholate','3-dehydrocholate','cholate',
                '7,12-diketolithocholate','12-dehydrocholate','hyodeoxycholate','taurocholate',
                'alpha-muricholate','tauroursodeoxycholate','taurocholenate sulfate')
bile_metabolome <- metabolome[, colnames(metabolome) %in% bile_acids]  
carbs <- c('fructose', 'ribulose/xylulose', 'mannitol/sorbitol', 'arabitol/xylitol', 
           'ribitol', 'sucrose')
carb_metabolome <- metabolome[, colnames(metabolome) %in% carbs]  
amino_acids <- c('glycine', 'proline', 'cysteine', 'isoleucine', 'valine', 
                 'leucine', 'tryptophan')
amino_metabolome <- metabolome[, colnames(metabolome) %in% amino_acids]  

# Separate groups
bile_metabolome <- clean_merge(metadata, bile_metabolome)
abx_bile_metabolome <- subset(bile_metabolome, abx != 'germfree')
abx_bile_metabolome$infection <- NULL
abx_bile_metabolome$susceptibility <- NULL
abx_bile_metabolome$abx <- factor(abx_bile_metabolome$abx)
bile_metabolome$infection <- NULL
bile_metabolome$abx <- NULL
bile_metabolome$abx <- factor(bile_metabolome$susceptibility)
noabx_bile_metabolome <- subset(abx_bile_metabolome, abx == 'none')
noabx_bile_metabolome$abx <- NULL
cef_abx_bile_metabolome <- subset(abx_bile_metabolome, abx == 'cefoperazone')
cef_abx_bile_metabolome$abx <- NULL
clinda_abx_bile_metabolome <- subset(abx_bile_metabolome, abx == 'clindamycin')
clinda_abx_bile_metabolome$abx <- NULL
strep_abx_bile_metabolome <- subset(abx_bile_metabolome, abx == 'streptomycin')
strep_abx_bile_metabolome$abx <- NULL
rm(abx_bile_metabolome, bile_metabolome)
# Find significant differences
bile_abx_pvalues1 <- c()
for (i in 1:ncol(cef_abx_bile_metabolome)){bile_abx_pvalues1[i] <- wilcox.test(noabx_bile_metabolome[,i], strep_abx_bile_metabolome[,i], exact=FALSE)$p.value}
bile_abx_pvalues2 <- c()
for (i in 1:ncol(cef_abx_bile_metabolome)){bile_abx_pvalues2[i] <- wilcox.test(noabx_bile_metabolome[,i], cef_abx_bile_metabolome[,i], exact=FALSE)$p.value}
bile_abx_pvalues3 <- c()
for (i in 1:ncol(clinda_abx_bile_metabolome)){bile_abx_pvalues3[i] <- wilcox.test(noabx_bile_metabolome[,i], clinda_abx_bile_metabolome[,i], exact=FALSE)$p.value}
bile_abx_pvalues1 <- p.adjust(bile_abx_pvalues1, method='BH')
bile_abx_pvalues2 <- p.adjust(bile_abx_pvalues2, method='BH')
bile_abx_pvalues3 <- p.adjust(bile_abx_pvalues3, method='BH')
bile_abx_pvalues1[is.na(bile_abx_pvalues1)] <- 1
bile_abx_pvalues2[is.na(bile_abx_pvalues2)] <- 1
bile_abx_pvalues3[is.na(bile_abx_pvalues3)] <- 1

carb_metabolome <- clean_merge(metadata, carb_metabolome)
abx_carb_metabolome <- subset(carb_metabolome, abx != 'germfree')
abx_carb_metabolome$infection <- NULL
abx_carb_metabolome$susceptibility <- NULL
abx_carb_metabolome$abx <- factor(abx_carb_metabolome$abx)
carb_metabolome$infection <- NULL
carb_metabolome$abx <- NULL
carb_metabolome$abx <- factor(carb_metabolome$susceptibility)
noabx_carb_metabolome <- subset(abx_carb_metabolome, abx == 'none')
noabx_carb_metabolome$abx <- NULL
cef_abx_carb_metabolome <- subset(abx_carb_metabolome, abx == 'cefoperazone')
cef_abx_carb_metabolome$abx <- NULL
clinda_abx_carb_metabolome <- subset(abx_carb_metabolome, abx == 'clindamycin')
clinda_abx_carb_metabolome$abx <- NULL
strep_abx_carb_metabolome <- subset(abx_carb_metabolome, abx == 'streptomycin')
strep_abx_carb_metabolome$abx <- NULL
rm(abx_carb_metabolome, carb_metabolome)
# Find significant differences
carb_abx_pvalues1 <- c()
for (i in 1:ncol(cef_abx_carb_metabolome)){carb_abx_pvalues1[i] <- wilcox.test(noabx_carb_metabolome[,i], strep_abx_carb_metabolome[,i], exact=FALSE)$p.value}
carb_abx_pvalues2 <- c()
for (i in 1:ncol(cef_abx_carb_metabolome)){carb_abx_pvalues2[i] <- wilcox.test(noabx_carb_metabolome[,i], cef_abx_carb_metabolome[,i], exact=FALSE)$p.value}
carb_abx_pvalues3 <- c()
for (i in 1:ncol(clinda_abx_carb_metabolome)){carb_abx_pvalues3[i] <- wilcox.test(noabx_carb_metabolome[,i], clinda_abx_carb_metabolome[,i], exact=FALSE)$p.value}
carb_abx_pvalues1 <- p.adjust(carb_abx_pvalues1, method='BH')
carb_abx_pvalues2 <- p.adjust(carb_abx_pvalues2, method='BH')
carb_abx_pvalues3 <- p.adjust(carb_abx_pvalues3, method='BH')
carb_abx_pvalues1[is.na(carb_abx_pvalues1)] <- 1
carb_abx_pvalues2[is.na(carb_abx_pvalues2)] <- 1
carb_abx_pvalues3[is.na(carb_abx_pvalues3)] <- 1

amino_metabolome <- clean_merge(metadata, amino_metabolome)
abx_amino_metabolome <- subset(amino_metabolome, abx != 'germfree')
abx_amino_metabolome$infection <- NULL
abx_amino_metabolome$susceptibility <- NULL
abx_amino_metabolome$abx <- factor(abx_amino_metabolome$abx)
amino_metabolome$infection <- NULL
amino_metabolome$abx <- NULL
amino_metabolome$abx <- factor(amino_metabolome$susceptibility)
noabx_amino_metabolome <- subset(abx_amino_metabolome, abx == 'none')
noabx_amino_metabolome$abx <- NULL
cef_abx_amino_metabolome <- subset(abx_amino_metabolome, abx == 'cefoperazone')
cef_abx_amino_metabolome$abx <- NULL
clinda_abx_amino_metabolome <- subset(abx_amino_metabolome, abx == 'clindamycin')
clinda_abx_amino_metabolome$abx <- NULL
strep_abx_amino_metabolome <- subset(abx_amino_metabolome, abx == 'streptomycin')
strep_abx_amino_metabolome$abx <- NULL
rm(abx_amino_metabolome, amino_metabolome)
# Find significant differences
amino_abx_pvalues1 <- c()
for (i in 1:ncol(cef_abx_amino_metabolome)){amino_abx_pvalues1[i] <- wilcox.test(noabx_amino_metabolome[,i], strep_abx_amino_metabolome[,i], exact=FALSE)$p.value}
amino_abx_pvalues2 <- c()
for (i in 1:ncol(cef_abx_amino_metabolome)){amino_abx_pvalues2[i] <- wilcox.test(noabx_amino_metabolome[,i], cef_abx_amino_metabolome[,i], exact=FALSE)$p.value}
amino_abx_pvalues3 <- c()
for (i in 1:ncol(clinda_abx_amino_metabolome)){amino_abx_pvalues3[i] <- wilcox.test(noabx_amino_metabolome[,i], clinda_abx_amino_metabolome[,i], exact=FALSE)$p.value}
amino_abx_pvalues1 <- p.adjust(amino_abx_pvalues1, method='BH')
amino_abx_pvalues2 <- p.adjust(amino_abx_pvalues2, method='BH')
amino_abx_pvalues3 <- p.adjust(amino_abx_pvalues3, method='BH')
amino_abx_pvalues1[is.na(amino_abx_pvalues1)] <- 1
amino_abx_pvalues2[is.na(amino_abx_pvalues2)] <- 1
amino_abx_pvalues3[is.na(amino_abx_pvalues3)] <- 1

#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_a, width=4, height=ncol(strep_abx_bile_metabolome)*2)
layout(matrix(c(1:(ncol(strep_abx_bile_metabolome)+2)), nrow=(ncol(strep_abx_bile_metabolome)+2), ncol=1, byrow = TRUE))

par(mar=c(0.2, 0, 0, 2), mgp=c(2.3, 0.75, 0), xpd=FALSE)
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
text(x=-10.2, y=-3, labels='A', cex=2.4, font=2, xpd=TRUE)
legend('bottomright', legend=c('No Antibiotics','Streptomycin', 'Cefoperazone', 'Clindamycin'), bty='n',
       pt.bg=c(noabx_col,strep_col,cef_col,clinda_col), pch=21, cex=1.2, pt.cex=2, ncol=2)

par(mar=c(0.2, 2, 0.2, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE, yaxs='i')
for(i in c(1:(ncol(strep_abx_bile_metabolome)))){
  xmax <- ceiling(max(c(max(strep_abx_bile_metabolome[,i]), max(cef_abx_bile_metabolome[,i]), 
                        max(clinda_abx_bile_metabolome[,i]), max(noabx_bile_metabolome[,i]))))
  while(xmax %% 5 != 0 ){xmax <- xmax + 1}
  plot(0, type='n', xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,xmax), ylim=c(0.3,2.8))
  stripchart(at=2.1, jitter(noabx_bile_metabolome[,i], amount=1e-5),
             pch=21, bg=noabx_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  stripchart(at=1.65, jitter(strep_abx_bile_metabolome[,i], amount=1e-5),
             pch=21, bg=strep_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  stripchart(at=1.2, jitter(cef_abx_bile_metabolome[,i], amount=1e-5), 
             pch=21, bg=cef_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  stripchart(at=0.66, jitter(clinda_abx_bile_metabolome[,i], amount=1e-5), 
             pch=21, bg=clinda_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  legend('topright', legend=colnames(strep_abx_bile_metabolome)[i], pch=1, cex=1.3, pt.cex=0, bty='n')
  
  if (bile_abx_pvalues1[i] < 0.001){
    text(x=xmax, y=1.65 , '***', cex=1.5, font=2, srt=90)
  } else if (bile_abx_pvalues1[i] <= 0.01){
    text(x=xmax, y=1.65 , '**', cex=1.5, font=2, srt=90)
  } else if (bile_abx_pvalues1[i] <= 0.05){
    text(x=xmax, y=1.65 , '*', cex=1.5, font=2, srt=90)
  } else {
    text(x=xmax, y=1.65 , 'n.s.', cex=1, font=2, srt=90)
  }
  if (bile_abx_pvalues2[i] < 0.001){
    text(x=xmax, y=1.2 , '***', cex=1.5, font=2, srt=90)
  } else if (bile_abx_pvalues2[i] <= 0.01){
    text(x=xmax, y=1.2 , '**', cex=1.5, font=2, srt=90)
  } else if (bile_abx_pvalues2[i] <= 0.05){
    text(x=xmax, y=1.2 , '*', cex=1.5, font=2, srt=90)
  } else {
    text(x=xmax, y=1.2 , 'n.s.', cex=1, font=2, srt=90)
  }
  if (bile_abx_pvalues3[i] < 0.001){
    text(x=xmax, y=0.66 , '***', cex=1.5, font=2, srt=90)
  } else if (bile_abx_pvalues3[i] <= 0.01){
    text(x=xmax, y=0.66 , '**', cex=1.5, font=2, srt=90)
  } else if (bile_abx_pvalues3[i] <= 0.05){
    text(x=xmax, y=0.66 , '*', cex=1.5, font=2, srt=90)
  } else {
    text(x=xmax, y=0.66 , 'n.s.', cex=1, font=2, srt=90)
  }
  if (xmax <= 10) {
    text(x=seq(0,xmax,1), y=0.42, labels=seq(0,xmax,1), cex=1)
    axis(1, at=seq(0,5,1), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 2000){
    text(x=seq(0,xmax,500), y=0.42, labels=seq(0,xmax,500), cex=1)
    axis(1, at=seq(0,xmax,500), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 1000){
    text(x=seq(0,xmax,200), y=0.42, labels=seq(0,xmax,200), cex=1)
    axis(1, at=seq(0,xmax,200), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 500){
    text(x=seq(0,xmax,100), y=0.42, labels=seq(0,xmax,100), cex=1)
    axis(1, at=seq(0,xmax,100), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 100){
    text(x=seq(0,xmax,50), y=0.42, labels=seq(0,xmax,50), cex=1)
    axis(1, at=seq(0,xmax,50), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 50){
    text(x=seq(0,xmax,10), y=0.42, labels=seq(0,xmax,10), cex=1)
    axis(1, at=seq(0,xmax,10), NA, cex.axis=0.8, tck=0.015)
  } else {
    text(x=seq(0,xmax,5), y=0.42, labels=seq(0,xmax,5), cex=1)
    axis(1, at=seq(0,xmax,5), NA, cex.axis=0.8, tck=0.015)
  }
  
  segments(median(noabx_bile_metabolome[,i]), 1.93, median(noabx_bile_metabolome[,i]), 2.27, lwd=2.5)
  segments(median(strep_abx_bile_metabolome[,i]), 1.48, median(strep_abx_bile_metabolome[,i]), 1.82, lwd=2.5)
  segments(median(cef_abx_bile_metabolome[,i]), 1.03, median(cef_abx_bile_metabolome[,i]), 1.37, lwd=2.5)
  segments(median(clinda_abx_bile_metabolome[,i]), 0.49, median(clinda_abx_bile_metabolome[,i]), 0.83, lwd=2.5)
}
par(mar=c(0, 0, 0, 0))
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
text(x=0, y=4, labels=expression(paste('Scaled Intensity (',log[10],')')), cex=1.4)

dev.off()


pdf(file=plot_b, width=4, height=ncol(strep_abx_carb_metabolome)*2)
layout(matrix(c(1:(ncol(strep_abx_carb_metabolome)+2)), nrow=(ncol(strep_abx_carb_metabolome)+2), ncol=1, byrow = TRUE))

par(mar=c(0.2, 0, 0, 2), mgp=c(2.3, 0.75, 0), xpd=FALSE)
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
text(x=-10.2, y=-3, labels='B', cex=2.4, font=2, xpd=TRUE)
legend('bottomright', legend=c('No Antibiotics','Streptomycin', 'Cefoperazone', 'Clindamycin'), bty='n',
       pt.bg=c(noabx_col,strep_col,cef_col,clinda_col), pch=21, cex=1.2, pt.cex=2, ncol=2)

par(mar=c(0.2, 2, 0.2, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE, yaxs='i')
for(i in c(1:(ncol(strep_abx_carb_metabolome)))){
  xmax <- ceiling(max(c(max(strep_abx_carb_metabolome[,i]), max(cef_abx_carb_metabolome[,i]), 
                        max(clinda_abx_carb_metabolome[,i]), max(noabx_carb_metabolome[,i]))))
  while(xmax %% 5 != 0 ){xmax <- xmax + 1}
  plot(0, type='n', xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,xmax), ylim=c(0.3,2.8))
  stripchart(at=2.1, jitter(noabx_carb_metabolome[,i], amount=1e-5),
             pch=21, bg=noabx_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  stripchart(at=1.65, jitter(strep_abx_carb_metabolome[,i], amount=1e-5),
             pch=21, bg=strep_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  stripchart(at=1.2, jitter(cef_abx_carb_metabolome[,i], amount=1e-5), 
             pch=21, bg=cef_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  stripchart(at=0.66, jitter(clinda_abx_carb_metabolome[,i], amount=1e-5), 
             pch=21, bg=clinda_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  legend('topright', legend=colnames(strep_abx_carb_metabolome)[i], pch=1, cex=1.3, pt.cex=0, bty='n')
  
  if (carb_abx_pvalues1[i] < 0.001){
    text(x=xmax, y=1.65 , '***', cex=1.5, font=2, srt=90)
  } else if (carb_abx_pvalues1[i] <= 0.01){
    text(x=xmax, y=1.65 , '**', cex=1.5, font=2, srt=90)
  } else if (carb_abx_pvalues1[i] <= 0.05){
    text(x=xmax, y=1.65 , '*', cex=1.5, font=2, srt=90)
  } else {
    text(x=xmax, y=1.65 , 'n.s.', cex=1, font=2, srt=90)
  }
  if (carb_abx_pvalues2[i] < 0.001){
    text(x=xmax, y=1.2 , '***', cex=1.5, font=2, srt=90)
  } else if (carb_abx_pvalues2[i] <= 0.01){
    text(x=xmax, y=1.2 , '**', cex=1.5, font=2, srt=90)
  } else if (carb_abx_pvalues2[i] <= 0.05){
    text(x=xmax, y=1.2 , '*', cex=1.5, font=2, srt=90)
  } else {
    text(x=xmax, y=1.2 , 'n.s.', cex=1, font=2, srt=90)
  }
  if (carb_abx_pvalues3[i] < 0.001){
    text(x=xmax, y=0.66 , '***', cex=1.5, font=2, srt=90)
  } else if (carb_abx_pvalues3[i] <= 0.01){
    text(x=xmax, y=0.66 , '**', cex=1.5, font=2, srt=90)
  } else if (carb_abx_pvalues3[i] <= 0.05){
    text(x=xmax, y=0.66 , '*', cex=1.5, font=2, srt=90)
  } else {
    text(x=xmax, y=0.66 , 'n.s.', cex=1, font=2, srt=90)
  }
  if (xmax <= 10) {
    text(x=seq(0,xmax,1), y=0.42, labels=seq(0,xmax,1), cex=1)
    axis(1, at=seq(0,5,1), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 2000){
    text(x=seq(0,xmax,500), y=0.42, labels=seq(0,xmax,500), cex=1)
    axis(1, at=seq(0,xmax,500), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 1000){
    text(x=seq(0,xmax,200), y=0.42, labels=seq(0,xmax,200), cex=1)
    axis(1, at=seq(0,xmax,200), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 500){
    text(x=seq(0,xmax,100), y=0.42, labels=seq(0,xmax,100), cex=1)
    axis(1, at=seq(0,xmax,100), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 100){
    text(x=seq(0,xmax,50), y=0.42, labels=seq(0,xmax,50), cex=1)
    axis(1, at=seq(0,xmax,50), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 50){
    text(x=seq(0,xmax,10), y=0.42, labels=seq(0,xmax,10), cex=1)
    axis(1, at=seq(0,xmax,10), NA, cex.axis=0.8, tck=0.015)
  } else {
    text(x=seq(0,xmax,5), y=0.42, labels=seq(0,xmax,5), cex=1)
    axis(1, at=seq(0,xmax,5), NA, cex.axis=0.8, tck=0.015)
  }
  
  segments(median(noabx_carb_metabolome[,i]), 1.93, median(noabx_carb_metabolome[,i]), 2.27, lwd=2.5)
  segments(median(strep_abx_carb_metabolome[,i]), 1.48, median(strep_abx_carb_metabolome[,i]), 1.82, lwd=2.5)
  segments(median(cef_abx_carb_metabolome[,i]), 1.03, median(cef_abx_carb_metabolome[,i]), 1.37, lwd=2.5)
  segments(median(clinda_abx_carb_metabolome[,i]), 0.49, median(clinda_abx_carb_metabolome[,i]), 0.83, lwd=2.5)
}
par(mar=c(0, 0, 0, 0))
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
text(x=0, y=4, labels=expression(paste('Scaled Intensity (',log[10],')')), cex=1.4)

dev.off()



pdf(file=plot_c, width=4, height=ncol(strep_abx_amino_metabolome)*2)
layout(matrix(c(1:(ncol(strep_abx_amino_metabolome)+2)), nrow=(ncol(strep_abx_amino_metabolome)+2), ncol=1, byrow = TRUE))

par(mar=c(0.2, 0, 0, 2), mgp=c(2.3, 0.75, 0), xpd=FALSE)
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
text(x=-10.2, y=-3, labels='C', cex=2.4, font=2, xpd=TRUE)
legend('bottomright', legend=c('No Antibiotics','Streptomycin', 'Cefoperazone', 'Clindamycin'), bty='n',
       pt.bg=c(noabx_col,strep_col,cef_col,clinda_col), pch=21, cex=1.2, pt.cex=2, ncol=2)

par(mar=c(0.2, 2, 0.2, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE, yaxs='i')
for(i in c(1:(ncol(strep_abx_amino_metabolome)))){
  xmax <- ceiling(max(c(max(strep_abx_amino_metabolome[,i]), max(cef_abx_amino_metabolome[,i]), 
                        max(clinda_abx_amino_metabolome[,i]), max(noabx_amino_metabolome[,i]))))
  while(xmax %% 5 != 0 ){xmax <- xmax + 1}
  plot(0, type='n', xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,xmax), ylim=c(0.3,2.8))
  stripchart(at=2.1, jitter(noabx_amino_metabolome[,i], amount=1e-5),
             pch=21, bg=noabx_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  stripchart(at=1.65, jitter(strep_abx_amino_metabolome[,i], amount=1e-5),
             pch=21, bg=strep_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  stripchart(at=1.2, jitter(cef_abx_amino_metabolome[,i], amount=1e-5), 
             pch=21, bg=cef_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  stripchart(at=0.66, jitter(clinda_abx_amino_metabolome[,i], amount=1e-5), 
             pch=21, bg=clinda_col, method='jitter', jitter=0.12, cex=2, lwd=0.5, add=TRUE)
  legend('topright', legend=colnames(strep_abx_amino_metabolome)[i], pch=1, cex=1.3, pt.cex=0, bty='n')
  
  if (amino_abx_pvalues1[i] < 0.001){
    text(x=xmax, y=1.65 , '***', cex=1.5, font=2, srt=90)
  } else if (amino_abx_pvalues1[i] <= 0.01){
    text(x=xmax, y=1.65 , '**', cex=1.5, font=2, srt=90)
  } else if (amino_abx_pvalues1[i] <= 0.05){
    text(x=xmax, y=1.65 , '*', cex=1.5, font=2, srt=90)
  } else {
    text(x=xmax, y=1.65 , 'n.s.', cex=1, font=2, srt=90)
  }
  if (amino_abx_pvalues2[i] < 0.001){
    text(x=xmax, y=1.2 , '***', cex=1.5, font=2, srt=90)
  } else if (amino_abx_pvalues2[i] <= 0.01){
    text(x=xmax, y=1.2 , '**', cex=1.5, font=2, srt=90)
  } else if (amino_abx_pvalues2[i] <= 0.05){
    text(x=xmax, y=1.2 , '*', cex=1.5, font=2, srt=90)
  } else {
    text(x=xmax, y=1.2 , 'n.s.', cex=1, font=2, srt=90)
  }
  if (amino_abx_pvalues3[i] < 0.001){
    text(x=xmax, y=0.66 , '***', cex=1.5, font=2, srt=90)
  } else if (amino_abx_pvalues3[i] <= 0.01){
    text(x=xmax, y=0.66 , '**', cex=1.5, font=2, srt=90)
  } else if (amino_abx_pvalues3[i] <= 0.05){
    text(x=xmax, y=0.66 , '*', cex=1.5, font=2, srt=90)
  } else {
    text(x=xmax, y=0.66 , 'n.s.', cex=1, font=2, srt=90)
  }
  if (xmax <= 10) {
    text(x=seq(0,xmax,1), y=0.42, labels=seq(0,xmax,1), cex=1)
    axis(1, at=seq(0,5,1), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 2000){
    text(x=seq(0,xmax,500), y=0.42, labels=seq(0,xmax,500), cex=1)
    axis(1, at=seq(0,xmax,500), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 1000){
    text(x=seq(0,xmax,200), y=0.42, labels=seq(0,xmax,200), cex=1)
    axis(1, at=seq(0,xmax,200), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 500){
    text(x=seq(0,xmax,100), y=0.42, labels=seq(0,xmax,100), cex=1)
    axis(1, at=seq(0,xmax,100), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 100){
    text(x=seq(0,xmax,50), y=0.42, labels=seq(0,xmax,50), cex=1)
    axis(1, at=seq(0,xmax,50), NA, cex.axis=0.8, tck=0.015)
  } else if (xmax > 50){
    text(x=seq(0,xmax,10), y=0.42, labels=seq(0,xmax,10), cex=1)
    axis(1, at=seq(0,xmax,10), NA, cex.axis=0.8, tck=0.015)
  } else {
    text(x=seq(0,xmax,5), y=0.42, labels=seq(0,xmax,5), cex=1)
    axis(1, at=seq(0,xmax,5), NA, cex.axis=0.8, tck=0.015)
  }
  
  segments(median(noabx_amino_metabolome[,i]), 1.93, median(noabx_amino_metabolome[,i]), 2.27, lwd=2.5)
  segments(median(strep_abx_amino_metabolome[,i]), 1.48, median(strep_abx_amino_metabolome[,i]), 1.82, lwd=2.5)
  segments(median(cef_abx_amino_metabolome[,i]), 1.03, median(cef_abx_amino_metabolome[,i]), 1.37, lwd=2.5)
  segments(median(clinda_abx_amino_metabolome[,i]), 0.49, median(clinda_abx_amino_metabolome[,i]), 0.83, lwd=2.5)
}
par(mar=c(0, 0, 0, 0))
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
text(x=0, y=4, labels=expression(paste('Scaled Intensity (',log[10],')')), cex=1.4)

dev.off()


#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()
