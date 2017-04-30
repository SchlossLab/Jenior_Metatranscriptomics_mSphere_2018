
# Start with clean environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define input files
metabolome <- 'data/wetlab_assays/metabolomics.tsv'
cef_importances <- 'data/metabolic_models/cefoperazone_630.bipartite.files/importances.tsv'
strep_importances <- 'data/metabolic_models/streptomycin_630.bipartite.files/importances.tsv'
clinda_importances <- 'data/metabolic_models/clindamycin_630.bipartite.files/importances.tsv'
metadata <- 'data/metadata.tsv'

# Define output files
supptable_S5A_file <'results/supplement/tables/table_S5strep.tsv'
supptable_S5B_file <'results/supplement/tables/table_S5cef.tsv'
supptable_S5C_file <'results/supplement/tables/table_S5clinda.tsv'
supptable_S5D_file <'results/supplement/tables/table_S5all.tsv'
plot_file <- 'results/figures/figure_6.pdf'


# Metabolomic data
metabolome <- read.delim(metabolome, sep='\t', header=T)

# Remove germfree samples

annotation <- metabolome[,1:5]
metabolome$SUPER_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$PUBCHEM <- NULL
metabolome <- subset(metabolome, KEGG != 'NA')
metabolome <- metabolome[match(unique(metabolome$KEGG), metabolome$KEGG),]
rownames(metabolome) <- metabolome$KEGG
metabolome$KEGG <- NULL
metabolome$BIOCHEMICAL <- NULL
annotation <- subset(annotation, KEGG != 'NA')
annotation <- annotation[match(unique(annotation$KEGG), annotation$KEGG),]
rownames(annotation) <- annotation$KEGG
annotation$KEGG <- NULL
annotation$SUB_PATHWAY <- NULL
annotation$PUBCHEM <- NULL
annotation$BIOCHEMICAL <- NULL
annotation$SUPER_PATHWAY <- gsub('_',' ', annotation$SUPER_PATHWAY)

# Combine with metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)
metabolome <- merge(t(metabolome), metadata, by='row.names')
rm(metadata)
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
metabolome <- subset(metabolome, abx != 'none')
metabolome_630 <- subset(metabolome, infection == '630')
metabolome_630$cage <- NULL
metabolome_630$mouse <- NULL
metabolome_630$gender <- NULL
metabolome_630$type <- NULL
metabolome_630$infection <- NULL
metabolome_630 <- aggregate(metabolome_630[, 1:399], list(metabolome_630$abx), median)
rownames(metabolome_630) <- metabolome_630$Group.1
metabolome_630$Group.1 <- NULL
metabolome_630 <- as.data.frame(t(metabolome_630))
colnames(metabolome_630) <- c('cef_630','clinda_630','strep_630')
metabolome_mock <- subset(metabolome, infection == 'mock')
metabolome_mock$cage <- NULL
metabolome_mock$mouse <- NULL
metabolome_mock$gender <- NULL
metabolome_mock$type <- NULL
metabolome_mock$infection <- NULL
metabolome_mock <- aggregate(metabolome_mock[, 1:399], list(metabolome_mock$abx), median)
rownames(metabolome_mock) <- metabolome_mock$Group.1
metabolome_mock$Group.1 <- NULL
metabolome_mock <- as.data.frame(t(metabolome_mock))
colnames(metabolome_mock) <- c('cef_mock','clinda_mock','strep_mock')

# Calculate ratio of mock to infected concentrations
metabolome_mock <- 10 ^ metabolome_mock
metabolome_630 <- 10 ^ metabolome_630
metabolome <- metabolome_mock / metabolome_630
metabolome <- log10(metabolome)
metabolome <- merge(metabolome, annotation, by='row.names')
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
colnames(metabolome) <- c('cefoperazone_conc', 'clindamycin_conc', 'streptomycin_conc', 'pathway')
rm(annotation, metabolome_mock, metabolome_630)

# Merge importances to one table
cef_importances <- read.delim(cef_importances, sep='\t', header=T, row.names=1)
cef_importances$p_value <- NULL
clinda_importances <- read.delim(clinda_importances, sep='\t', header=T, row.names=1)
clinda_importances$p_value <- NULL
clinda_importances$Compound_name <- NULL
strep_importances <- read.delim(strep_importances, sep='\t', header=T, row.names=1)
strep_importances$p_value <- NULL
strep_importances$Compound_name <- NULL
importances <- merge(cef_importances, clinda_importances, by='row.names')
rownames(importances) <- importances$Row.names
importances$Row.names <- NULL
importances <- merge(importances, strep_importances, by='row.names')
rownames(importances) <- importances$Row.names
importances$Row.names <- NULL
colnames(importances) <- c('Compound_name', 'cefoperazone_score', 'clindamycin_score', 'streptomycin_score')
importances$Compound_name <- gsub('_',' ', importances$Compound_name)
rm(cef_importances, clinda_importances, strep_importances)

# Merge metabolome medians and importance values
combined <- merge(importances, metabolome, by='row.names')
rownames(combined) <- combined$Row.names
combined$Row.names <- NULL
rm(importances, metabolome)

# Separate treatment groups
cef <- as.data.frame(cbind(combined$cefoperazone_score, combined$cefoperazone_conc, combined$pathway))
rownames(cef) <- rownames(combined)
colnames(cef) <- c('score', 'conc', 'pathway')
cef$name <- combined$Compound_name
cef$score <- as.numeric(as.vector(cef$score))
cef$conc <- as.numeric(as.vector(cef$conc))
cef <- subset(cef, cef[,1] != 0) # remove metabolites with 0 importance
clinda <- as.data.frame(cbind(combined$clindamycin_score, combined$clindamycin_conc, combined$pathway))
rownames(clinda) <- rownames(combined)
colnames(clinda) <- c('score', 'conc', 'pathway')
clinda$name <- combined$Compound_name
clinda$score <- as.numeric(as.vector(clinda$score))
clinda$conc <- as.numeric(as.vector(clinda$conc))
clinda <- subset(clinda, clinda[,1] != 0) # remove metabolites with 0 importance
strep <- as.data.frame(cbind(combined$streptomycin_score, combined$streptomycin_conc, combined$pathway))
rownames(strep) <- rownames(combined)
colnames(strep) <- c('score', 'conc', 'pathway')
strep$name <- combined$Compound_name
strep$score <- as.numeric(as.vector(strep$score))
strep$conc <- as.numeric(as.vector(strep$conc))
strep <- subset(strep, strep[,1] != 0) # remove metabolites with 0 importance
combined <- rbind(cef, clinda, strep)

# Fit to general linear models and identify outliers (L1 regression)
strep_fit <- glm(conc ~ score, data=strep)
strep$residuals <- residuals(strep_fit)
strep$residuals <- (strep$residuals / sd(strep$residuals))^2
strep_outliers <- strep[strep$residuals > 1.5, ]
cef_fit <- glm(conc ~ score, data=cef)
cef$residuals <- residuals(cef_fit)
cef$residuals <- (cef$residuals / sd(cef$residuals))^2
cef_outliers <- cef[cef$residuals > 1.5, ]
clinda_fit <- glm(conc ~ score, data=clinda)
clinda$residuals <- residuals(clinda_fit)
clinda$residuals <- (clinda$residuals / sd(clinda$residuals))^2
clinda_outliers <- clinda[clinda$residuals > 1.5, ]
combined_fit <- glm(conc ~ score, data=combined)
combined$residuals <- residuals(combined_fit)
combined$residuals <- (combined$residuals / sd(combined$residuals))^2
combined_outliers <- combined[combined$residuals > 1.5, ]

# Calculate stats
test <- as.data.frame(cbind(round(c(strep_fit$coefficients[[2]],
                                    cef_fit$coefficients[[2]],
                                    clinda_fit$coefficients[[2]],
                                    combined_fit$coefficients[[2]]), digits=3),
                            round(c(cor.test(strep[,1], strep[,2], method='spearman', exact=FALSE)$estimate,
                                    cor.test(cef[,1], cef[,2], method='spearman', exact=FALSE)$estimate,
                                    cor.test(clinda[,1], clinda[,2], method='spearman', exact=FALSE)$estimate,
                                    cor.test(combined[,1], combined[,2], method='spearman', exact=FALSE)$estimate), digits=3),
                            round(c(cor.test(strep[,1], strep[,2], method='spearman', exact=FALSE)$p.value,
                                             cor.test(cef[,1], cef[,2], method='spearman', exact=FALSE)$p.value,
                                             cor.test(clinda[,1], clinda[,2], method='spearman', exact=FALSE)$p.value,
                                             cor.test(combined[,1], combined[,2], method='spearman', exact=FALSE)$p.value), digits=3)))
rownames(test) <- c('streptomycin','cefoperazone','clindamycin','combined')
colnames(test) <- c('m','r', 'p')

# Write outliers to supplementary table
write.table(strep_outliers, file=supptable_S5A_file, quote=FALSE, sep='\t', row.names=TRUE)
write.table(cef_outliers, file=supptable_S5B_file, quote=FALSE, sep='\t', row.names=TRUE)
write.table(clinda_outliers, file=supptable_S5C_file, quote=FALSE, sep='\t', row.names=TRUE)
write.table(combined_outliers, file=supptable_S5D_file, quote=FALSE, sep='\t', row.names=TRUE)
# Assemble into multi-paneled Excel table downstream

# Add column for colors to combined outliers
index <- as.vector(unique(combined_outliers$pathway))
values <- c('chartreuse1', 'darkorchid2', 'gold')
combined_outliers$color <- values[match(combined_outliers$pathway, index)]

#----------------------------------------#

# Set up multi-panel figure
pdf(file=plot_file, width=9, height=9)
layout(matrix(c(1,2,
                3,4), 
              nrow=2, ncol=2, byrow=TRUE))
par(las=1, mar=c(3,3,1,1), mgp=c(1.8,0.7,0))

# Plot the data and correlations
plot(combined[,1], combined[,2], xlab='Importance Score', ylab=expression(paste(Delta,' Median Scaled Intensity')), 
     pch=19, cex=1.1, xlim=c(-10,10), ylim=c(-1,15), col='gray20')
filledrectangle(wx=24, wy=3.4, col='gray90', mid=c(0,1.2), angle=6)
abline(v=0, lty=2, col='gray60')
box()
points(combined[,1], combined[,2], pch=19, col='gray20', cex=1.1)
abline(combined_fit, col='black', lwd=2)
mtext('A', side=2, line=2, las=2, adj=0.8, padj=-8, cex=1.2)
points(combined_outliers[,1], combined_outliers[,2], pch=21, col='gray20', bg=combined_outliers$color, cex=2, lwd=2)
legend('topleft', legend=as.vector(unique(combined_outliers$pathway)), 
       pt.bg=c('chartreuse1', 'darkorchid2', 'gold'), col='gray25',
       pch=21, pt.lwd=2, pt.cex=2, cex=1.2)
legend('topright', legend=c('All Treatment Groups'), pt.cex=0, bty='n', cex=1.2, col='black')
text(x=-8.75, y=9.3, as.expression(bquote(atop(paste(italic('rho'),' = ',.(test$r[5])),
                                          paste(italic('P'),' = ',.(test$p[5]),'*  ')))), cex=1.2)

# streptomycin alone
plot(strep[,1], strep[,2], xlab='Importance Score', ylab=expression(paste(Delta,' Median Scaled Intensity')), 
     pch=19, xlim=c(-10,10), ylim=c(-1,13), col=wes_palette("FantasticFox")[1])
filledrectangle(wx=24, wy=3.4, col='gray90', mid=c(0,1), angle=6)
abline(v=0, lty=2, col='gray60')
box()
points(strep[,1], strep[,2], pch=19, col=wes_palette("FantasticFox")[1])
abline(strep_fit, col='black', lwd=2)
mtext('B', side=2, line=2, las=2, adj=0.8, padj=-8, cex=1.2)
legend('topright', legend='Streptomycin', bty='n', cex=1.1, col='black')
legend('topleft', legend=c(as.expression(bquote(paste(italic('rho'),' = ',.(test$r[1])))), 
                           as.expression(bquote(paste(italic('P'),' = ',.(test$p[1]))))), pt.cex=0, bty='n', cex=1.1)
points(strep_outliers[,1], strep_outliers[,2], pch=21, bg=wes_palette("FantasticFox")[1], cex=1.7, lwd=2)
text(x=c(5.5,3.7,5.6), 
     y=c(10.6,7.875043,3.8), 
     strep_outliers$name, cex=0.9, col='gray20')

# cefoperazone alone
plot(cef[,1], cef[,2], xlab='Importance Score', ylab=expression(paste(Delta,' Median Scaled Intensity')), 
     pch=19, xlim=c(-10,10), ylim=c(-1,11), col=wes_palette("FantasticFox")[3]) 
filledrectangle(wx=24, wy=2.4, col='gray90', mid=c(0,1.4), angle=3)
abline(v=0, lty=2, col='gray60')
box()
points(cef[,1], cef[,2], pch=19, col=wes_palette("FantasticFox")[3])
abline(cef_fit, col='black', lwd=2)
mtext('C', side=2, line=2, las=2, adj=0.8, padj=-8, cex=1.2)
legend('topleft', legend=c(as.expression(bquote(paste(italic('rho'),' = ',.(test$r[2])))), 
                           as.expression(bquote(paste(italic('P'),' = ',.(test$p[2]))))), pt.cex=0, bty='n', cex=1.1)
legend('topright', legend='Cefoperazone', bty='n', cex=1.1, col='black')
points(cef_outliers[,1], cef_outliers[,2], pch=21, bg=wes_palette("FantasticFox")[3], cex=1.7, lwd=2)
text(x=c(8,5,-6,7,5.7,5,-5.5), 
     y=c(2.8829322,6,2.6819865,7.9,4.7,-0.3,3.7199615), 
     cef_outliers$name, cex=0.9, col='gray20')
segments(x0=c(2.8,0.9,5.1), y0=c(2.9,3.5,3.8), 
         x1=c(5.3,3,5.6), y1=c(2.9,5.7,4.4), col='gray20')

# clindamycin alone
plot(clinda[,1], clinda[,2], xlab='Importance Score', ylab=expression(paste(Delta,' Median Scaled Intensity')), 
     pch=19, xlim=c(-10,10), ylim=c(-1,5), col=wes_palette("FantasticFox")[5])
filledrectangle(wx=24, wy=1.3, col='gray90', mid=c(0,1), angle=3)
abline(v=0, lty=2, col='gray60')
box()
points(clinda[,1], clinda[,2], pch=19, col=wes_palette("FantasticFox")[5])
abline(clinda_fit, col='black', lwd=2)
mtext('D', side=2, line=2, las=2, adj=0.8, padj=-8, cex=1.2)
legend('topleft', legend=c(as.expression(bquote(paste(italic('rho'),' = ',.(test$r[3])))), 
                           as.expression(bquote(paste(italic('P'),' = ',.(test$p[3]),'*')))), pt.cex=0, bty='n', cex=1.1)
legend('topright', legend='Clindamycin', bty='n', cex=1.1, col='black')
points(clinda_outliers[,1], clinda_outliers[,2], pch=21, bg=wes_palette("FantasticFox")[5], cex=1.7, lwd=2)
text(x=c(6.4,-3.8,4.9,7.6,6.3,-4,-7,5.8,5.2,6.1), 
     y=c(-0.1,1.8708661,2.2279556,4,1.85,-1,2.6433240,-0.6,3.2,0.4666667), 
     clinda_outliers$name, cex=0.9, col='gray20')
segments(x0=c(-4.3,1.1,3.2), y0=c(-0.8,0,0.3), 
         x1=c(-3.3,1.6,4.2), y1=c(-0.05,-0.4,0), col='gray20')


dev.off()

#----------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()
