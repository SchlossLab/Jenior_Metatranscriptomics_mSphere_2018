
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Output plot name
plot_file <- 'results/figures/figure_5.pdf'

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
metabolome <- metabolome[,!colnames(metabolome) %in% c('GfC1M1','GfC1M2','GfC1M3', 
                                                       'GfC2M1','GfC2M2','GfC2M3', 
                                                       'GfC3M1','GfC3M2','GfC3M3', 
                                                       'GfC4M1','GfC4M2','GfC4M3', 
                                                       'GfC5M1','GfC5M2','GfC5M3', 
                                                       'GfC6M1','GfC6M2','GfC6M3')] # Germfree
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.data.frame(t(metabolome))

#-------------------------------------------------------------------------------------------------------------------------#

# Separate groups
metabolome <- clean_merge(metadata, metabolome)
abx_metabolome <- subset(metabolome, abx != 'none')
abx_metabolome$infection <- NULL
abx_metabolome$susceptibility <- NULL
abx_metabolome$abx <- factor(abx_metabolome$abx)
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
metabolome$infection <- NULL
metabolome$abx <- NULL
metabolome$abx <- factor(metabolome$susceptibility)
rm(metadata)


# Random Forest feature selection
all_rf <- featureselect_RF(metabolome, 'susceptibility')
abx_rf <- featureselect_RF(abx_metabolome, 'abx')
cef_rf <- featureselect_RF(cef_metabolome, 'infection')
clinda_rf <- featureselect_RF(clinda_metabolome, 'infection')
strep_rf <- featureselect_RF(strep_metabolome, 'infection')
all_rf$feature <- gsub('_', ' ', all_rf$feature)
abx_rf$feature <- gsub('_', ' ', abx_rf$feature)
cef_rf$feature <- gsub('_', ' ', cef_rf$feature)
clinda_rf$feature <- gsub('_', ' ', clinda_rf$feature)
strep_rf$feature <- gsub('_', ' ', strep_rf$feature)


# Sort and subset top hits
all_rf <- all_rf[order(-all_rf$final_features_RF),][1:20,]
all_rf <- all_rf[order(all_rf$final_features_RF),]
abx_rf <- abx_rf[order(-abx_rf$final_features_RF),][1:20,]
abx_rf <- abx_rf[order(abx_rf$final_features_RF),]
cef_rf <- cef_rf[order(-cef_rf$final_features_RF),][1:10,]
cef_rf <- cef_rf[order(cef_rf$final_features_RF),]
clinda_rf <- clinda_rf[order(-clinda_rf$final_features_RF),][1:10,]
clinda_rf <- clinda_rf[order(clinda_rf$final_features_RF),]
strep_rf <- strep_rf[order(-strep_rf$final_features_RF),][1:10,]
strep_rf <- strep_rf[order(strep_rf$final_features_RF),]


# Subset concentrations
metabolome <- metabolome[, c(1, which(colnames(metabolome) %in% all_rf$feature))]
res_metabolome <- subset(metabolome, susceptibility == 'resistant')
res_metabolome$susceptibility <- NULL
sus_metabolome <- subset(metabolome, susceptibility == 'susceptible')
sus_metabolome$susceptibility <- NULL

abx_metabolome <- abx_metabolome[, c(1, which(colnames(abx_metabolome) %in% abx_rf$feature))]
cef_abx_metabolome <- subset(abx_metabolome, abx == 'cefoperazone')
cef_abx_metabolome$abx <- NULL
clinda_abx_metabolome <- subset(abx_metabolome, abx == 'clindamycin')
clinda_abx_metabolome$abx <- NULL
strep_abx_metabolome <- subset(abx_metabolome, abx == 'streptomycin')
strep_abx_metabolome$abx <- NULL

cef_metabolome <- cef_metabolome[, c(1, which(colnames(cef_metabolome) %in% cef_rf$feature))]
inf_cef_metabolome <- subset(cef_metabolome, infection == '630')
inf_cef_metabolome$infection<- NULL
mock_cef_metabolome <- subset(cef_metabolome, infection == 'mock')
mock_cef_metabolome$infection<- NULL

clinda_metabolome <- clinda_metabolome[, c(1, which(colnames(clinda_metabolome) %in% clinda_rf$feature))]
inf_clinda_metabolome <- subset(clinda_metabolome, infection == '630')
inf_clinda_metabolome$infection<- NULL
mock_clinda_metabolome <- subset(clinda_metabolome, infection == 'mock')
mock_clinda_metabolome$infection<- NULL

strep_metabolome <- strep_metabolome[, c(1, which(colnames(strep_metabolome) %in% strep_rf$feature))]
inf_strep_metabolome <- subset(strep_metabolome, infection == '630')
inf_strep_metabolome$infection<- NULL
mock_strep_metabolome <- subset(strep_metabolome, infection == 'mock')
mock_strep_metabolome$infection<- NULL

rm(all_rf,abx_rf,cef_rf,clinda_rf,strep_rf,
   metabolome,abx_metabolome,clinda_metabolome,strep_metabolome)


# Find significant differences
resistant_pvalues <- c()
for (i in 1:ncol(res_metabolome)){resistant_pvalues[i] <- wilcox.test(res_metabolome[,i], sus_metabolome[,i], exact=FALSE)$p.value}
resistant_pvalues <- p.adjust(pvalues, method='BH')

abx_pvalues1 <- c()
for (i in 1:ncol(cef_abx_metabolome)){abx_pvalues1[i] <- wilcox.test(cef_abx_metabolome[,i], clinda_abx_metabolome[,i], exact=FALSE)$p.value}
abx_pvalues2 <- c()
for (i in 1:ncol(cef_abx_metabolome)){abx_pvalues2[i] <- wilcox.test(cef_abx_metabolome[,i], strep_abx_metabolome[,i], exact=FALSE)$p.value}
abx_pvalues3 <- c()
for (i in 1:ncol(clinda_abx_metabolome)){abx_pvalues3[i] <- wilcox.test(clinda_abx_metabolome[,i], strep_abx_metabolome[,i], exact=FALSE)$p.value}
abx_pvalues <- c(abx_pvalues1, abx_pvalues2, abx_pvalues3)
abx_pvalues <- p.adjust(abx_pvalues, method='BH')
abx_pvalues1 <- abx_pvalues[1:17]
abx_pvalues2 <- abx_pvalues[18:34]
abx_pvalues3 <- abx_pvalues[35:51]
rm(abx_pvalues)

cef_pvalues <- c()
for (i in 1:ncol(inf_cef_metabolome)){cef_pvalues[i] <- wilcox.test(inf_cef_metabolome[,i], mock_cef_metabolome[,i], exact=FALSE)$p.value}
cef_pvalues <- p.adjust(cef_pvalues, method='BH')
clinda_pvalues <- c()
for (i in 1:ncol(inf_clinda_metabolome)){clinda_pvalues[i] <- wilcox.test(inf_clinda_metabolome[,i], mock_clinda_metabolome[,i], exact=FALSE)$p.value}
pvalues <- p.adjust(clinda_pvalues, method='BH')
strep_pvalues <- c()
for (i in 1:ncol(inf_strep_metabolome)){strep_pvalues[i] <- wilcox.test(inf_strep_metabolome[,i], mock_strep_metabolome[,i], exact=FALSE)$p.value}
strep_pvalues <- p.adjust(strep_pvalues, method='BH')


#-------------------------------------------------------------------------------------------------------------------------#


# Set up multi-panel figure
pdf(file=plot_file, width=12, height=10)
layout(matrix(c(1,2,2,
                3,3,4,
                5,6,7),
              nrow=3, ncol=3, byrow = TRUE))



# Streptomycin plot
par(mar=c(4, 13, 2, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE, yaxs='i')
plot(1, type='n', ylim=c(0,nrow(strep_lefse)*2), xlim=c(0,log10(strep_size)), 
     ylab='', xlab='Relative Abundance %', xaxt='n', yaxt='n')
title('Streptomycin-pretreated', line=0.5, cex.main=1.3, col.main=strep_col, font.main=1)
index <- 1
for(i in c(1:ncol(strep_mock_otu))){
  stripchart(at=index-0.35, jitter(strep_mock_otu[,i], amount=1e-5), 
             pch=21, bg='mediumseagreen', method='jitter', jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  stripchart(at=index+0.35, jitter(strep_infected_otu[,i], amount=1e-5), 
             pch=21, bg='mediumorchid3', method='jitter', jitter=0.12, cex=1.5, lwd=0.5, add=TRUE)
  if (i != ncol(strep_mock_otu)){
    abline(h=index+1, lty=2)
  }
  segments(median(strep_mock_otu[,i]), index-0.6, median(strep_mock_otu[,i]), index, lwd=2) #adds line for median
  segments(median(strep_infected_otu[,i]), index+0.6, median(strep_infected_otu[,i]), index, lwd=2)
  if (wilcox.test(strep_mock_otu[,i], strep_infected_otu[,i], exact=FALSE)$p.value <= 0.05){
    text(x=log10(strep_size)+0.25, y=index, labels='*', font=2, cex=1.8, xpd=TRUE)
  }
  index <- index + 2
}
axis(1, at=c(0,(log10(strep_size/1000)),(log10(strep_size/100)),(log10(strep_size/10)),log10(strep_size)), labels=c('0','0.1','1','10','100')) 
legend('topright', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))), 'Mock-infected'),
       pch=c(21, 21), pt.bg=c('mediumorchid3','mediumseagreen'), bg='white', pt.cex=1.4, cex=0.9)
formatted_names <- lapply(1:nrow(strep_lefse), function(i) bquote(paste(italic(.(strep_lefse$genus[i])), ' ', .(as.vector(strep_lefse$OTU)[i]), sep='')))
axis(2, at=seq(1,index-2,2)+0.4, labels=do.call(expression, formatted_names), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 
formatted_p <- lapply(1:nrow(strep_lefse), function(i) bquote(paste(.(as.vector(strep_lefse$phylum)[i]), '; ', italic('p'), ' = ', .(strep_lefse$pValue[i]), sep='')))
axis(2, at=seq(1,index-2,2)-0.4, labels=do.call(expression, formatted_p), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 
axis(side=1, at=short_ticks+0.03, label=rep('',length(short_ticks)), tck=-0.01)
axis(side=1, at=minor_ticks+0.42, label=rep('',length(minor_ticks)), tck=-0.01)
axis(side=1, at=minor_ticks+1.42, label=rep('',length(minor_ticks)), tck=-0.01)
axis(side=1, at=minor_ticks+2.42, label=rep('',length(minor_ticks)), tck=-0.01)

mtext('d', side=2, line=2, las=2, adj=13, padj=-11, cex=1.0, font=2)







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




