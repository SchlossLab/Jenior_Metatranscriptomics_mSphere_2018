# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Output plot name
plot_file1 <- 'results/supplement/figures/figure_S5.pdf'
plot_file2 <- 'results/supplement/figures/figure_S6.pdf'

# Input Metabolomes
metabolome_file <- 'data/metabolome/metabolomics.tsv'

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
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2'), ] # Remove contaminated sample
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL

# Metabolomes
metabolome <- metabolome[,!colnames(metabolome) %in% c('CefC5M2')] # Contaminated sample
metabolome <- metabolome[,!colnames(metabolome) %in% c('GfC1M1','GfC1M2','GfC1M3', 
                                                       'GfC2M1','GfC2M2','GfC2M3', 
                                                       'GfC3M1','GfC3M2','GfC3M3', 
                                                       'GfC4M1','GfC4M2','GfC4M3', 
                                                       'GfC5M1','GfC5M2','GfC5M3', 
                                                       'GfC6M1','GfC6M2','GfC6M3')] # Germfree samples 
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL

# Prepare to merge and combine data
metabolites <- rownames(metabolome)
metabolites <- as.data.frame(cbind(metabolites, rep('metabolite', length(metabolites)), seq(1, length(metabolites), 1)))
rownames(metabolites) <- metabolites$metabolites
metabolites$label <- paste0(as.character(metabolites$V2),as.character(metabolites$V3))
metabolites$V2 <- NULL
metabolites$V3 <- NULL
metabolome <- clean_merge(metabolites, metabolome)
rownames(metabolome) <- metabolome$label
metabolome$label <- NULL
metabolome$metabolites <- NULL
metabolome <- t(metabolome)
metabolome_metadata <- clean_merge(metadata, metabolome)
metabolome_susceptibility <- subset(metabolome_metadata, infection == 'mock')
metabolome_susceptibility$infection <- NULL
metabolome_susceptibility$abx <- NULL
metabolome_susceptibility$susceptibility <- factor(metabolome_susceptibility$susceptibility)
metabolome_infection <- subset(metabolome_metadata, abx != 'none')
metabolome_infection$abx <- NULL
metabolome_infection$susceptibility <- NULL
metabolome_infection$infection <- factor(metabolome_infection$infection)
rm(metadata, metabolome, metabolome_metadata)

#--------------------------------------------------------------------#

# susceptibility - mock infection only
# Determine optimal ntree and mtry
factor1 <- as.vector(levels(metabolome_susceptibility$susceptibility))[1]
factor2 <- as.vector(levels(metabolome_susceptibility$susceptibility))[2]
num1 <- round(length(which(metabolome_susceptibility$susceptibility == factor1)) * 0.623)
num2 <- round(length(which(metabolome_susceptibility$susceptibility == factor2)) * 0.623)
ntree_multiplier <- max(c((num1/num2), (num2/num2))) * 3 
trees <- round(ncol(metabolome_susceptibility) - 1) * ntree_multiplier
tries <- round(sqrt(ncol(metabolome_susceptibility) - 1))
rm(factor1, factor2, num1, num2, ntree_multiplier)

# Run random forest and assess predictive value
susceptibility_rf <- randomForest(metabolome_susceptibility$susceptibility ~ ., 
                                  data=metabolome_susceptibility, importance=TRUE, replace=FALSE, 
                                  do.trace=FALSE, err.rate=TRUE, ntree=trees, mtry=tries)
print(susceptibility_rf)
rm(trees, tries)

# Retreive importance and overall error rate
susceptibility_importances <- as.data.frame(importance(susceptibility_rf, type=1))
susceptibility_importances$label <- rownames(susceptibility_importances)
susceptibility_accuracy <- paste('Accuracy = ',as.character((100-round((median(susceptibility_rf$err.rate[,1])*100), 2))),'%',sep='')
rm(susceptibility_rf)

# Subset to the most important OTUs and sort, subset to top 10
susceptibility_importances <- subset(susceptibility_importances, susceptibility_importances$MeanDecreaseAccuracy > abs(min(susceptibility_importances$MeanDecreaseAccuracy)))
susceptibility_importances <- as.data.frame(susceptibility_importances[order(-susceptibility_importances$MeanDecreaseAccuracy),])
susceptibility_importances <- susceptibility_importances[1:10,]

# Reassociate metabolite names - importances
susceptibility_importances <- merge(metabolites, susceptibility_importances, by='label')
susceptibility_importances$metabolites <- gsub('_', ' ', susceptibility_importances$metabolites)
rownames(susceptibility_importances) <- susceptibility_importances$metabolites
susceptibility_importances$metabolites <- NULL

# Break into experimental groups
metabolome_susceptibility <- metabolome_susceptibility[, c(1,which(colnames(metabolome_susceptibility) %in% susceptibility_importances$label))]
metabolome_susceptible <- subset(metabolome_susceptibility, susceptibility == 'susceptible')
metabolome_susceptible$susceptibility <- NULL
metabolome_resistant <- subset(metabolome_susceptibility, susceptibility == 'resistant')
metabolome_resistant$susceptibility <- NULL
rm(metabolome_susceptibility)

#------------------#

# infection - between all antibiotic treatments
# Determine optimal ntree and mtry
factor1 <- as.vector(levels(metabolome_infection$infection))[1]
factor2 <- as.vector(levels(metabolome_infection$infection))[2]
num1 <- round(length(which(metabolome_infection$infection == factor1)) * 0.623)
num2 <- round(length(which(metabolome_infection$infection == factor2)) * 0.623)
ntree_multiplier <- max(c((num1/num2), (num2/num2))) * 3 
trees <- round(ncol(metabolome_infection) - 1) * ntree_multiplier
tries <- round(sqrt(ncol(metabolome_infection) - 1))
rm(factor1, factor2, num1, num2, ntree_multiplier)

# Run random forest and assess predictive value
infection_rf <- randomForest(metabolome_infection$infection ~ ., 
                             data=metabolome_infection, importance=TRUE, replace=FALSE, 
                             do.trace=FALSE, err.rate=TRUE, ntree=trees, mtry=tries)
print(infection_rf)
rm(trees, tries)

# Retreive importance and overall error rate
infection_importances <- as.data.frame(importance(infection_rf, type=1))
infection_importances$label <- rownames(infection_importances)
infection_accuracy <- paste('Accuracy = ',as.character((100-round((median(infection_rf$err.rate[,1])*100), 2))),'%',sep='')
rm(infection_rf)

# Subset to the most important OTUs and sort
infection_importances <- subset(infection_importances, infection_importances$MeanDecreaseAccuracy > abs(min(infection_importances$MeanDecreaseAccuracy)))
infection_importances <- as.data.frame(infection_importances[order(-infection_importances$MeanDecreaseAccuracy),])
infection_importances <- infection_importances[1:10,]

# Reassociate metabolite names - importances
infection_importances <- merge(metabolites, infection_importances, by='label')
infection_importances$metabolites <- gsub('_', ' ', infection_importances$metabolites)
rownames(infection_importances) <- infection_importances$metabolites
infection_importances$metabolites <- NULL

# Break into experimental groups
metabolome_infection <- metabolome_infection[, c(1,which(colnames(metabolome_infection) %in% infection_importances$label))]
metabolome_infected <- subset(metabolome_infection, infection == '630')
metabolome_infected$infection <- NULL
metabolome_mock <- subset(metabolome_infection, infection == 'mock')
metabolome_mock$infection <- NULL
rm(metabolome_infection)

#--------------------------------------------------------------------#

# Testing for significant change abundance
# susceptibility
susceptibility_pvalues <- c()
for (index in 1:ncol(metabolome_susceptible)){
  susceptibility_pvalues[index] <- wilcox.test(metabolome_susceptible[,index], metabolome_resistant[,index], exact=FALSE)$p.value
}
susceptibility_pvalues <- round(p.adjust(susceptibility_pvalues, method='BH'), 3)
for (index in 1:length(susceptibility_pvalues)) {
  if (susceptibility_pvalues[index] == 0) {
    susceptibility_pvalues[index] <- '< 0.001 ***'
  }
  else if (susceptibility_pvalues[index] == 0.001) {
    susceptibility_pvalues[index] <- paste('= ', as.character(susceptibility_pvalues[index]), ' ***', sep='')
  }
  else if (susceptibility_pvalues[index] <= 0.01) {
    susceptibility_pvalues[index] <- paste('= ', as.character(susceptibility_pvalues[index]), ' **', sep='')
  }
  else if (susceptibility_pvalues[index] <= 0.05) {
    susceptibility_pvalues[index] <- paste('= ', as.character(susceptibility_pvalues[index]), ' *', sep='')
  }
  else {
    susceptibility_pvalues[index] <- paste('= ', as.character(susceptibility_pvalues[index]), ' n.s.', sep='')
  }
}
susceptibility_importances$pvalues <- susceptibility_pvalues
rm(susceptibility_pvalues)

# infection
infection_pvalues <- c()
for (index in 1:ncol(metabolome_infected)){
  infection_pvalues[index] <- wilcox.test(metabolome_susceptible[,index], metabolome_resistant[,index], exact=FALSE)$p.value
}
infection_pvalues <- round(p.adjust(infection_pvalues, method='BH'), 3)
for (index in 1:length(infection_pvalues)) {
  if (infection_pvalues[index] == 0) {
    infection_pvalues[index] <- '< 0.001 ***'
  }
  else if (infection_pvalues[index] == 0.001) {
    infection_pvalues[index] <- paste('= ', as.character(infection_pvalues[index]), ' ***', sep='')
  }
  else if (infection_pvalues[index] <= 0.01) {
    infection_pvalues[index] <- paste('= ', as.character(infection_pvalues[index]), ' **', sep='')
  }
  else if (infection_pvalues[index] <= 0.05) {
    infection_pvalues[index] <- paste('= ', as.character(infection_pvalues[index]), ' *', sep='')
  }
  else {
    infection_pvalues[index] <- paste('= ', as.character(infection_pvalues[index]), ' n.s.', sep='')
  }
}
infection_importances$pvalues <- infection_pvalues
rm(infection_pvalues)

#--------------------------------------------------------------------#

# Sort groups according to decreasing importance
susceptibility_importances <- as.data.frame(susceptibility_importances[order(susceptibility_importances$MeanDecreaseAccuracy),])
metabolome_susceptible <- metabolome_susceptible[susceptibility_importances$label]
metabolome_resistant <- metabolome_resistant[susceptibility_importances$label]
infection_importances <- as.data.frame(infection_importances[order(infection_importances$MeanDecreaseAccuracy),])
metabolome_infection <- metabolome_infection[infection_importances$label]
metabolome_mock <- metabolome_mock[infection_importances$label]

# Reassociate metabolite names - concentrations
metabolome_susceptible <- t(metabolome_susceptible)
metabolome_susceptible <- merge(metabolome_susceptible, metabolites, by.y='label', by.x='row.names', all.x=TRUE)
metabolome_susceptible$metabolites <- gsub('_', ' ', metabolome_susceptible$metabolites)
rownames(metabolome_susceptible) <- metabolome_susceptible$metabolites
metabolome_susceptible$metabolites <- NULL
metabolome_susceptible$label <- NULL
metabolome_susceptible$Row.names <- NULL
metabolome_susceptible <- as.data.frame(t(metabolome_susceptible))
metabolome_resistant <- t(metabolome_resistant)
metabolome_resistant <- merge(metabolome_resistant, metabolites, by.y='label', by.x='row.names', all.x=TRUE)
metabolome_resistant$metabolites <- gsub('_', ' ', metabolome_resistant$metabolites)
rownames(metabolome_resistant) <- metabolome_resistant$metabolites
metabolome_resistant$metabolites <- NULL
metabolome_resistant$label <- NULL
metabolome_resistant$Row.names <- NULL
metabolome_resistant <- as.data.frame(t(metabolome_resistant))
metabolome_infected <- t(metabolome_infected)
metabolome_infected <- merge(metabolome_infected, metabolites, by.y='label', by.x='row.names', all.x=TRUE)
metabolome_infected$metabolites <- gsub('_', ' ', metabolome_infected$metabolites)
rownames(metabolome_infected) <- metabolome_infected$metabolites
metabolome_infected$metabolites <- NULL
metabolome_infected$label <- NULL
metabolome_infected$Row.names <- NULL
metabolome_infected <- as.data.frame(t(metabolome_infected))
metabolome_mock <- t(metabolome_mock)
metabolome_mock <- merge(metabolome_mock, metabolites, by.y='label', by.x='row.names', all.x=TRUE)
metabolome_mock$metabolites <- gsub('_', ' ', metabolome_mock$metabolites)
rownames(metabolome_mock) <- metabolome_mock$metabolites
metabolome_mock$metabolites <- NULL
metabolome_mock$label <- NULL
metabolome_mock$Row.names <- NULL
metabolome_mock <- as.data.frame(t(metabolome_mock))
rm(metabolites)

#--------------------------------------------------------------------#

# Set up plotting environment
pdf(file=plot_file1, width=11, height=6)
layout(matrix(c(1,2,2), nrow=1, ncol=3, byrow=TRUE))

#---------------------#

# Susceptibility
# RF mean decrease accuracy
#par(mar=c(1.8,3,1,1), xaxs='i', xaxt='n', xpd=FALSE, mgp=c(2,0.2,0))
#dotchart(susceptibility_importances$MeanDecreaseAccuracy, labels=rownames(susceptibility_importances),
#         lcolor=NA, cex=1.2, color='black', 
#         xlab='', xlim=c(3,5), pch=19, lwd=3)
#segments(x0=rep(3, 15), y0=c(1:10), x1=rep(5, 15), y1=c(1:10), lty=2) # Dotted lines
#legend('bottomright', legend=susceptibility_accuracy, pt.cex=0, cex=1.2, bty='n')
#par(xaxt='s')
#axis(side=1, at=c(3:5), labels=c(0,4,5), cex.axis=1.2, tck=-0.025)
#axis.break(1, 3.2, style='slash')
#mtext('Mean Decrease Accuracy', side=1, padj=1.8, cex=0.9)
#mtext('A', side=2, line=2, las=2, adj=1, padj=-13.2, cex=1.7)



pdf(file='results/supplement/figures/figure_rf.pdf', width=9, height=5)
par(mar=c(1.8,3,1,1), xaxs='i', xaxt='n', xpd=FALSE, mgp=c(2.5,0.5,0))
dotchart(infection_importances$MeanDecreaseAccuracy, labels=rownames(infection_importances),
         lcolor=NA, cex=1.2, color='black', 
         xlab='', xlim=c(7,14), pch=19, lwd=3)
segments(x0=rep(7, 10), y0=c(1:10), x1=rep(14, 10), y1=c(1:10), lty=2)
legend('bottomright', legend=infection_accuracy, pt.cex=0, cex=1.2, bty='n')
par(xaxt='s')
axis(side=1, at=c(7:14), labels=c(0,8:14), cex.axis=1.2, tck=-0.05)
axis.break(1, 7.5, style='slash')
mtext('Mean Decrease Accuracy', side=1, padj=1.8, cex=0.9)
dev.off()




# OTU abundance differences
#par(mar=c(3,19,1,1), xaxs='r', mgp=c(2,1,0))
#plot(1, type='n', ylim=c(0.8, (ncol(metabolome_susceptible)*2)-0.8), xlim=c(0,3), 
#     ylab='', xlab='Abundance', xaxt='n', yaxt='n', cex.lab=1.4)
#index <- 1
#for(i in colnames(metabolome_susceptible)){
#  stripchart(at=index+0.35, metabolome_susceptible[,i], 
#             pch=21, bg='firebrick1', method='jitter', jitter=0.15, cex=1.7, lwd=0.5, add=TRUE)
#  stripchart(at=index-0.35, colonized_preabx_shared[,i], 
#             pch=21, bg='dodgerblue1', method='jitter', jitter=0.15, cex=1.7, lwd=0.5, add=TRUE)
#  if (i != colnames(metabolome_susceptible)[length(colnames(metabolome_susceptible))]){
#    abline(h=index+1, lty=2)
#  }
#  segments(median(metabolome_susceptible[,i]), index+0.6, median(metabolome_susceptible[,i]), index+0.1, lwd=2.5) #adds line for median
#  segments(median(colonized_preabx_shared[,i]), index-0.6, median(colonized_preabx_shared[,i]), index-0.1, lwd=2.5)
#  index <- index + 2
#}
#axis(side=1, at=c(0:3), label=c('0','10','100','1000'), cex.axis=1.2, tck=-0.02)
#minors <- c(0.1,0.28,0.44,0.58,0.7,0.8,0.88,0.94,0.98)
#axis(side=1, at=minors, label=rep('',length(minors)), tck=-0.01)
#axis(side=1, at=minors+1, label=rep('',length(minors)), tck=-0.01)
#axis(side=1, at=minors+2, label=rep('',length(minors)), tck=-0.01)
#legend('topright', legend=c('Cleared', 'Colonized'),
#       pch=c(21, 21), pt.bg=c('firebrick1','dodgerblue1'), bg='white', pt.cex=1.7, cex=1.2)
#axis(2, at=seq(1,index-2,2)+0.6, labels=colnames(metabolome_susceptible), las=1, line=-0.5, tick=F, cex.axis=1.4)

#dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()



