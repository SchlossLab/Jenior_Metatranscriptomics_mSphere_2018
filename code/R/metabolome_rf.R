# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Output plot name
plot_file <- 'results/supplement/figures/figure_S4.pdf'

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
rm(metadata, metabolome_metadata)

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

# Subset to the most important OTUs and sort
susceptibility_importances <- subset(susceptibility_importances, susceptibility_importances$MeanDecreaseAccuracy > abs(min(susceptibility_importances$MeanDecreaseAccuracy)))
susceptibility_importances <- as.data.frame(susceptibility_importances[order(-susceptibility_importances$MeanDecreaseAccuracy),])
susceptibility_importances <- susceptibility_importances[1:50,]
susceptibility_importances <- as.data.frame(susceptibility_importances[order(susceptibility_importances$MeanDecreaseAccuracy),])

# Reassociate metabolite names
susceptibility_importances <- merge(metabolites, susceptibility_importances, by='label')
rownames(susceptibility_importances) <- susceptibility_importances$metabolites
susceptibility_importances$metabolites <- NULL

# Subset important metabolite relative concentrations, replace metabolite names
metabolome_susceptibility <- metabolome_susceptibility[, c(1,which(colnames(metabolome_susceptibility) %in% susceptibility_importances$label))]
temp_vector <- metabolome_susceptibility$susceptibility
metabolome_susceptibility$susceptibility <- NULL
metabolome_susceptibility <- t(metabolome_susceptibility)
metabolome_susceptibility <- merge(metabolites, metabolome_susceptibility, by.x='label', by.y='row.names', all.y=TRUE)
rownames(metabolome_susceptibility) <- metabolome_susceptibility$metabolites
metabolome_susceptibility$metabolites <- NULL
metabolome_susceptibility$label <- NULL
metabolome_susceptibility <- as.data.frame(t(metabolome_susceptibility))
metabolome_susceptibility$susceptibility <- temp_vector
rm(temp_vector)

# Break into experimental groups
metabolome_susceptible <- subset(metabolome_susceptibility, susceptibility == 'susceptible')
metabolome_susceptible$susceptibility <- NULL
metabolome_susceptible <- metabolome_susceptible[,match(rownames(susceptibility_importances), colnames(metabolome_susceptible))] 
metabolome_resistant <- subset(metabolome_susceptibility, susceptibility == 'resistant')
metabolome_resistant$susceptibility <- NULL
metabolome_resistant <- metabolome_resistant[,match(rownames(susceptibility_importances), colnames(metabolome_resistant))]
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
infection_importances <- infection_importances[1:50,]
infection_importances <- as.data.frame(infection_importances[order(infection_importances$MeanDecreaseAccuracy),])

# Reassociate metabolite names
infection_importances <- merge(metabolites, infection_importances, by='label')
rownames(infection_importances) <- infection_importances$metabolites
infection_importances$metabolites <- NULL

# Subset important metabolite relative concentrations
metabolome_infection <- metabolome_infection[, c(1,which(colnames(metabolome_infection) %in% infection_importances$label))]

# Subset important metabolite relative concentrations, replace metabolite names
metabolome_infection <- metabolome_infection[, c(1,which(colnames(metabolome_infection) %in% infection_importances$label))]
temp_vector <- metabolome_infection$infection
metabolome_infection$infection <- NULL
metabolome_infection <- t(metabolome_infection)
metabolome_infection <- merge(metabolites, metabolome_infection, by.x='label', by.y='row.names', all.y=TRUE)
rownames(metabolome_infection) <- metabolome_infection$metabolites
metabolome_infection$metabolites <- NULL
metabolome_infection$label <- NULL
metabolome_infection <- as.data.frame(t(metabolome_infection))
metabolome_infection$infection <- temp_vector
rm(temp_vector)

# Break into experimental groups
metabolome_infected <- subset(metabolome_infection, infection == '630')
metabolome_infected$infection <- NULL
metabolome_infected <- metabolome_infected[,match(rownames(infection_importances), colnames(metabolome_infected))] 
metabolome_mock <- subset(metabolome_infection, infection == 'mock')
metabolome_mock$infection <- NULL
metabolome_mock <- metabolome_mock[,match(rownames(infection_importances), colnames(metabolome_mock))]
rm(metabolome, metabolome_infection)

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

# Subset everything to top 15 to make plotting easier
# susceptibility
metabolome_susceptible <- metabolome_susceptible[,1:15]
metabolome_resistant <- metabolome_resistant[,1:15]
susceptibility_importances <- susceptibility_importances[1:15,]
# infection
metabolome_infected <- metabolome_infected[,1:15]
metabolome_mock <- metabolome_mock[,1:15]
infection_importances <- infection_importances[1:15,]

#--------------------------------------------------------------------#

# Set up plotting environment
pdf(file=plot_file, width=11, height=6)
layout(matrix(c(1,2,2), nrow=1, ncol=3, byrow=TRUE))

#---------------------#

# Susceptibility
# RF mean decrease accuracy
par(mar=c(1.8,3,1,1), xaxs='i', xaxt='n', xpd=FALSE, mgp=c(2,0.2,0))
dotchart(susceptibility_importances$MeanDecreaseAccuracy, labels=rownames(susceptibility_importances),
         lcolor=NA, cex=1.2, color='black', 
         xlab='', xlim=c(8,15), pch=19, lwd=3)
segments(x0=rep(8, 10), y0=c(1:15), x1=rep(15, 10), y1=c(1:15), lty=2) # Dotted lines
legend('bottomright', legend=susceptibility_accuracy, pt.cex=0, cex=1.2, bty='n')
par(xaxt='s')
axis(side=1, at=c(8:15), labels=c(0,9:15), cex.axis=1.2, tck=-0.025)
axis.break(1, 8.5, style='slash')
mtext('Mean Decrease Accuracy', side=1, padj=1.8, cex=0.9)
mtext('A', side=2, line=2, las=2, adj=1, padj=-13.2, cex=1.7)

# OTU abundance differences
par(mar=c(3,19,1,1), xaxs='r', mgp=c(2,1,0))
plot(1, type='n', ylim=c(0.8, (ncol(cleared_preabx_shared)*2)-0.8), xlim=c(0,3), 
     ylab='', xlab='Abundance', xaxt='n', yaxt='n', cex.lab=1.4)
index <- 1
for(i in colnames(cleared_preabx_shared)){
  stripchart(at=index+0.35, cleared_preabx_shared[,i], 
             pch=21, bg='firebrick1', method='jitter', jitter=0.15, cex=1.7, lwd=0.5, add=TRUE)
  stripchart(at=index-0.35, colonized_preabx_shared[,i], 
             pch=21, bg='dodgerblue1', method='jitter', jitter=0.15, cex=1.7, lwd=0.5, add=TRUE)
  if (i != colnames(cleared_preabx_shared)[length(colnames(cleared_preabx_shared))]){
    abline(h=index+1, lty=2)
  }
  segments(median(cleared_preabx_shared[,i]), index+0.6, median(cleared_preabx_shared[,i]), index+0.1, lwd=2.5) #adds line for median
  segments(median(colonized_preabx_shared[,i]), index-0.6, median(colonized_preabx_shared[,i]), index-0.1, lwd=2.5)
  index <- index + 2
}
axis(side=1, at=c(0:3), label=c('0','10','100','1000'), cex.axis=1.2, tck=-0.02)
minors <- c(0.1,0.28,0.44,0.58,0.7,0.8,0.88,0.94,0.98)
axis(side=1, at=minors, label=rep('',length(minors)), tck=-0.01)
axis(side=1, at=minors+1, label=rep('',length(minors)), tck=-0.01)
axis(side=1, at=minors+2, label=rep('',length(minors)), tck=-0.01)
legend('topright', legend=c('Cleared', 'Colonized'),
       pch=c(21, 21), pt.bg=c('firebrick1','dodgerblue1'), bg='white', pt.cex=1.7, cex=1.2)
axis(2, at=seq(1,index-2,2)+0.6, labels=rownames(preabx_importances), las=1, line=-0.5, tick=F, cex.axis=1.4)
formatted_taxa <- lapply(1:nrow(preabx_importances), function(x) bquote(paste(.(preabx_importances$phylum[x]),'; ',italic(.(preabx_importances$genus[x])), sep='')))
axis(2, at=seq(1,index-2,2), labels=do.call(expression, formatted_taxa), las=1, line=-0.5, tick=F, cex.axis=1.1, font=3) 
italic_p <- lapply(1:length(preabx_importances$pvalues), function(x) bquote(paste(italic('p'), .(preabx_importances$pvalues[x]), sep=' ')))
axis(2, at=seq(1,index-2,2)-0.6, labels=do.call(expression, italic_p), las=1, line=-0.5, tick=F, cex.axis=1.2, font=3) 
mtext('B', side=2, line=2, las=2, adj=13, padj=-13, cex=1.7)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()



