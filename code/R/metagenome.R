
# Read in data and format
cef_annotation <- read.delim('/Users/schloss/Desktop/Jenior_812/data/cefoperazone.metagenome.pathways', sep='\t', header=T)
strep_annotation <- read.delim('/Users/schloss/Desktop/Jenior_812/data/streptomycin.metagenome.pathways', sep='\t', header=T)
clinda_annotation <- read.delim('/Users/schloss/Desktop/Jenior_812/data/clindamycin.metagenome.pathways', sep='\t', header=T)

cef_annotation <- cef_annotation[order(cef_annotation$count, decreasing = TRUE),] 
cef_annotation <- subset(cef_annotation, category != 'Metabolic pathways')
cef_annotation$count <- (cef_annotation$count / sum(cef_annotation$count)) * 100
cef_annotation <- head(cef_annotation, 10)

strep_annotation <- strep_annotation[order(strep_annotation$count, decreasing = TRUE),] 
strep_annotation <- subset(strep_annotation, category != 'Metabolic pathways')
strep_annotation$count <- (strep_annotation$count / sum(strep_annotation$count)) * 100
strep_annotation <- head(strep_annotation, 10)

clinda_annotation <- clinda_annotation[order(clinda_annotation$count, decreasing = TRUE),] 
clinda_annotation <- subset(clinda_annotation, category != 'Metabolic pathways')
clinda_annotation$count <- (clinda_annotation$count / sum(clinda_annotation$count)) * 100
clinda_annotation <- head(clinda_annotation, 10)


# Plot the data



# need to make into stacked barplot



par(las=1, mar=c(3,18,1,1), mgp=c(1.7,0.7,0))
barplot(rev(cef_annotation$count), names.arg=rev(cef_annotation$category), 
        col='dodgerblue', horiz=TRUE, xlim=c(0,15), xaxt='n', xlab='Percent of Annotated Genes')
axis(side=1, at=c(0,5,10,15), labels=c('0%','5%','10%','15%'), tick=TRUE)
abline(v=c(0,15), lwd=2)
abline(v=c(5,10), lty=2)
abline(h=12.46, lwd=1.5)

par(las=1, mar=c(3,18,1,1), mgp=c(1.7,0.7,0))
barplot(rev(strep_annotation$count), names.arg=rev(strep_annotation$category), 
        col='darkorange', horiz=TRUE, xlim=c(0,15), xaxt='n', xlab='Percent of Annotated Genes')
axis(side=1, at=c(0,5,10,15), labels=c('0%','5%','10%','15%'), tick=TRUE)
abline(v=c(0,15), lwd=2)
abline(v=c(5,10), lty=2)
abline(h=12.46, lwd=1.5)

par(las=1, mar=c(3,18,1,1), mgp=c(1.7,0.7,0))
barplot(rev(clinda_annotation$count), names.arg=rev(clinda_annotation$category), 
        col='forestgreen', horiz=TRUE, xlim=c(0,15), xaxt='n', xlab='Percent of Annotated Genes')
axis(side=1, at=c(0,5,10,15), labels=c('0%','5%','10%','15%'), tick=TRUE)
abline(v=c(0,15), lwd=2)
abline(v=c(5,10), lty=2)
abline(h=12.46, lwd=1.5)

