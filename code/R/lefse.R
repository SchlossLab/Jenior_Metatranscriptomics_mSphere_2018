

cef_lefse_file <- '~/Desktop/16S/cefoperazone.final.0.03.Cdifficile630.0.03.lefse_summary'
clinda_lefse_file <- '~/Desktop/16S/clindamycin.final.0.03.Cdifficile630.0.03.lefse_summary'
strep_lefse_file <- '~/Desktop/16S/streptomycin.final.0.03.Cdifficile630.0.03.subsample.0.03.lefse_summary'
plot_file <- '~/Desktop/lefse.pdf'

# Load in data
cef_lefse <- read.delim(cef_lefse_file, sep='\t', header=T, row.names=1)
clinda_lefse <- read.delim(clinda_lefse_file, sep='\t', header=T, row.names=1)
strep_lefse <- read.delim(strep_lefse_file, sep='\t', header=T, row.names=1)

# Format data
cef_lefse <- cef_lefse[order(cef_lefse$Class, cef_lefse$LDA),]
clinda_lefse <- clinda_lefse[order(clinda_lefse$Class, clinda_lefse$LDA),]
strep_lefse <- strep_lefse[order(strep_lefse$Class, strep_lefse$LDA),]

# Plot it
pdf(file=plot_file, width=18, height=10)
layout(matrix(1:3, nrow=1, byrow=T))
par(mar=c(5,12,3,1))
barplot(cef_lefse$LDA, horiz=TRUE, xlim=c(0,3.5), xlab='LDA', main='Cefoperazone', col='firebrick', names.arg=cef_lefse$Tax, las=1, cex.main=2, cex.lab=1.2)
box()
par(mar=c(5,12,3,1))
barplot(clinda_lefse$LDA, horiz=TRUE, xlim=c(0,5), xlab='LDA', main='Clindamycin', col=c('firebrick','firebrick','firebrick','firebrick','firebrick','firebrick','firebrick','firebrick','darkblue','darkblue','darkblue'), names.arg=clinda_lefse$Tax, las=1, cex.main=2, cex.lab=1.2)
box()
par(mar=c(5,12,3,1), xpd=TRUE)
barplot(strep_lefse$LDA, horiz=TRUE, xlim=c(0,5.5), xlab='LDA', main='Streptomycin', col=c('firebrick','firebrick','firebrick','darkblue','darkblue','darkblue'), names.arg=strep_lefse$Tax, las=1, cex.main=2, cex.lab=1.2)
box()
dev.off()
