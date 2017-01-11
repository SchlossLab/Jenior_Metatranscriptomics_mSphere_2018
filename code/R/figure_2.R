
# Load dependencies
deps <- c('wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Define variables
nmds_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.0.03.unique_list.thetayc.0.03.lt.ave.nmds.axes'
summary_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.0.03.unique_list.groups.summary'
shared_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.family.subsample.shared'
taxonomy_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/16S_analysis/all_treatments.family.cons.taxonomy'
cef_lefse_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2015/data/lefse/cefoperazone.final.0.03.Cdifficile630.0.03.lefse_summary'
clinda_lefse_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2015/data/lefse/clindamycin.final.0.03.Cdifficile630.0.03.lefse_summary'
strep_lefse_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2015/data/lefse/streptomycin.final.0.03.Cdifficile630.0.03.subsample.0.03.lefse_summary'
metadata_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/metadata.tsv'

# Load in data
nmds <- read.delim(nmds_file, sep='\t', header=T, row.names=1)
nmds <- nmds[!rownames(nmds) %in% c('CefC5M2'), ] # Remove contaminated sample
metadata <- read.delim(metadata_file, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2'), ] # Remove contaminated sample
summary <- read.delim(summary_file, sep='\t', header=T, row.names=2)
summary <- summary[!rownames(summary) %in% c('CefC5M2'), ] # Remove contaminated sample
taxonomy <- read.delim(taxonomy_file, sep='\t', header=T, row.names=1)
shared <- read.delim(shared_file, sep='\t', header=T, row.names=2)
shared <- shared[!rownames(shared) %in% c('CefC5M2'), ]  # Remove contaminated sample
shared$numOtus <- NULL
shared$label <- NULL
cef_lefse <- read.delim(cef_lefse_file, sep='\t', header=T, row.names=1)
clinda_lefse <- read.delim(clinda_lefse_file, sep='\t', header=T, row.names=1)
strep_lefse <- read.delim(strep_lefse_file, sep='\t', header=T, row.names=1)

rm(nmds_file, summary_file, metadata_file, cef_lefse_file, clinda_lefse_file, strep_lefse_file)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Format data

# Combine the metadata with axes and format
metadata_axes <- merge(metadata, nmds, by='row.names')
mock_axes <- subset(metadata_axes, infection == 'mock')
mock_axes <- subset(mock_axes, type != 'germfree')
infected_axes <- subset(metadata_axes, infection == '630')
infected_axes <- subset(infected_axes, type != 'germfree')




# Format membership
filtered_shared$abx <- NULL
filtered_shared <- (filtered_shared/ rowSums(filtered_shared)) * 100
taxa <- gsub('_', ', ', rownames(t(rel_abund)))



# Format LEfSe results
cef_lefse <- cef_lefse[order(cef_lefse$Class, cef_lefse$LDA),]
clinda_lefse <- clinda_lefse[order(clinda_lefse$Class, clinda_lefse$LDA),]
strep_lefse <- strep_lefse[order(strep_lefse$Class, strep_lefse$LDA),]

rm(nmds, summary, metadata)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/supplement/figures/figure_2.pdf'
pdf(file=plot_file, width=12, height=7)
layout(matrix(c(1,2,2,
                4,5,6),
              nrow=2, ncol=3, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# NMDS of treatment groups
par(las=1, mar=c(5,5,1,1))
plot(metadata_axes$axis1, metadata_axes$axis2, pch=21, cex=0,
     xlim=c(-0.8,0.8), ylim=c(-0.8,0.8), cex.lab=2, cex.axis=1.7,
     xlab='NMDS Axis 1', ylab='NMDS Axis 2')

# add mock points
points(x=mock_axes$axis1, y=mock_axes$axis2, 
       col=c(wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen', 'black', wes_palette("FantasticFox")[1])[mock_axes$abx], 
       pch=1, lwd=4, cex=2.5)
# add infected points
points(x=infected_axes$axis1, y=infected_axes$axis2, 
       col=c(wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'forestgreen', 'black', wes_palette("FantasticFox")[1])[infected_axes$abx], 
       pch=2, lwd=4, cex=2.5)

# Add legends
legend('topleft', legend=c('Streptomycin-treated', 'Cefoperzone-treated', 'Clindamycin-treated', 'No Antibiotics'), 
       col=c(wes_palette("FantasticFox")[1], wes_palette("FantasticFox")[3], wes_palette("FantasticFox")[5], 'black'), 
       pch=15, cex=1.9, pt.cex=3, bty='n')
legend('bottomleft', legend=c('Mock Infected', '630 Infected'), 
       col='black', pch=c(16,17), cex=2, pt.cex=2.5, bty='n')

mtext('A', side=2, line=2, las=2, adj=1.7, padj=-18.1, cex=1.3)

#-----------------------#

# Family-level phylotype bar chart

par(las=1, mar=c(4,4.5,1,12), xpd=TRUE)

# When needed, use this pallete   
final_colors <- c("gold1", "orangered1", "aquamarine3", "firebrick", "forestgreen", "blue3", 
                  "mediumorchid2", "violetred4", "mediumpurple4", "dodgerblue3", "goldenrod3", "chartreuse3")


# Plot the final formatted table
barplot(t(filtered_shared), col=final_colors, yaxt='n', ylim=c(0,100), ylab='% Relative Abundance', font=2)
box()
axis(side=2, at=seq(0,100,20), tick=TRUE)
segments(x0=rep(0,4), y0=seq(20,80,20), x1=rep(5,4), y1=seq(20,80,20), lty=2)

# Create a figure legend in the margin
legend(5.025, 75, legend=taxa, pt.bg=final_colors, pch=22, pt.cex=1.3, cex=0.7)

mtext('B', side=2, line=2, las=2, adj=1.7, padj=-10.5, cex=1.3)

#-----------------------#

# LEfSe results

# Cefoperazone plot
par(mar=c(5,12,3,1))
barplot(cef_lefse$LDA, horiz=TRUE, xlim=c(0,3.5), xlab='LDA', main='Cefoperazone', 
        col='firebrick', names.arg=cef_lefse$Tax, las=1, cex.main=2, cex.lab=1.5)
box()
mtext('C', side=2, line=2, las=2, adj=1.7, padj=-8, cex=1.3)

# Clindamycin plot
par(mar=c(5,12,3,1))
barplot(clinda_lefse$LDA, horiz=TRUE, xlim=c(0,5), xlab='LDA', main='Clindamycin', 
        col=c('firebrick','firebrick','firebrick','firebrick','firebrick','firebrick','firebrick','firebrick','darkblue','darkblue','darkblue'), 
        names.arg=clinda_lefse$Tax, las=1, cex.main=2, cex.lab=1.5)
box()

# Streptomycin plot
par(mar=c(5,12,3,1), xpd=TRUE)
barplot(strep_lefse$LDA, horiz=TRUE, xlim=c(0,5.5), xlab='LDA', main='Streptomycin', 
        col=c('firebrick','firebrick','firebrick','darkblue','darkblue','darkblue'), 
        names.arg=strep_lefse$Tax, las=1, cex.main=2, cex.lab=1.5)
box()

#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
dev.off()


rm(plot_file, metadata_axes, mock_axes, infected_axes, metadata_summary, 
   strep_div, cef_div, clinda_div, conv_div)
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(dep, deps, pkg)
gc()