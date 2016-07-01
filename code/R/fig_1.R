


deps <- c('wesanderson', 'shape');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep, deps)

fox <- wes_palette("FantasticFox")
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/figure_1.pdf'

#-------------------------------------------------------------------------------------------------------------------------------------#

# Set up multi-panel figure
pdf(file=plot_file, width=15, height=9)
layout(matrix(c(1,2,
                3,3), 
              nrow=2, ncol=2, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------------------#

# Figure 1A.  Timeline of mouse experiments

# Create an empty plot
par(mar=c(0,0,0,0))
plot(0, type='n', axes=F, xlab='', ylab='', xlim=c(-4.75,4), ylim=c(-2,5))

# Abx in drinking water timeline
rect(xleft=-4, ybottom=2.8, xright=0, ytop=3.2, col=fox[3], border='black')
Arrows(x0=-4, y0=3, x1=3.5, y1=3, lwd=4, arr.type='triangle', arr.length=0.75, arr.width=0.4)
segments(x0=c(-4,0,2,2.75), y0=c(3.5,3.5,3.5,3.5), x1=c(-4,0,2,2.75), y1=c(2.5,2.5,2.5,2.5), lwd=4)
segments(x0=c(-4,-3,-2,-1,1), y0=c(3.25,3.25,3.25,3.25,3.25), x1=c(-4,-3,-2,-1,1), y1=c(2.75,2.75,2.75,2.75,2.75), lwd=2)
points(x=c(2,2.75), y=c(4,4), pch=25, bg=c('white','black'), col='black', cex=2.5)
text(x=c(-4,0,2,2.75), y=c(2.2,2.2,2.2,2.2), c('Day -7', 'Day -2', 'Day 0', '18 hrs'), cex=0.9)
text(x=-4.5, y=3.2, 'Cefoperazone', cex=0.8)
text(x=-4.5, y=2.95, 'or', font=2)
text(x=-4.5, y=2.7, 'Streptomycin', cex=0.8)

# IP injection abx timeline
Arrows(x0=-4, y0=0, x1=-1.5, y1=0, lwd=4, arr.type='triangle', arr.length=0.75, arr.width=0.4)
segments(x0=c(-4,-3,-2.25), y0=c(-0.5,-0.5,-0.5), x1=c(-4,-3,-2.25), y1=c(0.5,0.5,0.5), lwd=4)
points(x=c(-4,-3,-2.25), y=c(1,1,1), pch=c(25,25,25), bg=c(fox[5],'white','black'), col='black', cex=2.5)
text(x=c(-4,-3,-2.25), y=c(-0.8,-0.8,-0.8), c('Day -1', 'Day 0', '18 hrs'), cex=0.9)
text(x=-4.5, y=0, 'Clindamycin', cex=0.8)

# Legend
legend(x=0, y=0.7, legend=expression('Antibiotic in Drinking Water', 'IP Injection of Antibiotic',paste(italic('C. difficile'), ' Spore Gavage'), 'Sacrifice & Necropsy'), 
       pt.bg=c(fox[3],fox[5],'white','black'), pch=c(22,25,25,25), cex=1.2, pt.cex=c(3.2,2.2,2.2,2.2), bty='n')

# Plot label
legend('topleft', legend='A', cex=2, bty='n')

#-------------------------------------------------------------------------------------------------------------------------------------#



dev.off()














# Figure 1B.  NMDS of abx treatments

# Create empty plot to add data
par(las=1, mar=c(3,3,1,1), mgp=c(2,0.8,0))
plot(pick.abx.nmds.axes$axis1, pick.abx.nmds.axes$axis2, pch=21, cex=2.5,
     bg=c('black', 'darkorange3', 'forestgreen', 'blue3')[pick.abx.nmds.axes$abx],
     xlim=c(-0.8,0.8), ylim=c(-0.8,0.8), 
     xlab='NMDS Axis 1', ylab='NMDS Axis 2')

# Add legend
legend('bottomleft', legend=c('Cefoperzone', 'Clindamycin', 'Streptomycin', 'Conventional'), 
       col=c('darkorange3', 'forestgreen', 'blue3', 'black'), pch=15, cex=0.9, pt.cex=1.4, bty = "n")

# Plot label
mtext('B', side=2, cex=1.7, font=2, las=1, padj=-5.2, line=3)

#-------------------------------------------------------------------------------------------------------------------------------------#

# Figure 1C.  Bap plot of family level abundances
par(las=1, mar=c(4,4.5,1,12), xpd=TRUE)

# When needed, use this pallete   
final_colors <- c("gold1", "orangered1", "aquamarine3", "firebrick", "forestgreen", "blue3", 
                  "mediumorchid2", "violetred4", "mediumpurple4", "dodgerblue3", "goldenrod3", "chartreuse3")

filtered_shared$abx <- NULL
filtered_shared <- (filtered_shared/ rowSums(filtered_shared)) * 100

# Plot the final formatted table
barplot(t(filtered_shared), col=final_colors, yaxt='n', ylim=c(0,100), ylab='% Relative Abundance', font=2)
box()
axis(side=2, at=seq(0,100,20), tick=TRUE)
segments(x0=rep(0,4), y0=seq(20,80,20), x1=rep(5,4), y1=seq(20,80,20), lty=2)

# Create a figure legend in the margin
taxa <- gsub('_', ', ', rownames(t(rel_abund)))
legend(5.025, 75, legend=taxa, pt.bg=final_colors, pch=22, pt.cex=1.3, cex=0.7)

# Add figure label in margin
mtext('C', side=2, cex=1.7, font=2, las=1, padj=-5.7, line=2.5)





