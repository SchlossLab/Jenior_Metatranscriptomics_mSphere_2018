
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





