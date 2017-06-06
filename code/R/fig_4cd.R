
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define files
interaction <- 'data/metabolic_models/cefoperazone/infected/community.files/cd630_cefoperazone.bipartite.files.and.lactobacillus_630.bipartite.files.interaction.tsv'
community <- 'data/metabolic_models/cefoperazone/infected/community.files/community_importance.tsv'

# Output plot
plot_file <- 'results/figures/figure_4cd.pdf'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data
interaction <- read.delim(interaction, sep='\t', header=TRUE, row.names=1)
interaction$compound_name <- gsub('_', ' ', interaction$compound_name)

community <- read.delim(community, sep='\t', header=TRUE, row.names=1)
community$compound_name <- gsub('_', ' ', community$compound_name)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data
# Subset to only 90th percentile
interaction <- subset(interaction, percentile >= 80)
interaction$percentile <- NULL
interaction$ratio <- NULL
interaction$magnitude <- NULL
community <- subset(community, percentile == 90)
community$percentile <- NULL

# Single interaction
# Find top edges of competition & cooperation
competition <- interaction[order(-interaction$interaction_score),][c(1:10),]
competition <- competition[order(-competition$interaction_score),][c(1:10),]
cooperation <- interaction[order(interaction$interaction_score),][c(1:10),]
cooperation <- cooperation[order(-cooperation$interaction_score),]
interaction_scores <- rbind(competition, cooperation)
interaction_scores$cd630_cefoperazone_score <- NULL
interaction_scores$lactobacillus_630_score <- NULL
competition$interaction_score <- NULL
rownames(competition) <- competition$compound_name
competition$compound_name <- NULL
filler <- as.data.frame(matrix(0, ncol = 2, nrow = 10), row.names=rep('', 10))
colnames(filler) <- colnames(competition)
competition <- rbind(filler, competition)
competition <- as.matrix(t(competition))
cooperation$interaction_score <- NULL
rownames(cooperation) <- cooperation$compound_name
cooperation$compound_name <- NULL
filler <- as.data.frame(matrix(0, ncol = 2, nrow = 10), row.names=rep('', 10))
colnames(filler) <- colnames(cooperation)
cooperation <- rbind(cooperation, filler)
interaction_scores$compound_name <- gsub('beta-','',interaction_scores$compound_name)
interaction_scores$compound_name <- gsub('alpha,alpha\'-','',interaction_scores$compound_name)

# Whole community
# Find community-wide trends in production & consumption
community$compound_name <- gsub('beta-','',community$compound_name)
community$compound_name <- gsub('O-','',community$compound_name)

consumption <- community[order(-community$cumulative_metabolite_score),][c(1:10),]
rownames(consumption) <- consumption$compound_name
consumption$compound_name <- NULL
production <- community[order(community$cumulative_metabolite_score),][c(1:10),]
production <- production[order(-production$cumulative_metabolite_score),]
rownames(production) <- production$compound_name
production$compound_name <- NULL
filler <- c(0,0,0)
names(filler) <- ''
community <- rbind(consumption, filler, production)
community$consumption_score <- NULL
community$production_score <- NULL

#-------------------------------------------------------------------------------------------------------------------------#

# Generate figure
pdf(file=plot_file, width=7, height=10)
layout(matrix(c(1,
                2),
              nrow=2, ncol=1, byrow = TRUE))

#------------------#

# Single Interaction
par(mar=c(3,10,1,4), mgp=c(2.5, 0.6, 0), las=1)
barplot(competition, xaxt='n', horiz=TRUE, xlim=c(-15,15), yaxt='n', beside=TRUE, width=0.55,
        xlab='', ylab='', col=c('firebrick3','steelblue4'), cex.names=0.8, axis.lty=1.2)
barplot(cooperation[,1], xaxt='n', horiz=TRUE, space=0.6, xlim=c(-15,15), yaxt='n',
        xlab='', ylab='', col='firebrick3', cex.names=0.8, axis.lty=1.2, add=TRUE)
barplot(cooperation[,2], xaxt='n', horiz=TRUE, space=0.6, xlim=c(-15,15), yaxt='n',
        xlab='', ylab='', col='steelblue4', cex.names=0.8, axis.lty=1.2, add=TRUE)
text(x=-12, y=c(0.2,17.4), c('Cooperation', 'Competition'), cex=0.8)
box(lwd=1.3)
abline(v=0, lwd=1.5)
abline(h=16.5, lty=5, lwd=1.2)
axis(side=1, at=seq(-15,15,5), label=seq(-15,15,5))
mtext(expression(paste('Metabolite Score (',log[2],')')), side=1, padj=2.2, cex=1.1)
mtext(c('Interaction', 'Score:'), side=4, padj=c(-24,-22.5), adj=c(0,-0.3), cex=0.75, font=2)
mtext('c', side=2, padj=-13, adj=16, font=2, cex=1.3)
legend('topleft', legend=c(expression(italic('C. difficile'),italic('Lactobacillus'))), pt.bg=c('firebrick3','steelblue4'), 
       pch=22, pt.cex=1.8, cex=1)
axis(side=2, at=seq(17.6,32.9,1.7), labels=rev(interaction_scores$compound_name[1:10]), 
     tick=FALSE, cex.axis=0.8)
axis(side=2, at=seq(1,15.4,1.6), labels=interaction_scores$compound_name[11:20], 
     tick=FALSE, cex.axis=0.8)
axis(side=4, at=seq(0.5,32.8,1.7), labels=rev(interaction_scores$interaction_score), 
     tick=FALSE, cex.axis=0.8)

#------------------#

# Whole Community
par(mar=c(9.3,3,1,1), mgp=c(1.5, 0.6, 0), las=1)
barplot(community$cumulative_metabolite_score, xaxt='n', horiz=FALSE, space=0.4, ylim=c(-16,16), yaxt='n',
        xlab='', ylab=expression(paste('Cumulative Metabolite Score (',log[2],')')), col=c(rep('darkorchid4',10),'gray',rep('chocolate3',10)), axis.lty=1.2)
box(lwd=1.3)
abline(h=0, lwd=1.5)
axis(side=2, at=seq(-16,16,4), label=seq(-16,16,4))
text(x=seq(1.5,14.1,1.4), y=-17, rownames(community)[1:10], xpd=TRUE, srt=60, pos=2, cex=0.8)
text(x=seq(17.2,29.8,1.4), y=-17, rownames(community)[12:21], xpd=TRUE, srt=60, pos=2, cex=0.8)
text(x=15, y=c(15,-15.1), c('More Consumption', 'More Production'), cex=0.8)
mtext('d', side=2, padj=-10, adj=4.2, font=2, cex=1.3)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)}
setwd(starting_dir)
rm(list=ls())
gc()

