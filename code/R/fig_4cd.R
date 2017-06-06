

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
interaction <- subset(interaction, percentile == 90)
interaction$percentile <- NULL
community <- subset(community, percentile == 90)
community$percentile <- NULL

# Single interaction
# Find top edges of competition & cooperation
competition <- interaction[order(-interaction$interaction_score),][c(1:8),]
cooperation <- interaction[order(interaction$interaction_score),][c(1:8),]
interaction <- rbind(competition, cooperation)
interaction <- interaction[order(interaction$interaction_score),]
interaction$cd630_cefoperazone_score <- -interaction$cd630_cefoperazone_score
rownames(interaction) <- interaction$compound_name
interaction$compound_name <- NULL
interaction$cd630_cefoperazone_score <- NULL
interaction$lactobacillus_630_score <- NULL
interaction$ratio <- NULL
interaction$magnitude <- NULL
interaction$percentile <- NULL
rm(competition, cooperation)

# Whole community
# Find community-wide trends in production & consumption
consumption <- community[order(-community$cumulative_metabolite_score),][c(1:10),]
production <- community[order(community$cumulative_metabolite_score),][c(1:10),]
community <- rbind(consumption, production)
rownames(community) <- community$compound_name
community$compound_name <- NULL
community$consumption_score <- NULL
community$production_score <- NULL
rm(consumption, production)

#-------------------------------------------------------------------------------------------------------------------------#

# Generate figure
pdf(file=plot_file, width=7, height=7)
layout(matrix(c(1,
                2),
              nrow=2, ncol=1, byrow = TRUE))

#------------------#

# Single Interaction

par(mar=c(3,15,2,1), mgp=c(2.5, 0.75, 0), las=1, xaxs='i')
barplot(interaction[,1], xaxt='n', xlim=c(-15,15), horiz=TRUE, 
        xlab='', ylab='', col=c(rep('steelblue4',8),rep('firebrick3',8)), cex.names=0.8)
box()
abline(v=0, lwd=2)
abline(h=9.7, lty=5, lwd=1.5)


axis(1, at=seq(0,14,2), label=seq(0,14,2))
mtext('Interaction Score', side=1, padj=2.2, cex=0.75)
mtext('c', side=2, padj=-10, adj=18, font=2)
legend('topleft', legend=c(expression(italic('C. difficile')),italic('Lactobacillus')), pt.bg=c('black','white'), 
       pch=22, pt.cex=1.5, cex=0.9)




#------------------#

# Whole Community





dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)}
setwd(starting_dir)
rm(list=ls())
gc()