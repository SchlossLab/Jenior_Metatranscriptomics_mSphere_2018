
# Start with clean environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')
if ('igraph' %in% installed.packages()[,'Package'] == FALSE) {install.packages('igraph', quiet=TRUE)}
library('igraph', verbose=FALSE, character.only=TRUE)

# Define input files
metabolome <- 'data/metabolome/scaled_intensities.tsv'
strep_630_scores <- 'data/metabolic_models/streptomycin/infected/community.files/community_importance.tsv'
strep_mock_scores <- 'data/metabolic_models/streptomycin/mock/community.files/community_importance.tsv'
cef_630_scores <- 'data/metabolic_models/cefoperazone/infected/community.files/community_importance.tsv'
cef_mock_scores <- 'data/metabolic_models/cefoperazone/mock/community.files/community_importance.tsv'
clinda_630_scores <- 'data/metabolic_models/clindamycin/infected/community.files/community_importance.tsv'
clinda_mock_scores <- 'data/metabolic_models/clindamycin/mock/community.files/community_importance.tsv'
noabx_mock_scores <- 'data/metabolic_models/no_antibiotics/mock/community.files/community_importance.tsv'
metadata <- 'data/metadata.tsv'

# Define output files
plot_file <- 'results/figures/figure_5.pdf'

#--------------------------------------------------------------------------------#

# Read in data
# Study metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

# Metabolomic data
metabolome <- read.delim(metabolome, sep='\t', header=T)

# Cumulative Metaboloties Scores
strep_630_scores <- read.delim(strep_630_scores, sep='\t', header=T, row.names=1)
strep_mock_scores <- read.delim(strep_mock_scores, sep='\t', header=T, row.names=1)
cef_630_scores <- read.delim(cef_630_scores, sep='\t', header=T, row.names=1)
cef_mock_scores <- read.delim(cef_mock_scores, sep='\t', header=T, row.names=1)
clinda_630_scores <- read.delim(clinda_630_scores, sep='\t', header=T, row.names=1)
clinda_mock_scores <- read.delim(clinda_mock_scores, sep='\t', header=T, row.names=1)

#--------------------------------------------------------------------------------#

# Format data
# Metadata
metadata <- subset(metadata, type != 'germfree')
metadata$susceptibility <- NULL

# Untargeted metabolomics
metabolome <- subset(metabolome, KEGG != 'NA')
metabolome_annotation <- metabolome[,1:5]
metabolome_annotation <- subset(metabolome_annotation, !duplicated(metabolome_annotation$KEGG))
rownames(metabolome_annotation) <- metabolome_annotation$KEGG
metabolome_annotation$KEGG <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$PUBCHEM <- NULL
metabolome$BIOCHEMICAL <- NULL
metabolome <- metabolome[match(unique(metabolome$KEGG), metabolome$KEGG),]
rownames(metabolome) <- metabolome$KEGG
metabolome$KEGG <- NULL
metabolome <- as.data.frame(10 ^ metabolome) # undo log10 transformation

# Metabolite scores
# Take top percentiles 
strep_630_scores <- subset(strep_630_scores, percentile >= 70)
strep_630_scores$percentile <- NULL
strep_mock_scores <- subset(strep_mock_scores, percentile >= 70)
strep_mock_scores$percentile <- NULL
cef_630_scores <- subset(cef_630_scores, percentile >= 70)
cef_630_scores$percentile <- NULL
cef_mock_scores <- subset(cef_mock_scores, percentile >= 70)
cef_mock_scores$percentile <- NULL
clinda_630_scores <- subset(clinda_630_scores, percentile >= 70)
clinda_630_scores$percentile <- NULL
clinda_mock_scores <- subset(clinda_mock_scores, percentile >= 70)
clinda_mock_scores$percentile <- NULL

# Untransform data
strep_630_scores$consumption_score <- NULL
strep_630_scores$production_score <- NULL
temp_pos <- subset(strep_630_scores, cumulative_metabolite_score > 0)
temp_pos$cumulative_metabolite_score <- 2 ^ temp_pos$cumulative_metabolite_score
temp_neg <- subset(strep_630_scores, cumulative_metabolite_score < 0)
temp_neg$cumulative_metabolite_score <- -(2 ^ abs(temp_neg$cumulative_metabolite_score))
strep_630_scores <- rbind(temp_pos, temp_neg)
strep_mock_scores$production_score <- NULL
strep_mock_scores$consumption_score <- NULL
temp_pos <- subset(strep_mock_scores, cumulative_metabolite_score > 0)
temp_pos$cumulative_metabolite_score <- 2 ^ temp_pos$cumulative_metabolite_score
temp_neg <- subset(strep_mock_scores, cumulative_metabolite_score < 0)
temp_neg$cumulative_metabolite_score <- -(2 ^ abs(temp_neg$cumulative_metabolite_score))
strep_mock_scores <- rbind(temp_pos, temp_neg)
cef_630_scores$production_score <- NULL
cef_630_scores$consumption_score <- NULL
temp_pos <- subset(cef_630_scores, cumulative_metabolite_score > 0)
temp_pos$cumulative_metabolite_score <- 2 ^ temp_pos$cumulative_metabolite_score
temp_neg <- subset(cef_630_scores, cumulative_metabolite_score < 0)
temp_neg$cumulative_metabolite_score <- -(2 ^ abs(temp_neg$cumulative_metabolite_score))
cef_630_scores <- rbind(temp_pos, temp_neg)
cef_mock_scores$production_score <- NULL
cef_mock_scores$consumption_score <- NULL
temp_pos <- subset(cef_mock_scores, cumulative_metabolite_score > 0)
temp_pos$cumulative_metabolite_score <- 2 ^ temp_pos$cumulative_metabolite_score
temp_neg <- subset(cef_mock_scores, cumulative_metabolite_score < 0)
temp_neg$cumulative_metabolite_score <- -(2 ^ abs(temp_neg$cumulative_metabolite_score))
cef_mock_scores <- rbind(temp_pos, temp_neg)
clinda_630_scores$production_score <- NULL
clinda_630_scores$consumption_score <- NULL
temp_pos <- subset(clinda_630_scores, cumulative_metabolite_score > 0)
temp_pos$cumulative_metabolite_score <- 2 ^ temp_pos$cumulative_metabolite_score
temp_neg <- subset(clinda_630_scores, cumulative_metabolite_score < 0)
temp_neg$cumulative_metabolite_score <- -(2 ^ abs(temp_neg$cumulative_metabolite_score))
clinda_630_scores <- rbind(temp_pos, temp_neg)
clinda_mock_scores$production_score <- NULL
clinda_mock_scores$consumption_score <- NULL
temp_pos <- subset(clinda_mock_scores, cumulative_metabolite_score > 0)
temp_pos$cumulative_metabolite_score <- 2 ^ temp_pos$cumulative_metabolite_score
temp_neg <- subset(clinda_mock_scores, cumulative_metabolite_score < 0)
temp_neg$cumulative_metabolite_score <- -(2 ^ abs(temp_neg$cumulative_metabolite_score))
clinda_mock_scores <- rbind(temp_pos, temp_neg)
rm(temp_pos, temp_neg)

# Merge infected and uninfected within treatment groups, find amount of difference between conditions
strep_630_scores$compound_name <- NULL
strep_scores <- merge(strep_mock_scores, strep_630_scores, by='row.names')
rownames(strep_scores) <- strep_scores$Row.names
strep_scores$Row.names <- NULL
strep_scores$score_difference <- strep_scores$cumulative_metabolite_score.y - strep_scores$cumulative_metabolite_score.x
strep_scores$cumulative_metabolite_score.x <- NULL
strep_scores$cumulative_metabolite_score.y <- NULL
temp_pos <- subset(strep_scores, score_difference > 0)
temp_pos$score_difference <- log2(temp_pos$score_difference)
temp_neg <- subset(strep_scores, score_difference < 0)
temp_neg$score_difference <- -(log2(abs(temp_neg$score_difference)))
strep_scores <- rbind(temp_pos, temp_neg)
rm(strep_630_scores, strep_mock_scores)
cef_630_scores$compound_name <- NULL
cef_scores <- merge(cef_mock_scores, cef_630_scores, by='row.names')
rownames(cef_scores) <- cef_scores$Row.names
cef_scores$Row.names <- NULL
cef_scores$score_difference <- cef_scores$cumulative_metabolite_score.y - cef_scores$cumulative_metabolite_score.x
cef_scores$cumulative_metabolite_score.x <- NULL
cef_scores$cumulative_metabolite_score.y <- NULL
temp_pos <- subset(cef_scores, score_difference > 0)
temp_pos$score_difference <- log2(temp_pos$score_difference)
temp_neg <- subset(cef_scores, score_difference < 0)
temp_neg$score_difference <- -(log2(abs(temp_neg$score_difference)))
cef_scores <- rbind(temp_pos, temp_neg)
rm(cef_630_scores, cef_mock_scores)
clinda_630_scores$compound_name <- NULL
clinda_scores <- merge(clinda_mock_scores, clinda_630_scores, by='row.names')
rownames(clinda_scores) <- clinda_scores$Row.names
clinda_scores$Row.names <- NULL
clinda_scores$score_difference <- clinda_scores$cumulative_metabolite_score.y - clinda_scores$cumulative_metabolite_score.x
clinda_scores$cumulative_metabolite_score.x <- NULL
clinda_scores$cumulative_metabolite_score.y <- NULL
temp_pos <- subset(clinda_scores, score_difference > 0)
temp_pos$score_difference <- log2(temp_pos$score_difference)
temp_neg <- subset(clinda_scores, score_difference < 0)
temp_neg$score_difference <- -(log2(abs(temp_neg$score_difference)))
clinda_scores <- rbind(temp_pos, temp_neg)
rm(clinda_630_scores, clinda_mock_scores)
rm(temp_pos, temp_neg)

# Combine metabolomics with metadata, subset, aggregate, find differences
metabolome <- merge(metadata, t(metabolome), by='row.names')
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
metabolome$cage <- NULL
metabolome$mouse <- NULL
metabolome$gender <- NULL
metabolome$type <- NULL
metabolome_strep <- subset(metabolome, abx == 'streptomycin')
metabolome_strep$abx <- NULL
metabolome_strep <- aggregate(metabolome_strep[, 2:ncol(metabolome_strep)], list(metabolome_strep$infection), median)
rownames(metabolome_strep) <- metabolome_strep$Group.1
metabolome_strep$Group.1 <- NULL
metabolome_strep <- as.data.frame(t(metabolome_strep))
metabolome_strep$conc_difference <- metabolome_strep$`630` - metabolome_strep$mock
metabolome_strep$`630` <- NULL
metabolome_strep$mock <- NULL
temp_pos <- subset(metabolome_strep, conc_difference > 0)
temp_pos$conc_difference <- log10(temp_pos$conc_difference)
temp_neg <- subset(metabolome_strep, conc_difference < 0)
temp_neg$conc_difference <- -(log10(abs(temp_neg$conc_difference)))
metabolome_strep <- rbind(temp_pos, temp_neg)
metabolome_cef <- subset(metabolome, abx == 'streptomycin')
metabolome_cef$abx <- NULL
metabolome_cef <- aggregate(metabolome_cef[, 2:ncol(metabolome_cef)], list(metabolome_cef$infection), median)
rownames(metabolome_cef) <- metabolome_cef$Group.1
metabolome_cef$Group.1 <- NULL
metabolome_cef <- as.data.frame(t(metabolome_cef))
metabolome_cef$conc_difference <- metabolome_cef$`630` - metabolome_cef$mock
metabolome_cef$`630` <- NULL
metabolome_cef$mock <- NULL
temp_pos <- subset(metabolome_cef, conc_difference > 0)
temp_pos$conc_difference <- log10(temp_pos$conc_difference)
temp_neg <- subset(metabolome_cef, conc_difference < 0)
temp_neg$conc_difference <- -(log10(abs(temp_neg$conc_difference)))
metabolome_cef <- rbind(temp_pos, temp_neg)
metabolome_clinda <- subset(metabolome, abx == 'clindamycin')
metabolome_clinda$abx <- NULL
metabolome_clinda <- aggregate(metabolome_clinda[, 2:ncol(metabolome_clinda)], list(metabolome_clinda$infection), median)
rownames(metabolome_clinda) <- metabolome_clinda$Group.1
metabolome_clinda$Group.1 <- NULL
metabolome_clinda <- as.data.frame(t(metabolome_clinda))
metabolome_clinda$conc_difference <- metabolome_clinda$`630` - metabolome_clinda$mock
metabolome_clinda$`630` <- NULL
metabolome_clinda$mock <- NULL
temp_pos <- subset(metabolome_clinda, conc_difference > 0)
temp_pos$conc_difference <- log10(temp_pos$conc_difference)
temp_neg <- subset(metabolome_clinda, conc_difference < 0)
temp_neg$conc_difference <- -(log10(abs(temp_neg$conc_difference)))
metabolome_clinda <- rbind(temp_pos, temp_neg)
rm(metabolome, metadata)
rm(temp_pos, temp_neg)

# Merge metabolome medians and metabolite scores
combined_strep <- merge(strep_scores, metabolome_strep, by='row.names')
rownames(combined_strep) <- combined_strep$Row.names
combined_strep$Row.names <- NULL
rm(strep_scores, metabolome_strep)
combined_cef <- merge(cef_scores, metabolome_cef, by='row.names')
rownames(combined_cef) <- combined_cef$Row.names
combined_cef$Row.names <- NULL
rm(cef_scores, metabolome_cef)
combined_clinda <- merge(clinda_scores, metabolome_clinda, by='row.names')
rownames(combined_clinda) <- combined_clinda$Row.names
combined_clinda$Row.names <- NULL
rm(clinda_scores, metabolome_clinda)

#---------------------------------------------------------------------------------------#

# Fit to general linear models and calculate correlation
all_groups <- rbind(combined_strep, combined_cef, combined_clinda)
all_fit <- glm(conc_difference ~ score_difference, data=all_groups)
all_r <- round(cor.test(x=all_groups$score_difference, y=all_groups$conc_difference, method='spearman', exact=FALSE)$estimate, 3)
all_p <- round(cor.test(x=all_groups$score_difference, y=all_groups$conc_difference, method='spearman', exact=FALSE)$p.value, 3)
all_groups$residuals <- residuals(all_fit)
all_groups$residuals <- (all_groups$residuals / sd(all_groups$residuals))^2
all_outliers <- all_groups[all_groups$residuals > 1.5, ]

strep_fit <- glm(conc_difference ~ score_difference, data=combined_strep)
strep_r <- round(cor.test(x=combined_strep$score_difference, y=combined_strep$conc_difference, method='spearman', exact=FALSE)$estimate, 3)
strep_p <- round(cor.test(x=combined_strep$score_difference, y=combined_strep$conc_difference, method='spearman', exact=FALSE)$p.value, 3)
combined_strep$residuals <- residuals(strep_fit)
combined_strep$residuals <- (combined_strep$residuals / sd(combined_strep$residuals))^2
strep_outliers <- combined_strep[combined_strep$residuals > 1.5, ]

cef_fit <- glm(conc_difference ~ score_difference, data=combined_cef)
cef_r <- round(cor.test(x=combined_cef$score_difference, y=combined_cef$conc_difference, method='spearman', exact=FALSE)$estimate, 3)
cef_p <- round(cor.test(x=combined_cef$score_difference, y=combined_cef$conc_difference, method='spearman', exact=FALSE)$p.value, 3)
combined_cef$residuals <- residuals(cef_fit)
combined_cef$residuals <- (combined_cef$residuals / sd(combined_cef$residuals))^2
cef_outliers <- combined_cef[combined_cef$residuals > 1.5, ]

clinda_fit <- glm(conc_difference ~ score_difference, data=combined_clinda)
clinda_r <- round(cor.test(x=combined_clinda$score_difference, y=combined_clinda$conc_difference, method='spearman', exact=FALSE)$estimate, 3)
clinda_p <- round(cor.test(x=combined_clinda$score_difference, y=combined_clinda$conc_difference, method='spearman', exact=FALSE)$p.value, 3)
combined_clinda$residuals <- residuals(clinda_fit)
combined_clinda$residuals <- (combined_clinda$residuals / sd(combined_clinda$residuals))^2
clinda_outliers <- combined_clinda[combined_clinda$residuals > 1.5, ]

# Reassociate with annotation information with outliers
all_outliers <- merge(all_outliers, metabolome_annotation, by='row.names', all.x=TRUE)
rownames(all_outliers) <- all_outliers$Row.names
all_outliers$Row.names <- NULL
strep_outliers <- merge(strep_outliers, metabolome_annotation, by='row.names', strep.x=TRUE)
rownames(strep_outliers) <- strep_outliers$Row.names
strep_outliers$Row.names <- NULL
cef_outliers <- merge(cef_outliers, metabolome_annotation, by='row.names', cef.x=TRUE)
rownames(cef_outliers) <- cef_outliers$Row.names
cef_outliers$Row.names <- NULL
clinda_outliers <- merge(clinda_outliers, metabolome_annotation, by='row.names', clinda.x=TRUE)
rownames(clinda_outliers) <- clinda_outliers$Row.names
clinda_outliers$Row.names <- NULL
rm(metabolome_annotation)
all_outliers$SUPER_PATHWAY <- gsub('_', ' ', all_outliers$SUPER_PATHWAY)
strep_outliers$compound_name <- gsub('_', ' ', strep_outliers$compound_name)
cef_outliers$compound_name <- gsub('_', ' ', cef_outliers$compound_name)
clinda_outliers$compound_name <- gsub('_', ' ', clinda_outliers$compound_name)
all_pathways <- as.vector(unique(all_outliers$SUPER_PATHWAY))
all_pathways <- all_pathways[!is.na(all_pathways)]
all_outliers_col <- c('firebrick2','firebrick2',
                      'darkorange2','darkorange2','darkorange2','darkorange2',
                      'gold2','gold2',
                      'chartreuse3','chartreuse3',
                      'gold2',
                      'darkorange2','darkorange2',
                      'gold2',
                      'dodgerblue3',
                      'darkorchid3',
                      'chartreuse3','chartreuse3',
                      'gold2',
                      'darkorange2','darkorange2',
                      'chartreuse3')
all_outliers$all_outliers_col <- all_outliers_col

#----------------------------------------#

# Set up multi-panel figure
pdf(file=plot_file, width=3.5, height=12)
layout(matrix(c(1,
                2,
                3,
                4), 
              nrow=4, ncol=1, byrow=TRUE))
par(las=1, mar=c(3,3,1,1), mgp=c(1.8,0.7,0))

# Pooled analysis
plot(x=all_groups$score_difference, y=all_groups$conc_difference, 
     xlab=expression(paste(Delta,' Cumulative Metabolite Score')), 
     ylab=expression(paste(Delta,' Median Scaled Intensity')), 
     pch=19, cex=1.1, xlim=c(-15,15), ylim=c(-10,10), col='gray40')
filledrectangle(wx=34, wy=4.7, col='gray90', mid=c(0,all_fit$coefficients[1]), angle=2.5)
abline(v=0, lty=2, col='gray60')
abline(h=0, lty=2, col='gray60')
abline(all_fit, lwd=2)
box()
points(all_groups$score_difference, all_groups$conc_difference, pch=19, col='gray40', cex=1.2)
mtext('a', side=2, line=2, las=2, adj=1, padj=-8.5, cex=1.2, font=2)
text(x=-8.5, y=6, 'All Treatment Groups', cex=1.1)
legend('bottomright', legend=c(as.expression(bquote(paste(italic('rho'),' = ',.(all_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(all_p),'**')))), 
       pt.cex=0, bty='n', col='black')
points(all_outliers$score_difference, all_outliers$conc_difference, pch=21, col='black', bg=all_outliers$all_outliers_col, cex=2, lwd=2)
legend('topleft', legend=all_pathways, ncol=2,
       pt.bg=c('firebrick2', 
               'darkorange2',
               'gold2',
               'chartreuse3',
               'dodgerblue3',
               'darkorchid3'), col='black',
       pch=21, pt.lwd=1.5, pt.cex=1.5, cex=0.8)

# streptomycin alone
plot(x=combined_strep$score_difference, y=combined_strep$conc_difference, 
     xlab=expression(paste(Delta,' Cumulative Metabolite Score')), 
     ylab=expression(paste(Delta,' Median Scaled Intensity')), 
     pch=19, cex=1.1, xlim=c(-15,15), ylim=c(-10,10), col=strep_col)
filledrectangle(wx=34, wy=5, col='gray90', mid=c(0,0.608), angle=0.1)
abline(v=0, lty=2, col='gray60')
abline(h=0, lty=2, col='gray60')
abline(strep_fit, lwd=2)
box()
points(combined_strep$score_difference, combined_strep$conc_difference, pch=19, col=strep_col, cex=1.2)
points(strep_outliers$score_difference, strep_outliers$conc_difference, pch=21, col='black', bg=strep_col, cex=2, lwd=2)
mtext('b', side=2, line=2, las=2, adj=1, padj=-8.5, cex=1.2, font=2)
legend('topleft', legend=c('Streptomycin-pretreated'), pt.cex=0, bty='n', col='black')
legend('bottomright', legend=c(as.expression(bquote(paste(italic('rho'),' = ',.(strep_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(strep_p),'')))), 
       pt.cex=0, bty='n', col='black')
text(x=c(7.163, 4.5, 11.452, 8.170, 6), 
     y=c(-8, 3.8, 6, -3.7, 8), 
     strep_outliers$compound_name, cex=0.9, col='gray30')
segments(x0=6, y0=7.5, x1=7.8, y1=4.2, col='gray30')


# cefoperazone alone
plot(x=combined_cef$score_difference, y=combined_cef$conc_difference, 
     xlab=expression(paste(Delta,' Cumulative Metabolite Score')), 
     ylab=expression(paste(Delta,' Median Scaled Intensity')), 
     pch=19, cex=1.1, xlim=c(-15,15), ylim=c(-10,10), col=cef_col)
filledrectangle(wx=34, wy=5.3, col='gray90', mid=c(0,0.309), angle=4.5)
abline(v=0, lty=2, col='gray60')
abline(h=0, lty=2, col='gray60')
abline(cef_fit, lwd=2)
box()
points(combined_cef$score_difference, combined_cef$conc_difference, pch=19, col=cef_col, cex=1.2)
points(cef_outliers$score_difference, cef_outliers$conc_difference, pch=21, col='black', bg=cef_col, cex=2, lwd=2)
mtext('c', side=2, line=2, las=2, adj=1, padj=-8.5, cex=1.2, font=2)
legend('topleft', legend=c('Cefoperazone-pretreated'), pt.cex=0, bty='n', col='black')
legend('bottomright', legend=c(as.expression(bquote(paste(italic('rho'),' = ',.(cef_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(cef_p),'*')))), 
       pt.cex=0, bty='n', col='black')
text(x=c(-12.4, 12.2, -9.578, -9.025, 13, 6, -2.6, 5.5, 5.2), 
     y=c(-9.206, -3, 3.5, -5.2, 6, 6.4, 5, -3.4, 4), 
     cef_outliers$compound_name, cex=0.9, col='gray30')

# clindamycin alone
plot(x=combined_clinda$score_difference, y=combined_clinda$conc_difference, 
     xlab=expression(paste(Delta,' Cumulative Metabolite Score')), 
     ylab=expression(paste(Delta,' Median Scaled Intensity')), 
     pch=19, cex=1.1, xlim=c(-15,15), ylim=c(-10,10), col=clinda_col)
filledrectangle(wx=34, wy=3.9, col='gray90', mid=c(0,0.271), angle=1.5)
abline(v=0, lty=2, col='gray60')
abline(h=0, lty=2, col='gray60')
abline(clinda_fit, lwd=2)
box()
points(combined_clinda$score_difference, combined_clinda$conc_difference, pch=19, col=clinda_col, cex=1.2)
points(clinda_outliers$score_difference, clinda_outliers$conc_difference, pch=21, col='black', bg=clinda_col, cex=2, lwd=2)
mtext('d', side=2, line=2, las=2, adj=1, padj=-8.5, cex=1.2, font=2)
legend('topleft', legend=c('Clindamycin-pretreated'), pt.cex=0, bty='n', col='black')
legend('bottomright', legend=c(as.expression(bquote(paste(italic('rho'),' = ',.(clinda_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(clinda_p),'*')))), 
       pt.cex=0, bty='n', col='black')
text(x=c(-6, -11, 12.2, 7, 10.5, 13.7, 4.5), 
     y=c(3, 4.3, -3.780, 4, 6, 4, -6.4), 
     clinda_outliers$compound_name, cex=0.9, col='gray30')
segments(x0=10.5, y0=5.5, x1=10.9, y1=3.8, col='gray30')

dev.off()

#----------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()
