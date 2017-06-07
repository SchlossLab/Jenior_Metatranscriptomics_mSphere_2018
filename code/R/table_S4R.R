
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

# Define output file
supptable_S4 <'results/supplement/tables/table_S4.tsv'

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
noabx_mock_scores <- read.delim(noabx_mock_scores, sep='\t', header=T, row.names=1)

#--------------------------------------------------------------------------------#

# Format data
# Metadata
metadata <- subset(metadata, type != 'germfree')
metadata$susceptibility <- NULL

# Untargeted metabolomics
metabolome <- subset(metabolome, KEGG != 'NA')
metabolome_annotation <- metabolome[,1:5]
metabolome$SUPER_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$PUBCHEM <- NULL
metabolome$BIOCHEMICAL <- NULL
metabolome <- metabolome[match(unique(metabolome$KEGG), metabolome$KEGG),]
rownames(metabolome) <- metabolome$KEGG
metabolome$KEGG <- NULL
#metabolome <- as.data.frame(10 ^ metabolome) # undo log10 transformation

# Metabolite scores
strep_630_scores <- subset(strep_630_scores, percentile >= 70)
strep_630_scores$consumption_score <- NULL
strep_630_scores$production_score <- NULL
strep_630_scores$percentile <- NULL
#temp_pos <- subset(strep_630_scores, cumulative_metabolite_score > 0)
#temp_pos$cumulative_metabolite_score <- 2 ^ temp_pos$cumulative_metabolite_score
#temp_neg <- subset(strep_630_scores, cumulative_metabolite_score < 0)
#temp_neg$cumulative_metabolite_score <- -(2 ^ abs(temp_neg$cumulative_metabolite_score))
#strep_630_scores <- rbind(temp_pos, temp_neg)
strep_mock_scores <- subset(strep_mock_scores, percentile >= 70)

strep_mock_scores$consumption_score <- NULL
strep_mock_scores$production_score <- NULL
strep_mock_scores$percentile <- NULL
#temp_pos <- subset(strep_mock_scores, cumulative_metabolite_score > 0)
#temp_pos$cumulative_metabolite_score <- 2 ^ temp_pos$cumulative_metabolite_score
#temp_neg <- subset(strep_mock_scores, cumulative_metabolite_score < 0)
#temp_neg$cumulative_metabolite_score <- -(2 ^ abs(temp_neg$cumulative_metabolite_score))
#strep_mock_scores <- rbind(temp_pos, temp_neg)
cef_630_scores <- subset(cef_630_scores, percentile >= 70)

cef_630_scores$consumption_score <- NULL
cef_630_scores$production_score <- NULL
cef_630_scores$percentile <- NULL
#temp_pos <- subset(cef_630_scores, cumulative_metabolite_score > 0)
#temp_pos$cumulative_metabolite_score <- 2 ^ temp_pos$cumulative_metabolite_score
#temp_neg <- subset(cef_630_scores, cumulative_metabolite_score < 0)
#temp_neg$cumulative_metabolite_score <- -(2 ^ abs(temp_neg$cumulative_metabolite_score))
#cef_630_scores <- rbind(temp_pos, temp_neg)
cef_mock_scores <- subset(cef_mock_scores, percentile >= 70)

cef_mock_scores$consumption_score <- NULL
cef_mock_scores$production_score <- NULL
cef_mock_scores$percentile <- NULL
#temp_pos <- subset(cef_mock_scores, cumulative_metabolite_score > 0)
#temp_pos$cumulative_metabolite_score <- 2 ^ temp_pos$cumulative_metabolite_score
#temp_neg <- subset(cef_mock_scores, cumulative_metabolite_score < 0)
#temp_neg$cumulative_metabolite_score <- -(2 ^ abs(temp_neg$cumulative_metabolite_score))
#cef_mock_scores <- rbind(temp_pos, temp_neg)
clinda_630_scores <- subset(clinda_630_scores, percentile >= 70)

clinda_630_scores$consumption_score <- NULL
clinda_630_scores$production_score <- NULL
clinda_630_scores$percentile <- NULL
#temp_pos <- subset(clinda_630_scores, cumulative_metabolite_score > 0)
#temp_pos$cumulative_metabolite_score <- 2 ^ temp_pos$cumulative_metabolite_score
#temp_neg <- subset(clinda_630_scores, cumulative_metabolite_score < 0)
#temp_neg$cumulative_metabolite_score <- -(2 ^ abs(temp_neg$cumulative_metabolite_score))
#clinda_630_scores <- rbind(temp_pos, temp_neg)
clinda_mock_scores <- subset(clinda_mock_scores, percentile >= 70)

clinda_mock_scores$consumption_score <- NULL
clinda_mock_scores$production_score <- NULL
clinda_mock_scores$percentile <- NULL
#temp_pos <- subset(clinda_mock_scores, cumulative_metabolite_score > 0)
#temp_pos$cumulative_metabolite_score <- 2 ^ temp_pos$cumulative_metabolite_score
#temp_neg <- subset(clinda_mock_scores, cumulative_metabolite_score < 0)
#temp_neg$cumulative_metabolite_score <- -(2 ^ abs(temp_neg$cumulative_metabolite_score))
#clinda_mock_scores <- rbind(temp_pos, temp_neg)
noabx_mock_scores <- subset(noabx_mock_scores, percentile >= 70)

noabx_mock_scores$consumption_score <- NULL
noabx_mock_scores$production_score <- NULL
noabx_mock_scores$percentile <- NULL
#temp_pos <- subset(noabx_mock_scores, cumulative_metabolite_score > 0)
#temp_pos$cumulative_metabolite_score <- 2 ^ temp_pos$cumulative_metabolite_score
#temp_neg <- subset(noabx_mock_scores, cumulative_metabolite_score < 0)
#temp_neg$cumulative_metabolite_score <- -(2 ^ abs(temp_neg$cumulative_metabolite_score))
#noabx_mock_scores <- rbind(temp_pos, temp_neg)

# Combine with metadata and subset
metabolome <- merge(metadata, t(metabolome), by='row.names')
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
metabolome$cage <- NULL
metabolome$mouse <- NULL
metabolome$gender <- NULL
metabolome$type <- NULL
metabolome_630 <- subset(metabolome, infection == '630')
metabolome_630$infection <- NULL
metabolome_630 <- aggregate(metabolome_630[, 2:ncol(metabolome_630)], list(metabolome_630$abx), median)
rownames(metabolome_630) <- metabolome_630$Group.1
metabolome_630$Group.1 <- NULL
metabolome_630 <- as.data.frame(t(metabolome_630))
metabolome_mock <- subset(metabolome, infection == 'mock')
metabolome_mock$infection <- NULL
metabolome_mock <- aggregate(metabolome_mock[, 2:ncol(metabolome_mock)], list(metabolome_mock$abx), median)
rownames(metabolome_mock) <- metabolome_mock$Group.1
metabolome_mock$Group.1 <- NULL
metabolome_mock <- as.data.frame(t(metabolome_mock))
rm(metabolome, metadata)

# Merge metabolome medians and metabolite scores
temp <- metabolome_630
temp$cefoperazone <- NULL
temp$clindamycin <- NULL
combined_strep_630 <- merge(strep_630_scores, temp, by='row.names')
rownames(combined_strep_630) <- combined_strep_630$Row.names
combined_strep_630$Row.names <- NULL
rm(strep_630_scores)
temp <- metabolome_mock
temp$cefoperazone <- NULL
temp$clindamycin <- NULL
temp$none <- NULL
combined_strep_mock <- merge(strep_mock_scores, temp, by='row.names')
rownames(combined_strep_mock) <- combined_strep_mock$Row.names
combined_strep_mock$Row.names <- NULL
rm(strep_mock_scores)
temp <- metabolome_630
temp$streptomycin <- NULL
temp$clindamycin <- NULL
combined_cef_630 <- merge(cef_630_scores, temp, by='row.names')
rownames(combined_cef_630) <- combined_cef_630$Row.names
combined_cef_630$Row.names <- NULL
rm(cef_630_scores)
temp <- metabolome_mock
temp$streptomycin <- NULL
temp$clindamycin <- NULL
temp$none <- NULL
combined_cef_mock <- merge(cef_mock_scores, temp, by='row.names')
rownames(combined_cef_mock) <- combined_cef_mock$Row.names
combined_cef_mock$Row.names <- NULL
rm(cef_mock_scores)
temp <- metabolome_630
temp$cefoperazone <- NULL
temp$streptomycin <- NULL
combined_clinda_630 <- merge(clinda_630_scores, temp, by='row.names')
rownames(combined_clinda_630) <- combined_clinda_630$Row.names
combined_clinda_630$Row.names <- NULL
rm(clinda_630_scores)
temp <- metabolome_mock
temp$cefoperazone <- NULL
temp$streptomycin <- NULL
temp$none <- NULL
combined_clinda_mock <- merge(clinda_mock_scores, temp, by='row.names')
rownames(combined_clinda_mock) <- combined_clinda_mock$Row.names
combined_clinda_mock$Row.names <- NULL
rm(clinda_mock_scores)
temp <- metabolome_mock
temp$cefoperazone <- NULL
temp$clindamycin <- NULL
temp$streptomycin <- NULL
combined_noabx_mock <- merge(noabx_mock_scores, temp, by='row.names')
rownames(combined_noabx_mock) <- combined_noabx_mock$Row.names
combined_noabx_mock$Row.names <- NULL
rm(noabx_mock_scores)
rm(temp,metabolome_630,metabolome_mock)

# Fit to general linear models and calculate correlation, assemble table
strep_630_fit <- glm(streptomycin ~ cumulative_metabolite_score, data=combined_strep_630)
strep_630_r <- cor.test(x=combined_strep_630$cumulative_metabolite_score, y=combined_strep_630$streptomycin, method='spearman', exact=FALSE)$estimate
strep_630_p <- cor.test(x=combined_strep_630$cumulative_metabolite_score, y=combined_strep_630$streptomycin, method='spearman', exact=FALSE)$p.value
strep_mock_fit <- glm(streptomycin ~ cumulative_metabolite_score, data=combined_strep_mock)
strep_mock_r <- cor.test(x=combined_strep_mock$cumulative_metabolite_score, y=combined_strep_mock$streptomycin, method='spearman', exact=FALSE)$estimate
strep_mock_p <- cor.test(x=combined_strep_mock$cumulative_metabolite_score, y=combined_strep_mock$streptomycin, method='spearman', exact=FALSE)$p.value

cef_630_fit <- glm(cefoperazone ~ cumulative_metabolite_score, data=combined_cef_630)
cef_630_r <- cor.test(x=combined_cef_630$cumulative_metabolite_score, y=combined_cef_630$cefoperazone, method='spearman', exact=FALSE)$estimate
cef_630_p <- cor.test(x=combined_cef_630$cumulative_metabolite_score, y=combined_cef_630$cefoperazone, method='spearman', exact=FALSE)$p.value
cef_mock_fit <- glm(cefoperazone ~ cumulative_metabolite_score, data=combined_cef_mock)
cef_mock_r <- cor.test(x=combined_cef_mock$cumulative_metabolite_score, y=combined_cef_mock$cefoperazone, method='spearman', exact=FALSE)$estimate
cef_mock_p <- cor.test(x=combined_cef_mock$cumulative_metabolite_score, y=combined_cef_mock$cefoperazone, method='spearman', exact=FALSE)$p.value

clinda_630_fit <- glm(clindamycin ~ cumulative_metabolite_score, data=combined_clinda_630)
clinda_630_r <- cor.test(x=combined_clinda_630$cumulative_metabolite_score, y=combined_clinda_630$clindamycin, method='spearman', exact=FALSE)$estimate
clinda_630_p <- cor.test(x=combined_clinda_630$cumulative_metabolite_score, y=combined_clinda_630$clindamycin, method='spearman', exact=FALSE)$p.value
clinda_mock_fit <- glm(clindamycin ~ cumulative_metabolite_score, data=combined_clinda_mock)
clinda_mock_r <- cor.test(x=combined_clinda_mock$cumulative_metabolite_score, y=combined_clinda_mock$clindamycin, method='spearman', exact=FALSE)$estimate
clinda_mock_p <- cor.test(x=combined_clinda_mock$cumulative_metabolite_score, y=combined_clinda_mock$clindamycin, method='spearman', exact=FALSE)$p.value

noabx_mock_fit <- glm(none ~ cumulative_metabolite_score, data=combined_noabx_mock)
noabx_mock_r <- cor.test(x=combined_noabx_mock$cumulative_metabolite_score, y=combined_noabx_mock$none, method='spearman', exact=FALSE)$estimate
noabx_mock_p <- cor.test(x=combined_noabx_mock$cumulative_metabolite_score, y=combined_noabx_mock$none, method='spearman', exact=FALSE)$p.value






# Write to a file


supptable_S4_file


#----------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()
