
# Start with a blank slate
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')

# Define files
metadata <- 'data/metadata.tsv'
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'
shared <- 'data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'

# Calculate IQR of inv-Simpson diversity for a table
tableDiv <- function(datatable) {
  div_vect <- c()
  for (x in 1:nrow(datatable)) {
    div_vect <- c(div_vect, diversity(datatable[x,], index='invsimpson'))
  }
  div_iqr <- round(as.vector(quantile(div_vect)[2:4]), 3)
  return(div_iqr)
}

# Test differences in inv-Simpson diversity for 2 tables
testDiv <- function(datatable1, datatable2) {
  div_vect1 <- c()
  div_vect2 <- c()
  for (x in 1:nrow(datatable1)) {
    div_vect1 <- c(div_vect1, diversity(datatable1[x,], index='invsimpson'))
    div_vect2 <- c(div_vect2, diversity(datatable2[x,], index='invsimpson'))
  }
  pval <- round(wilcox.test(div_vect1, div_vect2, exact=FALSE)$p.value, 4)
  return(pval)
}

# Find IQR of bray-curtis dissimilarity between groups
distBray <- function(datatable1, datatable2) {
  joined <- as.data.frame(rbind(datatable1, datatable2))
  bray_dist <- as.matrix(vegdist(joined, method='bray', diag=TRUE, upper=TRUE))
  bray_dist <- bray_dist[(nrow(datatable1)+1):(nrow(bray_dist)),]
  bray_dist <- bray_dist[,1:nrow(datatable1)]
  dist_vect <- c()
  for (x in 1:ncol(bray_dist)) {
    dist_vect <- c(bray_dist[,x], dist_vect)
  }
  brayIQR <- round(as.vector(quantile(dist_vect)[2:4]), 3)
  return(brayIQR)
}

# Significant difference in bray-curtis dissimilarity between groups
testBray <- function(datatable1, datatable2) {
  joined <- as.data.frame(rbind(datatable1, datatable2))
  group_vect <- c(rep('group1', nrow(datatable1)), rep('group2', nrow(datatable2)))
  bray_dist <- vegdist(joined, method='bray')
  pval <- adonis(bray_dist ~ group_vect, joined, perm=999)$aov.tab[[6]][1]
  pval <- round(pval, 3)
  return(pval)
}
#----------------------------------------#

# Read data
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2','StrepC4M1'), ] # Remove possible contamination
metabolome <- read.delim(metabolome, sep='\t', header=T, row.names=1)
metabolome <- metabolome[, !colnames(metabolome) %in% c('CefC5M2','StrepC4M1')] # Remove possible contamination
shared <- read.delim(shared, sep='\t', header=T, row.names=2)
shared <- shared[!rownames(shared) %in% c('CefC5M2','StrepC4M1'), ] # Remove possible contamination

# Format data
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$type <- NULL
metadata$susceptibility <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome <- as.data.frame(t(metabolome))
colnames(metabolome) <- make.names(colnames(metabolome))
shared$label <- NULL
shared$numOtus <- NULL

# Merge datasets
shared <- merge(metadata, shared, by='row.names')
rownames(shared) <- shared$Row.names
shared$Row.names <- NULL
metabolome <- merge(metadata, metabolome, by='row.names')
rownames(metabolome) <- metabolome$Row.names
metabolome$Row.names <- NULL
rm(metadata)

# Subset datasets
strep <- subset(metabolome, abx == 'streptomycin')
strep$abx <- NULL
strep_metabolome_mock <- subset(strep, infection == 'mock')
strep_metabolome_mock$infection <- NULL
strep_metabolome_630 <- subset(strep, infection == '630')
strep_metabolome_630$infection <- NULL
cef <- subset(metabolome, abx == 'cefoperazone')
cef$abx <- NULL
cef_metabolome_mock <- subset(cef, infection == 'mock')
cef_metabolome_mock$infection <- NULL
cef_metabolome_630 <- subset(cef, infection == '630')
cef_metabolome_630$infection <- NULL
clinda <- subset(metabolome, abx == 'clindamycin')
clinda$abx <- NULL
clinda_metabolome_mock <- subset(clinda, infection == 'mock')
clinda_metabolome_mock$infection <- NULL
clinda_metabolome_630 <- subset(clinda, infection == '630')
clinda_metabolome_630$infection <- NULL
conv <- subset(metabolome, abx == 'none')
conv$abx <- NULL
conv_metabolome_mock <- subset(conv, infection == 'mock')
conv_metabolome_mock$infection <- NULL
rm(metabolome)
strep <- subset(shared, abx == 'streptomycin')
strep$abx <- NULL
strep_shared_mock <- subset(strep, infection == 'mock')
strep_shared_mock$infection <- NULL
strep_shared_630 <- subset(strep, infection == '630')
strep_shared_630$infection <- NULL
cef <- subset(shared, abx == 'cefoperazone')
cef$abx <- NULL
cef_shared_mock <- subset(cef, infection == 'mock')
cef_shared_mock$infection <- NULL
cef_shared_630 <- subset(cef, infection == '630')
cef_shared_630$infection <- NULL
clinda <- subset(shared, abx == 'clindamycin')
clinda$abx <- NULL
clinda_shared_mock <- subset(clinda, infection == 'mock')
clinda_shared_mock$infection <- NULL
clinda_shared_630 <- subset(clinda, infection == '630')
clinda_shared_630$infection <- NULL
conv <- subset(shared, abx == 'none')
conv$abx <- NULL
conv_shared_mock <- subset(conv, infection == 'mock')
conv_shared_mock$infection <- NULL
rm(shared)
rm(strep, cef, clinda, conv)

#----------------------------------------#

# Calculate stats - antibiotics
# IQR of inv-Simpson diversity
tableDiv(strep_shared_mock)
tableDiv(cef_shared_mock)
tableDiv(clinda_shared_mock)
tableDiv(conv_shared_mock)
# Signifianct difference in lpha-diversity
testDiv(strep_shared_mock, conv_shared_mock)
testDiv(cef_shared_mock, conv_shared_mock)
testDiv(clinda_shared_mock, conv_shared_mock)

# Bray-Curtis - both data sets
distBray(strep_shared_mock, conv_shared_mock)
distBray(cef_shared_mock, conv_shared_mock)
distBray(clinda_shared_mock, conv_shared_mock)
distBray(strep_metabolome_mock, conv_metabolome_mock)
distBray(cef_metabolome_mock, conv_metabolome_mock)
distBray(clinda_metabolome_mock, conv_metabolome_mock)
# Signifianct difference in beta-diversity
testBray(strep_shared_mock, conv_shared_mock)
testBray(cef_shared_mock, conv_shared_mock)
testBray(clinda_shared_mock, conv_shared_mock)
testBray(strep_metabolome_mock, conv_metabolome_mock)
testBray(cef_metabolome_mock, conv_metabolome_mock)
testBray(clinda_metabolome_mock, conv_metabolome_mock)

#----------------------------------------#

# Calculate stats - infection
# IQR of inv-Simpson diversity
tableDiv(cef_shared_630)
tableDiv(strep_shared_630)
tableDiv(clinda_shared_630)
# Signifianct difference in lpha-diversity
testDiv(cef_shared_630, cef_shared_mock)
testDiv(strep_shared_630, strep_shared_mock)
testDiv(clinda_shared_630, clinda_shared_mock)

# Bray-Curtis - both data sets
distBray(cef_shared_630, cef_shared_mock)
distBray(strep_shared_630, strep_shared_mock)
distBray(clinda_shared_630, clinda_shared_mock)
distBray(cef_metabolome_630, cef_metabolome_mock)
distBray(strep_metabolome_630, strep_metabolome_mock)
distBray(clinda_metabolome_630, clinda_metabolome_mock)
# Signifianct difference in beta-diversity
testBray(cef_shared_630, cef_shared_mock)
testBray(strep_shared_630, strep_shared_mock)
testBray(clinda_shared_630, clinda_shared_mock)
testBray(cef_metabolome_630, cef_metabolome_mock)
testBray(strep_metabolome_630, strep_metabolome_mock)
testBray(clinda_metabolome_630, clinda_metabolome_mock)

#----------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()
