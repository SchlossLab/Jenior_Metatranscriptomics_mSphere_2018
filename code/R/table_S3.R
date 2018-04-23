
# Start with a blank slate
rm(list=ls())
gc()

# Load dependencies
deps <- c('wesanderson', 'plotrix');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Function for population variance of columns in a matrix
pop_var <- function(data) {
  vars <- c()
  for (x in 1:ncol(data)){
    vars[x] <- sum((data[,x] - mean(data[,x]))^2) / length(data[,x])
  }
  return(vars)
}

# Function for sample variance of columns in a matrix
samp_var <- function(data) {
  vars <- c()
  for (x in 1:ncol(data)){
    vars[x] <- var(data[,x])
  }
  return(vars)
}

# Define files
metadata <- '~/Desktop/Repositories/Jenior_Modeling_mSystems_2017/data/metadata.tsv'
metabolome <- '~/Desktop/Repositories/Jenior_Modeling_mSystems_2017/data/wetlab_assays/metabolomics.scaled_intensities.tsv'
shared <- '~/Desktop/Repositories/Jenior_Modeling_mSystems_2017/data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'
cfu <- '~/Desktop/Repositories/Jenior_Modeling_mSystems_2017/data/wetlab_assays/cfu.dat'

#----------------------------------------#

# Read data
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2','StrepC4M1'), ] # Remove possible contamination
metabolome <- read.delim(metabolome, sep='\t', header=T, row.names=1)
metabolome <- metabolome[, !colnames(metabolome) %in% c('CefC5M2','StrepC4M1')] # Remove possible contamination
shared <- read.delim(shared, sep='\t', header=T, row.names=2)
shared <- shared[!rownames(shared) %in% c('CefC5M2','StrepC4M1'), ] # Remove possible contamination
cfu <- read.delim(cfu, sep='\t', header=T)

# Format data
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$type <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome <- as.data.frame(t(metabolome))
shared$label <- NULL
shared$numOtus <- NULL
shared <- log10(shared + 10)
cfu$mouse <- NULL
cfu$cfu_spore <- NULL
cfu <- subset(cfu, cage < 4 ) # Remove uninfected controls
cfu$cage <- NULL
cfu <- subset(cfu, cfu$treatment != 'conventional')
cfu <- subset(cfu, cfu$treatment != 'germfree')
cfu$treatment <- factor(cfu$treatment, levels=c('streptomycin', 'cefoperazone', 'clindamycin'))
cfu$cfu_vegetative <- log10(cfu$cfu_vegetative)

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
cef_cfu <- as.numeric(cfu[cfu$treatment == 'cefoperazone', 2])
strep_cfu <- as.numeric(cfu[cfu$treatment == 'streptomycin', 2])
clinda_cfu <- as.numeric(cfu[cfu$treatment == 'clindamycin', 2])
rm(cfu)

#----------------------------------------#

# Calculate sample variance
strep_metabolome_mock <- samp_var(strep_metabolome_mock)
strep_metabolome_630 <- samp_var(strep_metabolome_630)
cef_metabolome_mock <- samp_var(cef_metabolome_mock)
cef_metabolome_630 <- samp_var(cef_metabolome_630)
clinda_metabolome_mock <- samp_var(clinda_metabolome_mock)
clinda_metabolome_630 <- samp_var(clinda_metabolome_630)
conv_metabolome_mock <- samp_var(conv_metabolome_mock)
strep_shared_mock <- samp_var(strep_shared_mock)
strep_shared_630 <- samp_var(strep_shared_630)
cef_shared_mock <- samp_var(cef_shared_mock)
cef_shared_630 <- samp_var(cef_shared_630)
clinda_shared_mock <- samp_var(clinda_shared_mock)
clinda_shared_630 <- samp_var(clinda_shared_630)
conv_shared_mock <- samp_var(conv_shared_mock)
cfu_var <- c(var(strep_cfu),var(cef_cfu),var(clinda_cfu)) 
rm(strep_cfu,cef_cfu,clinda_cfu) 

# Find qunatiles of variance
# 16S
strep_shared_mock <- as.vector(quantile(strep_shared_mock))[2:4]
strep_shared_630 <- as.vector(quantile(strep_shared_630))[2:4]
cef_shared_mock <- as.vector(quantile(cef_shared_mock))[2:4]
cef_shared_630 <- as.vector(quantile(cef_shared_630))[2:4]
clinda_shared_mock <- as.vector(quantile(clinda_shared_mock))[2:4]
clinda_shared_630 <- as.vector(quantile(clinda_shared_630))[2:4]
conv_shared_mock <- as.vector(quantile(conv_shared_mock))[2:4]
# Metabolome
strep_metabolome_mock <- as.vector(quantile(strep_metabolome_mock))[2:4]
strep_metabolome_630 <- as.vector(quantile(strep_metabolome_630))[2:4]
cef_metabolome_mock <- as.vector(quantile(cef_metabolome_mock))[2:4]
cef_metabolome_630 <- as.vector(quantile(cef_metabolome_630))[2:4]
clinda_metabolome_mock <- as.vector(quantile(clinda_metabolome_mock))[2:4]
clinda_metabolome_630 <- as.vector(quantile(clinda_metabolome_630))[2:4]
conv_metabolome_mock <- as.vector(quantile(conv_metabolome_mock))[2:4]

# Round to 3 decimal places
cfu_var <- round(cfu_var, 3)
strep_shared_mock <- round(strep_shared_mock, 3)
strep_shared_630 <- round(strep_shared_630, 3)
cef_shared_mock <- round(cef_shared_mock, 3)
cef_shared_630 <- round(cef_shared_630, 3)
clinda_shared_mock <- round(clinda_shared_mock, 3)
clinda_shared_630 <- round(clinda_shared_630, 3)
conv_shared_mock <- round(conv_shared_mock, 3)
strep_metabolome_mock <- round(strep_metabolome_mock, 3)
strep_metabolome_630 <- round(strep_metabolome_630, 3)
cef_metabolome_mock <- round(cef_metabolome_mock, 3)
cef_metabolome_630 <- round(cef_metabolome_630, 3)
clinda_metabolome_mock <- round(clinda_metabolome_mock, 3)
clinda_metabolome_630 <- round(clinda_metabolome_630, 3)
conv_metabolome_mock <- round(conv_metabolome_mock, 3)

# Assemble Table
sample_treatment <- c('All',
                      'Streptomycin','Streptomycin',
                      'Cefoperazone','Cefoperazone',
                      'Clindamycin','Clindamycin',
                      'None',
                      'Streptomycin','Streptomycin',
                      'Cefoperazone','Cefoperazone',
                      'Clindamycin','Clindamycin',
                      'None')
sample_infection <- c('Infected',
                      'Mock','Infected',
                      'Mock','Infected',
                      'Mock','Infected',
                      'Mock',
                      'Mock','Infected',
                      'Mock','Infected',
                      'Mock','Infected',
                      'Mock')
sample_type <- c('CFU',
                 '16S','16S','16S','16S','16S','16S','16S',
                 'Metabolomics','Metabolomics','Metabolomics','Metabolomics','Metabolomics','Metabolomics','Metabolomics')
metadata_table <- cbind(sample_treatment, sample_infection, sample_type)
rm(sample_treatment, sample_infection, sample_type)
var_table <- rbind(cfu_var, 
                   strep_shared_mock,strep_shared_630,cef_shared_mock,cef_shared_630,clinda_shared_mock,clinda_shared_630,conv_shared_mock,
                   strep_metabolome_mock,strep_metabolome_630,cef_metabolome_mock,cef_metabolome_630,clinda_metabolome_mock,clinda_metabolome_630,conv_metabolome_mock)
var_table <- cbind(metadata_table, var_table)
rm(cfu_var, 
   strep_shared_mock,strep_shared_630,cef_shared_mock,cef_shared_630,clinda_shared_mock,clinda_shared_630,conv_shared_mock,
   strep_metabolome_mock,strep_metabolome_630,cef_metabolome_mock,cef_metabolome_630,clinda_metabolome_mock,clinda_metabolome_630,conv_metabolome_mock)
colnames(var_table) <- c('Pretreatment','Intection','Data_type','Q1','Median_(Q2)','Q3')
rm(metadata_table)

# Write to a file
write.table(var_table, file='results/supplement/tables/Table_s5.tsv', sep='\t', row.names=FALSE, quote=FALSE)

#----------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()
