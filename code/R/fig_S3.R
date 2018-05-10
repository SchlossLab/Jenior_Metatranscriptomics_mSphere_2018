
# Set up environment
starting_dir <- getwd()
rm(list=ls())
gc()

# Load in functions
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_mSphere_2018/code/R/functions.R')

# Output names
plot_file <- 'results/supplement/figures/figure_S3.pdf'
table_file <- 'results/supplement/tables/table_S2.tsv'

# Input Metabolomes
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'

# Input Metadata
metadata <- 'data/metadata.tsv'

#----------------#

# Read in data

# Metabolomes
metabolome <- read.delim(metabolome, sep='\t', header=TRUE)

# Metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data

# Metadata
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL

# Metabolomes
metabolome$BIOCHEMICAL <- gsub('_', ' ', metabolome$BIOCHEMICAL)
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.data.frame(t(metabolome))
colnames(metabolome) <- make.names(colnames(metabolome))

#-------------------------------------------------------------------------------------------------------------------------#

# Stats
# Prep data
metabolome <- clean_merge(metadata, metabolome)
metabolome <- subset(metabolome, abx != 'germfree')
resistant_metabolome <- subset(metabolome, abx == 'none')
resistant_metabolome$susceptibility <- NULL
resistant_metabolome$abx <- NULL
resistant_metabolome$infection <- NULL
metabolome <- subset(metabolome, abx != 'none')
metabolome$susceptibility <- NULL

# Stickland metabolites
stick <- c('proline','trans.4.hydroxyproline','glycine','X5.aminovalerate')
other_aa <- c('glycine', 'alanine', 'valine', 'cysteine', 'proline','leucine','isoleucine','methionine','tryptophan','phenylalanine',
              'lysine','arginine',
              'serine','threonine','tyrosine','asparagine','glutamine',
              'glutamate','aspartate')



# Get metabolites in Stickland fermentation pathway
stickland <- metabolome[, c('abx','infection', stick)]
otherAA <- metabolome[, c('abx','infection',other_aa)]
resistant_stickland <- resistant_metabolome[, stick]
resistant_otherAA <- resistant_metabolome[, other_aa]
rm(metabolome, resistant_metabolome, stick, other_aa)

# Reformat metabolite names
colnames(stickland) <- c('abx','infection','Proline','4-Hydroxyproline','Glycine','5-Aminovalerate')
colnames(resistant_stickland) <- c('Proline','4-Hydroxyproline','Glycine','5-Aminovalerate')

#------------#

# Subset groups to mock and infected
stickland_mock <- subset(stickland, infection == 'mock')
stickland_mock$infection <- NULL
stickland_mock_strep <- subset(stickland_mock, abx == 'streptomycin')
stickland_mock_strep$abx <- NULL
stickland_mock_cef <- subset(stickland_mock, abx == 'cefoperazone')
stickland_mock_cef$abx <- NULL
stickland_mock_clinda <- subset(stickland_mock, abx == 'clindamycin')
stickland_mock_clinda$abx <- NULL
rm(stickland_mock)
stickland_infected <- subset(stickland, infection == '630')
stickland_infected$infection <- NULL
stickland_infected_strep <- subset(stickland_infected, abx == 'streptomycin')
stickland_infected_strep$abx <- NULL
stickland_infected_cef <- subset(stickland_infected, abx == 'cefoperazone')
stickland_infected_cef$abx <- NULL
stickland_infected_clinda <- subset(stickland_infected, abx == 'clindamycin')
stickland_infected_clinda$abx <- NULL
rm(stickland_infected)

otherAA_mock <- subset(otherAA, infection == 'mock')
otherAA_mock$infection <- NULL
otherAA_mock_strep <- subset(otherAA_mock, abx == 'streptomycin')
otherAA_mock_strep$abx <- NULL
otherAA_mock_cef <- subset(otherAA_mock, abx == 'cefoperazone')
otherAA_mock_cef$abx <- NULL
otherAA_mock_clinda <- subset(otherAA_mock, abx == 'clindamycin')
otherAA_mock_clinda$abx <- NULL
rm(otherAA_mock)
otherAA_infected <- subset(otherAA, infection == '630')
otherAA_infected$infection <- NULL
otherAA_infected <- subset(otherAA, infection == '630')
otherAA_infected$infection <- NULL
otherAA_infected_strep <- subset(otherAA_infected, abx == 'streptomycin')
otherAA_infected_strep$abx <- NULL
otherAA_infected_cef <- subset(otherAA_infected, abx == 'cefoperazone')
otherAA_infected_cef$abx <- NULL
otherAA_infected_clinda <- subset(otherAA_infected, abx == 'clindamycin')
otherAA_infected_clinda$abx <- NULL
rm(otherAA_infected)
rm(metadata, stickland, otherAA)

#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
alphabet <- c('A','B','C','D')
pdf(file=plot_file, width=10, height=8)
layout(matrix(c(1,2,
                3,4),
              nrow=2, ncol=2, byrow=TRUE))
for (x in 1:ncol(resistant_stickland)) {
  metabolitePlot(resistant_stickland, 
                 stickland_mock_strep, stickland_infected_strep,
                 stickland_mock_cef, stickland_infected_cef,
                 stickland_mock_clinda, stickland_infected_clinda,
                 x, alphabet[x])
}
dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Assemble supplemental table
summary_table <- c('strep_infected_q25','strep_infected_median','strep_infected_q75',
                   'strep_mock_q25','strep_mock_median','strep_mock_q75',
                   'strep_mock_vs_infected_pval',
                   'cef_infected_q25','cef_infected_median','cef_infected_q75',
                   'cef_mock_q25','cef_mock_median','cef_mock_q75',
                   'cef_mock_vs_infected_pval',
                   'clinda_infected_q25','clinda_infected_median','clinda_infected_q75',
                   'clinda_mock_q25','clinda_mock_median','clinda_mock_q75',
                   'clinda_mock_vs_infected_pval',
                   'resistant_q25','resistant_median','resistant_q75',
                   'resistant_vs_strep_infected_pval','resistant_vs_strep_mock_pval',
                   'resistant_vs_cef_infected_pval','resistant_vs_cef_mock_pval',
                   'resistant_vs_clinda_infected_pval','resistant_vs_clinda_mock_pval')

for (x in 1:ncol(otherAA_mock_strep)) {
  
  # Statistics
  infected_strep_iqr <- as.vector(quantile(otherAA_mock_strep[,x])[2:4])
  mock_strep_iqr <- as.vector(quantile(otherAA_infected_strep[,x])[2:4])
  strep_pval <- wilcox.test(otherAA_mock_strep[,x], otherAA_infected_strep[,x], exact=FALSE)$p.value
  strep_pval <- round(strep_pval, 4)
  
  infected_cef_iqr <- as.vector(quantile(otherAA_mock_cef[,x])[2:4])
  mock_cef_iqr <- as.vector(quantile(otherAA_infected_cef[,x])[2:4])
  cef_pval <- wilcox.test(otherAA_mock_cef[,x], otherAA_infected_cef[,x], exact=FALSE)$p.value
  cef_pval <- round(cef_pval, 4)
  
  infected_clinda_iqr <- as.vector(quantile(otherAA_mock_clinda[,x])[2:4])
  mock_clinda_iqr <- as.vector(quantile(otherAA_infected_clinda[,x])[2:4])
  clinda_pval <- wilcox.test(otherAA_mock_clinda[,x], otherAA_infected_clinda[,x], exact=FALSE)$p.value
  clinda_pval <- round(clinda_pval, 4)
  
  resistant_quartile <- as.vector(quantile(resistant_otherAA[,x])[2:4])
  resistant_pvals <- p.adjust(c(wilcox.test(resistant_otherAA[,x], otherAA_infected_strep[,x], exact=FALSE)$p.value,
                              wilcox.test(resistant_otherAA[,x], otherAA_mock_strep[,x], exact=FALSE)$p.value,
                              wilcox.test(resistant_otherAA[,x], otherAA_infected_cef[,x], exact=FALSE)$p.value,
                              wilcox.test(resistant_otherAA[,x], otherAA_mock_cef[,x], exact=FALSE)$p.value,
                              wilcox.test(resistant_otherAA[,x], otherAA_infected_clinda[,x], exact=FALSE)$p.value,
                              wilcox.test(resistant_otherAA[,x], otherAA_mock_clinda[,x], exact=FALSE)$p.value),
                              method='BH')
  resistant_pvals <- round(resistant_pvals, 4)
  
  # Aggregate column
  summary_vect <- c(resistant_quartile, resistant_pvals,
                    infected_strep_iqr, mock_strep_iqr, strep_pval,
                    infected_cef_iqr, mock_cef_iqr, cef_pval,
                    infected_clinda_iqr, mock_clinda_iqr, clinda_pval)
  summary_table <- as.data.frame(cbind(summary_table, summary_vect))

}
colnames(summary_table) <- c('Statistic',
                             'Glycine', 'Alanine', 'Valine', 'Cysteine', 'Proline','Leucine','Isoleucine','Methionine','Tryptophan','Phenylalanine',
                             'Lysine','Arginine',
                             'Serine','Threonine','Tyrosine','Asparagine','Glutamine',
                             'Glutamate','Aspartate')

# Write supplementary table to file
write.table(summary_table, file=table_file, sep='\t', row.names=FALSE, quote=FALSE)

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()
