
# Set up environment
starting_dir <- getwd()
rm(list=ls())
gc()

# Load in functions
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')

# Output plot name
fig3_plot <- 'results/supplement/figures/figure_3.pdf'

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
stickland_aa <- c('proline','trans.4.hydroxyproline',
                  'aspartate','serine','threonine','cysteine',
                  'arginine','glutamate','asparagine','glycine','valine','leucine','isoleucine',
                  'methionine','phenylalanine','alanine')
stickland_carboxy <- c('X5.aminovalerate','alanine')

# Get metabolites in Stickland fermentation pathway
metabolome_aa <- metabolome[, c('abx','infection', stickland_aa)]
metabolome_carboxy <- metabolome[, c('abx','infection',stickland_carboxy)]
resistant_metabolome_aa <- resistant_metabolome[, stickland_aa]
resistant_metabolome_carboxy <- resistant_metabolome[, stickland_carboxy]
rm(metabolome, resistant_metabolome, stickland_aa, stickland_carboxy)

#------------#

# Reformat metabolite names
colnames(metabolome_aa) <- c('abx','infection',
                             'Proline','4-Hydroxyproline',
                             'Aspartate','Serine','Threonine','Cysteine',
                             'Arginine','Glutamate','Asparagine','Glycine','Valine','Leucine','Isoleucine',
                             'Methionine','Phenylalanine','Alanine')
colnames(resistant_metabolome_aa) <- c('Proline','4-Hydroxyproline',
                                       'Aspartate','Serine','Threonine','Cysteine',
                                       'Arginine','Glutamate','Asparagine','Glycine','Valine','Leucine','Isoleucine',
                                       'Methionine','Phenylalanine','Alanine')
colnames(metabolome_carboxy) <- c('abx','infection',
                                  '5-Aminovalerate')
colnames(resistant_metabolome_carboxy) <- c('5-Aminovalerate','alanine')
resistant_metabolome_carboxy$alanine <- NULL

#------------#

# Subset groups to mock and infected
metabolome_aa_mock <- subset(metabolome_aa, infection == 'mock')
metabolome_aa_mock$infection <- NULL
metabolome_aa_mock_strep <- subset(metabolome_aa_mock, abx == 'streptomycin')
metabolome_aa_mock_strep$abx <- NULL
metabolome_aa_mock_cef <- subset(metabolome_aa_mock, abx == 'cefoperazone')
metabolome_aa_mock_cef$abx <- NULL
metabolome_aa_mock_clinda <- subset(metabolome_aa_mock, abx == 'clindamycin')
metabolome_aa_mock_clinda$abx <- NULL
rm(metabolome_aa_mock)
metabolome_aa_infected <- subset(metabolome_aa, infection == '630')
metabolome_aa_infected$infection <- NULL
metabolome_aa_infected_strep <- subset(metabolome_aa_infected, abx == 'streptomycin')
metabolome_aa_infected_strep$abx <- NULL
metabolome_aa_infected_cef <- subset(metabolome_aa_infected, abx == 'cefoperazone')
metabolome_aa_infected_cef$abx <- NULL
metabolome_aa_infected_clinda <- subset(metabolome_aa_infected, abx == 'clindamycin')
metabolome_aa_infected_clinda$abx <- NULL
rm(metabolome_aa_infected)

metabolome_carboxy_mock <- subset(metabolome_carboxy, infection == 'mock')
metabolome_carboxy_mock$infection <- NULL
metabolome_carboxy_mock_strep <- subset(metabolome_carboxy_mock, abx == 'streptomycin')
metabolome_carboxy_mock_strep$abx <- NULL
metabolome_carboxy_mock_cef <- subset(metabolome_carboxy_mock, abx == 'cefoperazone')
metabolome_carboxy_mock_cef$abx <- NULL
metabolome_carboxy_mock_clinda <- subset(metabolome_carboxy_mock, abx == 'clindamycin')
metabolome_carboxy_mock_clinda$abx <- NULL
rm(metabolome_carboxy_mock)
metabolome_carboxy_infected <- subset(metabolome_carboxy, infection == '630')
metabolome_carboxy_infected$infection <- NULL
metabolome_carboxy_infected <- subset(metabolome_carboxy, infection == '630')
metabolome_carboxy_infected$infection <- NULL
metabolome_carboxy_infected_strep <- subset(metabolome_carboxy_infected, abx == 'streptomycin')
metabolome_carboxy_infected_strep$abx <- NULL
metabolome_carboxy_infected_cef <- subset(metabolome_carboxy_infected, abx == 'cefoperazone')
metabolome_carboxy_infected_cef$abx <- NULL
metabolome_carboxy_infected_clinda <- subset(metabolome_carboxy_infected, abx == 'clindamycin')
metabolome_carboxy_infected_clinda$abx <- NULL
rm(metabolome_carboxy_infected)

rm(metadata, metabolome_aa, metabolome_carboxy)

#-------------------------------------------------------------------------------------------------------------------------#
alphabet <- c('A','B','C','D','E','F','G','H','I','J','K','L','M',
              'N','O','P','Q','R','S','T','U','V','W','X','Y','Z')
aa_plot <- 'results/supplement/figures/figure_S3A.pdf'
carboxy_plot <- 'results/supplement/figures/figure_S3B.pdf'

# Plot the figure
pdf(file=aa_plot, width=20, height=16)
layout(matrix(c(1,2,3,4,
                5,6,7,8,
                9,10,11,12,
                13,14,15,16),
              nrow=4, ncol=4, byrow=TRUE))
for (x in 1:ncol(resistant_metabolome_aa)) {
  metabolitePlot(resistant_metabolome_aa, 
                 metabolome_aa_mock_strep, metabolome_aa_infected_strep,
                 metabolome_aa_mock_cef, metabolome_aa_infected_cef,
                 metabolome_aa_mock_clinda, metabolome_aa_infected_clinda,
                 x, '')
}
dev.off()

pdf(file=carboxy_plot, width=6, height=5)
metabolitePlot(resistant_metabolome_carboxy,
               metabolome_carboxy_mock_strep, metabolome_carboxy_infected_strep,
               metabolome_carboxy_mock_cef, metabolome_carboxy_infected_cef,
               metabolome_carboxy_mock_clinda, metabolome_carboxy_infected_clinda,
               1, '')
dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()
