
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define files
# Metadata
metadata <- 'data/metadata.tsv'

# Shared file
shared <- 'data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'
taxonomy <- 'data/16S_analysis/formatted.all_treatments.0.03.cons.taxonomy'

# Metabolomes
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'

# Output plot
strep_plot <- 'results/figures/figure_strep.pdf'
cef_plot <- 'results/figures/figure_cef.pdf'
clinda_plot <- 'results/figures/figure_clinda.pdf'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data
# Metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

# 16S data
shared <- read.delim(shared, sep='\t', header=TRUE, row.names=2)
taxonomy <- read.delim(taxonomy, sep='\t', header=TRUE, row.names=1)

# Metabolomes
metabolome <- read.delim(metabolome, sep='\t', header=TRUE)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data
# Metadata
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$susceptibility <- NULL

#----------------#

# 16S data
#Subset groups
shared$label <- NULL
shared$numOtus <- NULL
taxonomy$OTU.1 <- NULL
taxonomy$genus <- gsub('_', ' ', taxonomy$genus)
# Change OTU names...
shared <- clean_merge(metadata, shared)
strep_shared <- subset(shared, abx == 'streptomycin')
strep_shared$abx <- NULL
strep_630_shared <- subset(strep_shared, infection == '630')
strep_630_shared$infection <- NULL
strep_mock_shared <- subset(strep_shared, infection == 'mock')
strep_mock_shared$infection <- NULL
rm(strep_shared)
cef_shared <- subset(shared, abx == 'cefoperazone')
cef_shared$abx <- NULL
cef_630_shared <- subset(cef_shared, infection == '630')
cef_630_shared$infection <- NULL
cef_mock_shared <- subset(cef_shared, infection == 'mock')
cef_mock_shared$infection <- NULL
rm(cef_shared)
clinda_shared <- subset(shared, abx == 'clindamycin')
clinda_shared$abx <- NULL
clinda_630_shared <- subset(clinda_shared, infection == '630')
clinda_630_shared$infection <- NULL
clinda_mock_shared <- subset(clinda_shared, infection == 'mock')
clinda_mock_shared$infection <- NULL
rm(clinda_shared)

#----------------#

# Metabolomes
metabolome$PUBCHEM <- NULL
metabolome$BIOCHEMICAL <- gsub('_', ' ', metabolome$BIOCHEMICAL)
substr(metabolome$BIOCHEMICAL, 1, 1) <- toupper(substr(metabolome$BIOCHEMICAL, 1, 1))
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$KEGG <- NULL
metabolite_annotation <- metabolome[,c(1:2)]
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.data.frame(t(metabolome))
metabolome <- 10 ^ metabolome # Untransform
metabolome <- clean_merge(metadata, metabolome)
rm(metadata)

# Subset antibiotics
strep_metabolome <- subset(metabolome, abx == 'streptomycin')
strep_metabolome$abx <- NULL
strep_630_metabolome <- subset(strep_metabolome, infection == '630')
strep_630_metabolome$infection <- NULL
strep_mock_metabolome <- subset(strep_metabolome, infection == 'mock')
strep_mock_metabolome$infection <- NULL
rm(strep_metabolome)
cef_metabolome <- subset(metabolome, abx == 'cefoperazone')
cef_metabolome$abx <- NULL
cef_630_metabolome <- subset(cef_metabolome, infection == '630')
cef_630_metabolome$infection <- NULL
cef_mock_metabolome <- subset(cef_metabolome, infection == 'mock')
cef_mock_metabolome$infection <- NULL
rm(cef_metabolome)
clinda_metabolome <- subset(metabolome, abx == 'clindamycin')
clinda_metabolome$abx <- NULL
clinda_630_metabolome <- subset(clinda_metabolome, infection == '630')
clinda_630_metabolome$infection <- NULL
clinda_mock_metabolome <- subset(clinda_metabolome, infection == 'mock')
clinda_mock_metabolome$infection <- NULL
rm(clinda_metabolome)

#-------------------------------------------------------------------------------------------------------------------------#

# Correlation analysis
# Calculate differences
strep_shared <- (strep_mock_shared - strep_630_shared) * -1
cef_shared <- (cef_mock_shared - cef_630_shared) * -1
clinda_shared <- (clinda_mock_shared - clinda_630_shared) * -1
strep_metabolome <- (strep_mock_metabolome - strep_630_metabolome) * -1
cef_metabolome <- (cef_mock_metabolome - cef_630_metabolome) * -1
clinda_metabolome <- (clinda_mock_metabolome - clinda_630_metabolome) * -1
rm(strep_mock_shared, strep_630_shared, strep_mock_metabolome, strep_630_metabolome,
   cef_mock_shared, cef_630_shared, cef_mock_metabolome, cef_630_metabolome,
   clinda_mock_shared,clinda_630_shared, clinda_mock_metabolome, clinda_630_metabolome)

# Remove groups with low variance
strep_shared <- rm_lowVar(strep_shared)
cef_shared <- rm_lowVar(cef_shared)
clinda_shared <- rm_lowVar(clinda_shared)
strep_metabolome <- rm_lowVar(strep_metabolome)
cef_metabolome <- rm_lowVar(cef_metabolome)
clinda_metabolome <- rm_lowVar(clinda_metabolome)

# Spearman correlations
strep_correlation <- cor(x=strep_shared, y=strep_metabolome, method='spearman')
cef_correlation <- cor(x=cef_shared, y=cef_metabolome, method='spearman')
clinda_correlation <- cor(x=clinda_shared, y=clinda_metabolome, method='spearman')



rm(strep_shared, strep_metabolome, cef_shared, cef_metabolome, clinda_shared, clinda_metabolome)

#-------------------------------------------------------------------------------------------------------------------------#

pdf(file=strep_plot, width=30, height=20)
heatmap.2( strep_correlation,
           col=heat_palette,
           trace='none',
           scale='none',
           symm=FALSE,
           symbreaks=FALSE,
           cexRow=2, 
           cexCol=2,
           keysize=1,
           density.info='none')
dev.off()

pdf(file=cef_plot, width=30, height=20)
heatmap.2( cef_correlation,
           col=heat_palette,
           trace='none',
           scale='none',
           symm=FALSE,
           symbreaks=FALSE,
           cexRow=2, 
           cexCol=2,
           keysize=1,
           density.info='none')
dev.off()

pdf(file=clinda_plot, width=30, height=20)
heatmap.2( clinda_correlation,
           col=heat_palette,
           trace='none',
           scale='none',
           symm=FALSE,
           symbreaks=FALSE,
           cexRow=2, 
           cexCol=2,
           keysize=1,
           density.info='none')
dev.off()





