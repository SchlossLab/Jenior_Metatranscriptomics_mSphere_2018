
# Set up environment
rm(list=ls())
gc()

# Load in functions
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define files

# Normalized Metatranscriptomes
noabx_normalized_reads <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/conv_normalized.tsv'
cef_normalized_reads <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/cef_normalized.tsv'
clinda_normalized_reads <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/clinda_normalized.tsv'
strep_normalized_reads <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/strep_normalized.tsv'

# KEGG taxonomy IDs
kegg_tax <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/kegg_taxonomy.tsv'

# Taxonomy colors
tax_colors <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/taxonomy_color.tsv'

# Output plot
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/figures/figure_2.pdf'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data

# Normalized Metatranscriptomes
noabx_normalized_reads <- read.delim(noabx_normalized_reads, sep='\t', header=TRUE, row.names=5)
cef_normalized_reads <- read.delim(cef_normalized_reads, sep='\t', header=TRUE, row.names=6)
clinda_normalized_reads <- read.delim(clinda_normalized_reads, sep='\t', header=TRUE, row.names=6)
strep_normalized_reads <- read.delim(strep_normalized_reads, sep='\t', header=TRUE, row.names=6)

# KEGG organism file
kegg_tax <- read.delim(kegg_tax, sep='\t', header=TRUE)
kegg_tax[] <- lapply(kegg_tax, as.character)

# Taxonomy colors
tax_colors <- read.delim(tax_colors, sep='\t', header=TRUE)
tax_colors[] <- lapply(tax_colors, as.character)

#-------------------------------------------------------------------------------------------------------------------------#

# Screen for those genes that were able to be annotated
noabx_annotated <- noabx_normalized_reads[!rownames(noabx_normalized_reads) %in% rownames(noabx_normalized_reads[grep('unknown_\\d', noabx_normalized_reads$gene),]), ]
cef_annotated <- cef_normalized_reads[!rownames(cef_normalized_reads) %in% rownames(cef_normalized_reads[grep('unknown_\\d', cef_normalized_reads$gene),]), ]
clinda_annotated <- clinda_normalized_reads[!rownames(clinda_normalized_reads) %in% rownames(clinda_normalized_reads[grep('unknown_\\d', clinda_normalized_reads$gene),]), ]
strep_annotated <- strep_normalized_reads[!rownames(strep_normalized_reads) %in% rownames(strep_normalized_reads[grep('unknown_\\d', strep_normalized_reads$gene),]), ]
rm(noabx_normalized_reads, cef_normalized_reads, clinda_normalized_reads, strep_normalized_reads)

# Merge with abx-treated with untreated
cef_annotated$ko <- NULL
cef_annotated$gene <- NULL
cef_annotated$pathway <- NULL
cef_noabx_annotated <- clean_merge(cef_annotated, noabx_annotated)
clinda_annotated$ko <- NULL
clinda_annotated$gene <- NULL
clinda_annotated$pathway <- NULL
clinda_noabx_annotated <- clean_merge(clinda_annotated, noabx_annotated)
strep_annotated$ko <- NULL
strep_annotated$gene <- NULL
strep_annotated$pathway <- NULL
strep_noabx_annotated <- clean_merge(strep_annotated, noabx_annotated)
rm(cef_annotated, strep_annotated, clinda_annotated, noabx_annotated)




#-------------------------------------------------------------------------------------------------------------------------------------#

#Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()





