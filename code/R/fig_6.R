
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Define files




# Output plot
plot_file <- 'results/figures/figure_6.pdf'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data


#-------------------------------------------------------------------------------------------------------------------------#

# Format data


#-------------------------------------------------------------------------------------------------------------------------#

# Generate figure




dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)}
setwd(starting_dir)
rm(list=ls())
gc()
