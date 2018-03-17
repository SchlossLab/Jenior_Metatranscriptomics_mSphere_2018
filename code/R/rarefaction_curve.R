
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')
library('rtk')

# Calculate rarefaction curves
rarefaction <- function(read_abundances){
  read_abundances <- subset(read_abundances, read_abundances[,1] != 0)[,1]
  richness <- round(as.numeric(rarefy(read_abundances, sample=(sum(read_abundances)-1))))
  rarefy_vect <- c()
  for (x in 1:richness) {
    temp <- rrarefy(read_abundances, sample=x)
    rarefy_vect[x] <- sum(temp != 0)
  }
  return(rarefy_vact)
}

# Read and format data
reForm <- function(read_files){
  readAbund <- read.delim(read_files[1], sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
  for (x in c(2:length(read_files))) {
    y <- read.delim(read_files[x], sep='\t', header=FALSE, row.names=1, na.strings=c('','NA'))
    readAbund <- merge(readAbund, y, by='row.names', all=TRUE)
    rownames(readAbund) <- readAbund$Row.names
    readAbund$Row.names <- NULL
  }

  readAbund <-  readAbund[rowSums(readAbund) != 0, ]
  
  return(readAbund)
}

# Define files
# Metagenomes
cef_630_metagenome <- 'data/read_mapping/metagenome/cefoperazone_630.Cefoperazone.metaG.final.pool.norm.txt'
cef_mock_metagenome <- 'data/read_mapping/metagenome/cefoperazone_mock.Cefoperazone.metaG.final.pool.norm.txt'
clinda_630_metagenome <- 'data/read_mapping/metagenome/clindamycin_630.Clindamycin.metaG.final.pool.norm.txt'
clinda_mock_metagenome <- 'data/read_mapping/metagenome/clindamycin_mock.Clindamycin.metaG.final.pool.norm.txt'
strep_630_metagenome <- 'data/read_mapping/metagenome/streptomycin_630.Streptomycin.metaG.final.pool.norm.txt'
strep_mock_metagenome <- 'data/read_mapping/metagenome/streptomycin_mock.Streptomycin.metaG.final.pool.norm.txt'
noabx_mock_metagenome <- 'data/read_mapping/metagenome/conventional_mock.Conventional.metaG.final.pool.norm.txt'
metaG <- c(cef_630_metagenome,cef_mock_metagenome,clinda_630_metagenome,clinda_mock_metagenome,
           strep_630_metagenome,strep_mock_metagenome,noabx_mock_metagenome)
rm(cef_630_metagenome,cef_mock_metagenome,clinda_630_metagenome,clinda_mock_metagenome,
   strep_630_metagenome,strep_mock_metagenome,noabx_mock_metagenome)
# Metatranscriptomes
cef_630_metatranscriptome <- 'data/read_mapping/metatranscriptome/cefoperazone_630.Cefoperazone.metaT.final.pool.norm.txt'
cef_mock_metatranscriptome <- 'data/read_mapping/metatranscriptome/cefoperazone_mock.Cefoperazone.metaT.final.pool.norm.txt'
clinda_630_metatranscriptome <- 'data/read_mapping/metatranscriptome/clindamycin_630.Clindamycin.metaT.final.pool.norm.txt'
clinda_mock_metatranscriptome <- 'data/read_mapping/metatranscriptome/clindamycin_mock.Clindamycin.metaT.final.pool.norm.txt'
strep_630_metatranscriptome <- 'data/read_mapping/metatranscriptome/streptomycin_630.Streptomycin.metaT.final.pool.norm.txt'
strep_mock_metatranscriptome <- 'data/read_mapping/metatranscriptome/streptomycin_mock.Streptomycin.metaT.final.pool.norm.txt'
noabx_mock_metatranscriptome <- 'data/read_mapping/metatranscriptome/conventional.Conventional.metaT.final.pool.norm.txt'
metaT <- c(cef_630_metatranscriptome,cef_mock_metatranscriptome,clinda_630_metatranscriptome,clinda_mock_metatranscriptome,
           strep_630_metatranscriptome,strep_mock_metatranscriptome,noabx_mock_metatranscriptome)
rm(cef_630_metatranscriptome,cef_mock_metatranscriptome,clinda_630_metatranscriptome,clinda_mock_metatranscriptome,
   strep_630_metatranscriptome,strep_mock_metatranscriptome,noabx_mock_metatranscriptome)
#-------------------------------------------------------------------------------------------------------------------------#

# Read in data

# Metagenomes
metagenome <- reForm(metaG)
metatranscriptome <- reForm(metaT)




cef_mock_metagenome <- reForm(cef_mock_metagenome)
clinda_630_metagenome <- reForm(clinda_630_metagenome)
clinda_mock_metagenome <- reForm(clinda_mock_metagenome)
strep_630_metagenome <- reForm(strep_630_metagenome)
strep_mock_metagenome <- reForm(strep_mock_metagenome)
noabx_mock_metagenome <- reForm(noabx_mock_metagenome)
# Metatranscriptomes
cef_630_metatranscriptome <- reForm(cef_630_metatranscriptome)
cef_mock_metatranscriptome <- reForm(cef_mock_metatranscriptome)
clinda_630_metatranscriptome <- reForm(clinda_630_metatranscriptome)
clinda_mock_metatranscriptome <- reForm(clinda_mock_metatranscriptome)
strep_630_metatranscriptome <- reForm(strep_630_metatranscriptome)
strep_mock_metatranscriptome <- reForm(strep_mock_metatranscriptome)
noabx_mock_metatranscriptome <- reForm(noabx_mock_metatranscriptome)

#-------------------------------------------------------------------------------------------------------------------------#


# Probably need to merge MetaG and MetaT data respectively to work with rtk




rarecurve(cef_630_metagenome)



cef_630_metagenome <- rarefaction(cef_630_metagenome)




cef_mock_metagenome <- rarefaction(cef_mock_metagenome)








# Clean up
setwd(starting_dir)
rm(list=ls())
