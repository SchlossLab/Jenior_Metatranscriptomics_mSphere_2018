# Set up environment
rm(list=ls())
gc()

# Load dependencies
deps <- c('vegan', 'plotrix', 'reshape2', 'GMD')
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Set seed for RNG
set.seed(9861)

# Load in functions
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Output plot name
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/figures/figure_2.pdf'

# Input 0.03 OTU shared file
shared_otu_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'

# Input 0.03 OTU taxonomy file
tax_otu_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/16S_analysis/formatted.all_treatments.0.03.cons.taxonomy'

# Input Metabolomes
metabolome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metabolome/metabolomics.tsv'

# Input Metadata
metadata_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metadata.tsv'

#----------------#

# Read in data

# 16S data
shared_otu <- read.delim(shared_otu_file, sep='\t', header=TRUE, row.names=2)
tax_otu <- read.delim(tax_otu_file, sep='\t', header=TRUE, row.names=1)
rm(shared_otu_file, tax_otu_file)

# Metabolomes
metabolome <- read.delim(metabolome_file, sep='\t', header=TRUE)
rm(metabolome_file)

# Metadata
metadata <- read.delim(metadata_file, sep='\t', header=T, row.names=1)
rm(metadata_file)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data

# Metadata
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2'), ] # Remove contaminated sample
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL
metadata$susceptibility <- NULL

# Metabolomes
metabolome <- metabolome[,!colnames(metabolome) %in% c('CefC5M2')] # Contaminated sample
metabolome <- metabolome[,!colnames(metabolome) %in% c('GfC1M1','GfC1M2','GfC1M3',
                                                      'GfC2M1','GfC2M2','GfC2M3',
                                                      'GfC3M1','GfC3M2','GfC3M3',
                                                      'GfC4M1','GfC4M2','GfC4M3',
                                                      'GfC5M1','GfC5M2','GfC5M3',
                                                      'GfC6M1','GfC6M2','GfC6M3')] # Germfree samples
#conv_metabolome <- metabolome[,colnames(metabolome) %in% c('ConvC1M1','ConvC1M2','ConvC1M3','ConvC1M4',
#                                                       'ConvC2M1','ConvC2M2','ConvC2M3','ConvC2M4','ConvC2M5')] # Untreated SPF samples
metabolome <- metabolome[,!colnames(metabolome) %in% c('ConvC1M1','ConvC1M2','ConvC1M3','ConvC1M4',
                                                       'ConvC2M1','ConvC2M2','ConvC2M3','ConvC2M4','ConvC2M5')] # Untreated SPF samples
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome_annotation <- metabolome[,1:4]
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
#bile_metabolome <- rbind(subset(metabolome, SUB_PATHWAY == 'Primary_Bile_Acid_Metabolism'),
#                         subset(metabolome, SUB_PATHWAY == 'Secondary_Bile_Acid_Metabolism'))
#bile_metabolome$SUB_PATHWAY <- NULL
metabolome$SUB_PATHWAY <- NULL
#carb_metabolome <- subset(metabolome, SUPER_PATHWAY == 'Carbohydrate')
#carb_metabolome$SUPER_PATHWAY <- NULL
#carb_metabolome <- t(carb_metabolome)
#aa_metabolome <- subset(metabolome, SUPER_PATHWAY == 'Amino_Acid')
#aa_metabolome$SUPER_PATHWAY <- NULL
#aa_metabolome <- t(aa_metabolome)
metabolome$SUPER_PATHWAY <- NULL
metabolome <- t(metabolome)

# 16S data
shared_otu$label <- NULL
shared_otu$numOtus <- NULL
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('CefC5M2'), ] # Contaminated sample
shared_otu <- shared_otu[,!(names(shared_otu) == 'Otu0004')] # Remove residual C. difficile OTU
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('StrepC4M1','StrepC4M2','StrepC4M3'), ]
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('GfC1M1','GfC1M2','GfC1M3',
                                                      'GfC2M1','GfC2M2','GfC2M3',
                                                      'GfC3M1','GfC3M2','GfC3M3',
                                                      'GfC4M1','GfC4M2','GfC4M3',
                                                      'GfC5M1','GfC5M2','GfC5M3',
                                                      'GfC6M1','GfC6M2','GfC6M3'), ] # Germfree samples
conv_otu <- shared_otu[rownames(shared_otu) %in% c('ConvC1M1','ConvC1M2','ConvC1M3','ConvC1M4',
                                                    'ConvC2M1','ConvC2M2','ConvC2M3','ConvC2M4','ConvC2M5'), ]
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('ConvC1M1','ConvC1M2','ConvC1M3','ConvC1M4',
                                                       'ConvC2M1','ConvC2M2','ConvC2M3','ConvC2M4','ConvC2M5'), ] # Untreated SPF samples

# Format taxonomy
tax_otu$OTU_short <- NULL

#-------------------------------------------------------------------------------------------------------------------------#

# Calculate axes

# Metabolome
metabolome_nmds <- metaMDS(metabolome, k=2, trymax=100)$points
metabolome_nmds[,1] <- metabolome_nmds[,1] - 0.05
metabolome_nmds[,2] <- metabolome_nmds[,2] * -1
metabolome_nmds <- clean_merge(metadata, metabolome_nmds)

# Community structure
otu_nmds <- metaMDS(shared_otu, k=2, trymax=100)$points
otu_nmds[,1] <- otu_nmds[,1] - 0.2
otu_nmds[,2] <- otu_nmds[,2] + 0.1
otu_nmds <- clean_merge(metadata, otu_nmds)

# Subset NMDS axes to color points
metabolome_cefoperazone <- subset(metabolome_nmds, abx == 'cefoperazone')
metabolome_cefoperazone_630 <- subset(metabolome_cefoperazone, infection == '630')
metabolome_cefoperazone_mock <- subset(metabolome_cefoperazone, infection == 'mock')
rm(metabolome_cefoperazone)
metabolome_clindamycin <- subset(metabolome_nmds, abx == 'clindamycin')
metabolome_clindamycin_630 <- subset(metabolome_clindamycin, infection == '630')
metabolome_clindamycin_mock <- subset(metabolome_clindamycin, infection == 'mock')
rm(metabolome_clindamycin)
metabolome_streptomycin <- subset(metabolome_nmds, abx == 'streptomycin')
metabolome_streptomycin_630 <- subset(metabolome_streptomycin, infection == '630')
metabolome_streptomycin_mock <- subset(metabolome_streptomycin, infection == 'mock')
rm(metabolome_streptomycin)
otu_cefoperazone <- subset(otu_nmds, abx == 'cefoperazone')
otu_cefoperazone_630 <- subset(otu_cefoperazone, infection == '630')
otu_cefoperazone_mock <- subset(otu_cefoperazone, infection == 'mock')
rm(otu_cefoperazone)
otu_clindamycin <- subset(otu_nmds, abx == 'clindamycin')
otu_clindamycin_630 <- subset(otu_clindamycin, infection == '630')
otu_clindamycin_mock <- subset(otu_clindamycin, infection == 'mock')
rm(otu_clindamycin)
otu_streptomycin <- subset(otu_nmds, abx == 'streptomycin')
otu_streptomycin_630 <- subset(otu_streptomycin, infection == '630')
otu_streptomycin_mock <- subset(otu_streptomycin, infection == 'mock')
rm(otu_streptomycin)

#----------------#

# Combine 16S and metabolome for heatmaps

# Reverse transformation
metabolome <- 10 ^ metabolome

# Merge with metadata
metabolome <- clean_merge(metadata, metabolome)
shared_otu <- clean_merge(metadata, shared_otu)

# Subset original tables
infected_metabolome <- subset(metabolome, infection == '630')
infected_metabolome$infection <- NULL
cef_infected_metabolome <- subset(infected_metabolome, abx == 'cefoperazone')
cef_infected_metabolome$abx <- NULL
clinda_infected_metabolome <- subset(infected_metabolome, abx == 'clindamycin')
clinda_infected_metabolome$abx <- NULL
strep_infected_metabolome <- subset(infected_metabolome, abx == 'streptomycin')
strep_infected_metabolome$abx <- NULL
mock_metabolome <- subset(metabolome, infection == 'mock')
mock_metabolome$infection <- NULL
cef_mock_metabolome <- subset(mock_metabolome, abx == 'cefoperazone')
cef_mock_metabolome$abx <- NULL
clinda_mock_metabolome <- subset(mock_metabolome, abx == 'clindamycin')
clinda_mock_metabolome$abx <- NULL
strep_mock_metabolome <- subset(mock_metabolome, abx == 'streptomycin')
strep_mock_metabolome$abx <- NULL
rm(metabolome, infected_metabolome, mock_metabolome)
infected_otu <- subset(shared_otu, infection == '630')
infected_otu$infection <- NULL
cef_infected_otu <- subset(infected_otu, abx == 'cefoperazone')
cef_infected_otu$abx <- NULL
clinda_infected_otu <- subset(infected_otu, abx == 'clindamycin')
clinda_infected_otu$abx <- NULL
strep_infected_otu <- subset(infected_otu, abx == 'streptomycin')
strep_infected_otu$abx <- NULL
mock_otu <- subset(shared_otu, infection == 'mock')
mock_otu$infection <- NULL
cef_mock_otu <- subset(mock_otu, abx == 'cefoperazone')
cef_mock_otu$abx <- NULL
clinda_mock_otu <- subset(mock_otu, abx == 'clindamycin')
clinda_mock_otu$abx <- NULL
strep_mock_otu <- subset(mock_otu, abx == 'streptomycin')
strep_mock_otu$abx <- NULL
rm(shared_otu, infected_otu, mock_otu)

# Calculate change in datasets using median of mock-infected
delta_cef_metabolome <- sweep(cef_infected_metabolome, 2, apply(cef_mock_metabolome, 2, median))
delta_cef_metabolome <- delta_cef_metabolome[,apply(delta_cef_metabolome, 2, sum) > 0.0]
delta_clinda_metabolome <- sweep(clinda_infected_metabolome, 2, apply(clinda_mock_metabolome, 2, median))
delta_clinda_metabolome <- delta_clinda_metabolome[,apply(delta_clinda_metabolome, 2, sum) > 0.0]
delta_strep_metabolome <- sweep(strep_infected_metabolome, 2, apply(strep_mock_metabolome, 2, median))
delta_strep_metabolome <- delta_strep_metabolome[,apply(delta_strep_metabolome, 2, sum) > 0.0]
delta_cef_otu <- sweep(cef_infected_otu, 2, apply(cef_mock_otu, 2, median))
delta_cef_otu <- delta_cef_otu[,apply(delta_cef_otu, 2, sum) > 0.0]
delta_clinda_otu <- sweep(clinda_infected_otu, 2, apply(clinda_mock_otu, 2, median))
delta_clinda_otu <- delta_clinda_otu[,apply(delta_clinda_otu, 2, sum) > 0.0]
delta_strep_otu <- sweep(strep_infected_otu, 2, apply(strep_mock_otu, 2, median))
delta_strep_otu <- delta_strep_otu[,apply(delta_strep_otu, 2, sum) > 0.0]
rm(cef_infected_metabolome, cef_mock_metabolome, clinda_infected_metabolome, clinda_mock_metabolome, strep_infected_metabolome, strep_mock_metabolome,
   cef_infected_otu, cef_mock_otu, clinda_infected_otu, clinda_mock_otu, strep_infected_otu, strep_mock_otu)

# Correlate differences in 16S with metabolome
cef_corr <- cor(x=delta_cef_metabolome, y=delta_cef_otu, method='spearman')
clinda_corr <- cor(x=delta_clinda_metabolome, y=delta_clinda_otu, method='spearman')
strep_corr <- cor(x=delta_strep_metabolome, y=delta_strep_otu, method='spearman')
rm(delta_cef_metabolome, delta_cef_otu,
   delta_clinda_metabolome, delta_clinda_otu,
   delta_strep_metabolome, delta_strep_otu)

# Subset correlations for clustering and sort by decreasing row sum - metabolome
strep_corr <- clean_merge(metabolome_annotation, strep_corr)
strep_corr$SUB_PATHWAY <- NULL
strep_corr$PUBCHEM <- NULL
strep_corr$KEGG <- NULL
strep_corr_aa <- strep_corr[strep_corr$SUPER_PATHWAY %in% c('Amino_Acid','Peptide'), ]
strep_corr_aa$SUPER_PATHWAY <- NULL
strep_corr_aa <- strep_corr_aa[,order(-colSums(strep_corr_aa))]
strep_corr_carb <- subset(strep_corr, SUPER_PATHWAY == 'Carbohydrate')
strep_corr_carb$SUPER_PATHWAY <- NULL
strep_corr_carb <- strep_corr_carb[,order(-colSums(strep_corr_carb))]
strep_corr_vit <- subset(strep_corr, SUPER_PATHWAY == 'Cofactors_and_Vitamins')
strep_corr_vit$SUPER_PATHWAY <- NULL
strep_corr_vit <- strep_corr_vit[,order(-colSums(strep_corr_vit))]
strep_corr_en <- subset(strep_corr, SUPER_PATHWAY == 'Energy')
strep_corr_en$SUPER_PATHWAY <- NULL
strep_corr_en <- strep_corr_en[,order(-colSums(strep_corr_en))]
strep_corr_lip <- subset(strep_corr, SUPER_PATHWAY == 'Lipid')
strep_corr_lip$SUPER_PATHWAY <- NULL
strep_corr_lip <- strep_corr_lip[,order(-colSums(strep_corr_lip))]
strep_corr_nuc <- subset(strep_corr, SUPER_PATHWAY == 'Nucleotide')
strep_corr_nuc$SUPER_PATHWAY <- NULL
strep_corr_nuc <- strep_corr_nuc[,order(-colSums(strep_corr_nuc))]
strep_corr_xen <- subset(strep_corr, SUPER_PATHWAY == 'Xenobiotics')
strep_corr_xen$SUPER_PATHWAY <- NULL
strep_corr_xen <- strep_corr_xen[,order(-colSums(strep_corr_xen))]
strep_corr <- rbind(strep_corr_aa, strep_corr_carb, strep_corr_vit, strep_corr_en, strep_corr_lip, strep_corr_nuc, strep_corr_xen)
rm(strep_corr_aa, strep_corr_carb, strep_corr_vit, strep_corr_en, strep_corr_lip, strep_corr_nuc, strep_corr_xen)
cef_corr <- clean_merge(metabolome_annotation, cef_corr)
cef_corr$SUB_PATHWAY <- NULL
cef_corr$PUBCHEM <- NULL
cef_corr$KEGG <- NULL
cef_corr_aa <- cef_corr[cef_corr$SUPER_PATHWAY %in% c('Amino_Acid','Peptide'), ]
cef_corr_aa$SUPER_PATHWAY <- NULL
cef_corr_aa <- cef_corr_aa[,order(-colSums(cef_corr_aa))]
cef_corr_carb <- subset(cef_corr, SUPER_PATHWAY == 'Carbohydrate')
cef_corr_carb$SUPER_PATHWAY <- NULL
cef_corr_carb <- cef_corr_carb[,order(-colSums(cef_corr_carb))]
cef_corr_vit <- subset(cef_corr, SUPER_PATHWAY == 'Cofactors_and_Vitamins')
cef_corr_vit$SUPER_PATHWAY <- NULL
cef_corr_vit <- cef_corr_vit[,order(-colSums(cef_corr_vit))]
cef_corr_en <- subset(cef_corr, SUPER_PATHWAY == 'Energy')
cef_corr_en$SUPER_PATHWAY <- NULL
cef_corr_en <- cef_corr_en[,order(-colSums(cef_corr_en))]
cef_corr_lip <- subset(cef_corr, SUPER_PATHWAY == 'Lipid')
cef_corr_lip$SUPER_PATHWAY <- NULL
cef_corr_lip <- cef_corr_lip[,order(-colSums(cef_corr_lip))]
cef_corr_nuc <- subset(cef_corr, SUPER_PATHWAY == 'Nucleotide')
cef_corr_nuc$SUPER_PATHWAY <- NULL
cef_corr_nuc <- cef_corr_nuc[,order(-colSums(cef_corr_nuc))]
cef_corr_xen <- subset(cef_corr, SUPER_PATHWAY == 'Xenobiotics')
cef_corr_xen$SUPER_PATHWAY <- NULL
cef_corr_xen <- cef_corr_xen[,order(-colSums(cef_corr_xen))]
cef_corr <- rbind(cef_corr_aa, cef_corr_carb, cef_corr_vit, cef_corr_en, cef_corr_lip, cef_corr_nuc, cef_corr_xen)
rm(cef_corr_aa, cef_corr_carb, cef_corr_vit, cef_corr_en, cef_corr_lip, cef_corr_nuc, cef_corr_xen)
clinda_corr <- clean_merge(metabolome_annotation, clinda_corr)
clinda_corr$SUB_PATHWAY <- NULL
clinda_corr$PUBCHEM <- NULL
clinda_corr$KEGG <- NULL
clinda_corr_aa <- clinda_corr[clinda_corr$SUPER_PATHWAY %in% c('Amino_Acid','Peptide'), ]
clinda_corr_aa$SUPER_PATHWAY <- NULL
clinda_corr_aa <- clinda_corr_aa[,order(-colSums(clinda_corr_aa))]
clinda_corr_carb <- subset(clinda_corr, SUPER_PATHWAY == 'Carbohydrate')
clinda_corr_carb$SUPER_PATHWAY <- NULL
clinda_corr_carb <- clinda_corr_carb[,order(-colSums(clinda_corr_carb))]
clinda_corr_vit <- subset(clinda_corr, SUPER_PATHWAY == 'Cofactors_and_Vitamins')
clinda_corr_vit$SUPER_PATHWAY <- NULL
clinda_corr_vit <- clinda_corr_vit[,order(-colSums(clinda_corr_vit))]
clinda_corr_en <- subset(clinda_corr, SUPER_PATHWAY == 'Energy')
clinda_corr_en$SUPER_PATHWAY <- NULL
clinda_corr_en <- clinda_corr_en[,order(-colSums(clinda_corr_en))]
clinda_corr_lip <- subset(clinda_corr, SUPER_PATHWAY == 'Lipid')
clinda_corr_lip$SUPER_PATHWAY <- NULL
clinda_corr_lip <- clinda_corr_lip[,order(-colSums(clinda_corr_lip))]
clinda_corr_nuc <- subset(clinda_corr, SUPER_PATHWAY == 'Nucleotide')
clinda_corr_nuc$SUPER_PATHWAY <- NULL
clinda_corr_nuc <- clinda_corr_nuc[,order(-colSums(clinda_corr_nuc))]
clinda_corr_xen <- subset(clinda_corr, SUPER_PATHWAY == 'Xenobiotics')
clinda_corr_xen$SUPER_PATHWAY <- NULL
clinda_corr_xen <- clinda_corr_xen[,order(-colSums(clinda_corr_xen))]
clinda_corr <- rbind(clinda_corr_aa, clinda_corr_carb, clinda_corr_vit, clinda_corr_en, clinda_corr_lip, clinda_corr_nuc, clinda_corr_xen)
rm(clinda_corr_aa, clinda_corr_carb, clinda_corr_vit, clinda_corr_en, clinda_corr_lip, clinda_corr_nuc, clinda_corr_xen)

# Subset correlations for clustering and sort by decreasing row sum - 16S
strep_corr <- clean_merge(tax_otu, t(strep_corr))
strep_corr$genus <- NULL
strep_corr_bact <- subset(strep_corr, phylum == 'Bacteroidetes')
strep_corr_bact$phylum <- NULL
strep_corr_bact <- strep_corr_bact[,order(-colSums(strep_corr_bact))]
strep_corr_firm <- subset(strep_corr, phylum == 'Firmicutes')
strep_corr_firm$phylum <- NULL
strep_corr_firm <- strep_corr_firm[,order(-colSums(strep_corr_firm))]
strep_corr_prot <- subset(strep_corr, phylum == 'Proteobacteria')
strep_corr_prot$phylum <- NULL
strep_corr_prot <- strep_corr_prot[,order(-colSums(strep_corr_prot))]
strep_corr_actino <- subset(strep_corr, phylum == 'Actinobacteria')
strep_corr_actino$phylum <- NULL
strep_corr_actino <- strep_corr_actino[,order(-colSums(strep_corr_actino))]
strep_corr_other <- strep_corr[strep_corr$phylum %in% c('Bacteria_unclassified','Fusobacteria','Synergistetes'), ]
strep_corr_other$phylum <- NULL
strep_corr_other <- strep_corr_other[,order(-colSums(strep_corr_other))]
strep_corr <- rbind(strep_corr_bact, strep_corr_firm, strep_corr_prot, strep_corr_actino, strep_corr_other)
strep_corr <- t(strep_corr)
rm(strep_corr_bact, strep_corr_firm, strep_corr_prot, strep_corr_actino, strep_corr_other)
cef_corr <- clean_merge(tax_otu, t(cef_corr))
cef_corr$genus <- NULL
cef_corr_bact <- subset(cef_corr, phylum == 'Bacteroidetes')
cef_corr_bact$phylum <- NULL
cef_corr_bact <- cef_corr_bact[,order(-colSums(cef_corr_bact))]
cef_corr_firm <- subset(cef_corr, phylum == 'Firmicutes')
cef_corr_firm$phylum <- NULL
cef_corr_firm <- cef_corr_firm[,order(-colSums(cef_corr_firm))]
cef_corr_prot <- subset(cef_corr, phylum == 'Proteobacteria')
cef_corr_prot$phylum <- NULL
cef_corr_prot <- cef_corr_prot[,order(-colSums(cef_corr_prot))]
cef_corr_actino <- subset(cef_corr, phylum == 'Actinobacteria')
cef_corr_actino$phylum <- NULL
cef_corr_actino <- cef_corr_actino[,order(-colSums(cef_corr_actino))]
cef_corr_other <- cef_corr[cef_corr$phylum %in% c('Bacteria_unclassified','Fusobacteria','Synergistetes'), ]
cef_corr_other$phylum <- NULL
cef_corr_other <- cef_corr_other[,order(-colSums(cef_corr_other))]
cef_corr <- rbind(cef_corr_bact, cef_corr_firm, cef_corr_prot, cef_corr_actino, cef_corr_other)
cef_corr <- t(cef_corr)
rm(cef_corr_bact, cef_corr_firm, cef_corr_prot, cef_corr_actino, cef_corr_other)
clinda_corr <- clean_merge(tax_otu, t(clinda_corr))
clinda_corr$genus <- NULL
clinda_corr_bact <- subset(clinda_corr, phylum == 'Bacteroidetes')
clinda_corr_bact$phylum <- NULL
clinda_corr_bact <- clinda_corr_bact[,order(-colSums(clinda_corr_bact))]
clinda_corr_firm <- subset(clinda_corr, phylum == 'Firmicutes')
clinda_corr_firm$phylum <- NULL
clinda_corr_firm <- clinda_corr_firm[,order(-colSums(clinda_corr_firm))]
clinda_corr_prot <- subset(clinda_corr, phylum == 'Proteobacteria')
clinda_corr_prot$phylum <- NULL
clinda_corr_prot <- clinda_corr_prot[,order(-colSums(clinda_corr_prot))]
clinda_corr_actino <- subset(clinda_corr, phylum == 'Actinobacteria')
clinda_corr_actino$phylum <- NULL
clinda_corr_actino <- clinda_corr_actino[,order(-colSums(clinda_corr_actino))]
clinda_corr_other <- clinda_corr[clinda_corr$phylum %in% c('Bacteria_unclassified','Fusobacteria','Synergistetes'), ]
clinda_corr_other$phylum <- NULL
clinda_corr_other <- clinda_corr_other[,order(-colSums(clinda_corr_other))]
clinda_corr <- rbind(clinda_corr_bact, clinda_corr_firm, clinda_corr_prot, clinda_corr_actino, clinda_corr_other)
clinda_corr <- t(clinda_corr)
rm(clinda_corr_bact, clinda_corr_firm, clinda_corr_prot, clinda_corr_actino, clinda_corr_other)

#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_file, width=8.5, height=14)
layout(matrix(c(1,2,
                3,3,
                4,4,
                5,5),
              nrow=4, ncol=2, byrow=TRUE))

par(mar=c(4,4,1,1), las=1, mgp=c(3,0.75,0), xaxs='i', yaxs='i')

#-------------------#

# OTUs alone
plot(x=otu_nmds$MDS1, y=otu_nmds$MDS2, xlim=c(-1.2,1.2), ylim=c(-1.5,1.5),
     xlab='NMDS axis 1', ylab='NMDS axis 2', xaxt='n', yaxt='n', pch=19, cex=0.2)
axis(side=1, at=seq(-1.2,1.2,0.4), labels=c(-1.2,-0.8,-0.4,0,0.4,0.8,1.2))
axis(side=2, at=seq(-1.5,1.5,0.5), labels=seq(-1.5,1.5,0.5))
mtext('a', side=2, line=2, las=2, adj=2, padj=-10, cex=1.2, font=2)
legend('topright', legend=c('16S rRNA Gene Sequencing'), pch=21, cex=1.4, pt.cex=0, bty='n')
legend('bottomright', legend=c('Streptomycin-pretreated','Cefoperzone-pretreated','Clindamycin-pretreated'), 
       pt.bg=c(strep_col,cef_col,clinda_col), 
       pch=22, cex=1.2, pt.cex=2.4)
legend('bottomleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(16,17), cex=1.2, pt.cex=c(2.2,2))

points(x=otu_cefoperazone_630$MDS1, y=otu_cefoperazone_630$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=otu_clindamycin_630$MDS1, y=otu_clindamycin_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=otu_streptomycin_630$MDS1, y=otu_streptomycin_630$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=otu_cefoperazone_mock$MDS1, y=otu_cefoperazone_mock$MDS2, bg=cef_col, pch=24, cex=1.8, lwd=1.2)
points(x=otu_clindamycin_mock$MDS1, y=otu_clindamycin_mock$MDS2, bg=clinda_col, pch=24, cex=1.8, lwd=1.2)
points(x=otu_streptomycin_mock$MDS1, y=otu_streptomycin_mock$MDS2, bg=strep_col, pch=24, cex=1.8, lwd=1.2)

#-------------------#

# Metabolomics alone
plot(x=metabolome_nmds$MDS1, y=metabolome_nmds$MDS2, xlim=c(-0.25,0.2), ylim=c(-0.25,0.2),
     xlab='NMDS axis 1', ylab='NMDS axis 2', xaxt='n', yaxt='n', pch=19, cex=0.2)
axis(side=1, at=seq(-0.25,0.2,0.05), labels=seq(-0.25,0.2,0.05))
axis(side=2, at=seq(-0.25,0.2,0.05), labels=seq(-0.25,0.2,0.05))
mtext('b', side=2, line=2, las=2, adj=2, padj=-10, cex=1.2, font=2)
legend('topright', legend=c('Untargeted Metabolomics'), pch=21, cex=1.4, pt.cex=0, bty='n')
legend('bottomright', legend=c('Streptomycin-pretreated','Cefoperzone-pretreated','Clindamycin-pretreated'), 
       pt.bg=c(strep_col,cef_col,clinda_col), 
       pch=22, cex=1.2, pt.cex=2.4)
legend('bottomleft', legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(16,17), cex=1.2, pt.cex=c(2.2,2))

points(x=metabolome_cefoperazone_630$MDS1, y=metabolome_cefoperazone_630$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_clindamycin_630$MDS1, y=metabolome_clindamycin_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_streptomycin_630$MDS1, y=metabolome_streptomycin_630$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=metabolome_cefoperazone_mock$MDS1, y=metabolome_cefoperazone_mock$MDS2, bg=cef_col, pch=24, cex=1.8, lwd=1.2)
points(x=metabolome_clindamycin_mock$MDS1, y=metabolome_clindamycin_mock$MDS2, bg=clinda_col, pch=24, cex=1.8, lwd=1.2)
points(x=metabolome_streptomycin_mock$MDS1, y=metabolome_streptomycin_mock$MDS2, bg=strep_col, pch=24, cex=1.8, lwd=1.2)

#-------------------#

par(mar=c(1,1,1,1), las=1, mgp=c(3,0.75,0), xaxs='i', yaxs='i')

# Streptomycin heatmap
heatmap.3(t(strep_corr), dendrogram='none', key=FALSE, main='', labRow='', labCol='', cluster.by.row=FALSE, cluster.by.col=FALSE)
mtext('c', side=2, line=2, las=2, adj=2, padj=-9, cex=1.2, font=2)

#-------------------#

# Cefoperazone heatmap
heatmap.3(t(cef_corr), dendrogram='none', key=FALSE, main='', labRow='', labCol='', cluster.by.row=FALSE, cluster.by.col=FALSE)
mtext('d', side=2, line=2, las=2, adj=2, padj=-9, cex=1.2, font=2)

#-------------------#

# Clindamycin heatmap
heatmap.3(t(clinda_corr), dendrogram='none', key=FALSE, main='', labRow='', labCol='', cluster.by.row=FALSE, cluster.by.col=FALSE)
mtext('e', side=2, line=2, las=2, adj=2, padj=-9, cex=1.2, font=2)


dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
rm(list=ls())
gc()

