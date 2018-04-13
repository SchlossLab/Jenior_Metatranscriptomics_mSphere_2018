
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')

# Define files

# Metadata
metadata <- 'data/metadata.tsv'

# Normalized Metatranscriptomes
cef_normalized_reads <- 'data/read_mapping/cef_normalized_metaT.tsv'
clinda_normalized_reads <- 'data/read_mapping/clinda_normalized_metaT.tsv'
strep_normalized_reads <- 'data/read_mapping/strep_normalized_metaT.tsv'

# 16S abundances
genus_shared <- 'data/16S_analysis/all_treatments.0.03.unique_list.0.03.filter.0.03.subsample.shared'
genus_tax <- 'data/16S_analysis/all_treatments.genus.final.taxonomy'

# KEGG taxonomy IDs
kegg_tax <- 'data/kegg/kegg_taxonomy.tsv'

# Taxonomy colors
tax_colors <- 'data/taxonomy_color.tsv'

# Output plot
plot_file <- 'results/figures/figure_3.tiff'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data

# Metadata
metadata <- read.delim(metadata, sep='\t', header=TRUE, row.names=1)

# Normalized Metatranscriptomes
cef_normalized_reads <- read.delim(cef_normalized_reads, sep='\t', header=TRUE, row.names=7)
clinda_normalized_reads <- read.delim(clinda_normalized_reads, sep='\t', header=TRUE, row.names=7)
strep_normalized_reads <- read.delim(strep_normalized_reads, sep='\t', header=TRUE, row.names=7)
  
# KEGG organism file
kegg_tax <- read.delim(kegg_tax, sep='\t', header=TRUE)
kegg_tax[] <- lapply(kegg_tax, as.character)

# 16S data
genus_shared <- read.delim(genus_shared, sep='\t', header=TRUE, row.names=2)
genus_shared$label <- NULL
genus_shared$numOtus <- NULL
genus_tax <- read.delim(genus_tax, sep='\t', header=TRUE)
rownames(genus_tax) <- genus_tax$OTU
genus_tax$OTU <- NULL
 
# Taxonomy colors
tax_colors <- read.delim(tax_colors, sep='\t', header=TRUE)
tax_colors[] <- lapply(tax_colors, as.character)

#-------------------------------------------------------------------------------------------------------------------------#

# Format transcription
cef_normalized_reads$gene <- NULL
cef_normalized_reads[,c(1:2)] <- log2(cef_normalized_reads[,c(1:2)] + 1)
clinda_normalized_reads$gene <- NULL
clinda_normalized_reads[,c(1:2)] <- log2(clinda_normalized_reads[,c(1:2)] + 1)
strep_normalized_reads$gene <- NULL
strep_normalized_reads[,c(1:2)] <- log2(strep_normalized_reads[,c(1:2)] + 1)

# Remove those genes without an organism annotation
cef_annotated <- cef_normalized_reads[!is.na(cef_normalized_reads$organism),]
clinda_annotated <- clinda_normalized_reads[!is.na(clinda_normalized_reads$organism),]
strep_annotated <- strep_normalized_reads[!is.na(strep_normalized_reads$organism),]

# Remove any C. difficile genes with transcription only in infected
cdiff_omit <- c('cdf','pdc','cdc','cdl','pdf')
cef_annotated <- subset(cef_annotated, !(cef_annotated$organism %in% cdiff_omit))
clinda_annotated <- subset(clinda_annotated, !(clinda_annotated$organism %in% cdiff_omit))
strep_annotated <- subset(strep_annotated, !(strep_annotated$organism %in% cdiff_omit))
rm(cdiff_omit)

# Remove all possible mammalian genes
mamm_omit <- c('fab','cfa','ggo','hgl','hsa','mcc','mdo','pon','aml',
               'ptr','rno','shr','ssc','aml','bta','cge','ecb',
               'pps','fca','mmu','oaa','gga','ola','acs','aga')
strep_annotated <- subset(strep_annotated, !(strep_annotated$organism %in% mamm_omit))
cef_annotated <- subset(cef_annotated, !(cef_annotated$organism %in% mamm_omit))
clinda_annotated <- subset(clinda_annotated, !(clinda_annotated$organism %in% mamm_omit))
rm(mamm_omit)

# Calculate correlation coefficients
strep_corr <- as.character(round(cor.test(strep_annotated[,2], strep_annotated[,1], method='spearman', exact=FALSE)$estimate, digits=3))
cef_corr <- as.character(round(cor.test(cef_annotated[,2], cef_annotated[,1], method='spearman', exact=FALSE)$estimate, digits=3))
clinda_corr <- as.character(round(cor.test(clinda_annotated[,2], clinda_annotated[,1], method='spearman', exact=FALSE)$estimate, digits=3))

# Collate transcript abundances at genus-level 
genus_shared <- t(genus_shared)
genus_tax_shared <- merge(genus_tax, genus_shared, by='row.names')
genus_tax_shared$Row.names <- NULL
genus_tax_agg_shared <- aggregate(. ~ Genus, data=genus_tax_shared, FUN=sum)
rownames(genus_tax_agg_shared) <- genus_tax_agg_shared$Genus
genus_tax_agg_shared$Genus <- NULL
genus_tax_agg_shared <- t(genus_tax_agg_shared)
rm(genus_tax, genus_shared, genus_tax_shared)

# Combine with metadata and collate within treatment
genus_tax_agg_metadata_shared <- merge(metadata, genus_tax_agg_shared, by='row.names')
genus_tax_agg_metadata_shared$Row.names <- NULL
genus_tax_agg_metadata_shared$cage <- NULL
genus_tax_agg_metadata_shared$mouse <- NULL
genus_tax_agg_metadata_shared$gender <- NULL
genus_tax_agg_metadata_shared$type <- NULL
genus_tax_agg_metadata_shared$susceptibility <- NULL
genus_tax_agg_metadata_shared$infection <- NULL
genus_tax_agg_metadata_shared <- subset(genus_tax_agg_metadata_shared, genus_tax_agg_metadata_shared$abx %in% c('cefoperazone','clindamycin','streptomycin'))
genus_tax_agg_metadata_shared <- aggregate(. ~ abx, data=genus_tax_agg_metadata_shared, FUN=median)
rownames(genus_tax_agg_metadata_shared) <- genus_tax_agg_metadata_shared$abx
genus_tax_agg_metadata_shared$abx <- NULL
genus_tax_agg_metadata_shared <- (genus_tax_agg_metadata_shared / rowSums(genus_tax_agg_metadata_shared)) * 100
genus_tax_agg_metadata_shared <- as.data.frame(t(genus_tax_agg_metadata_shared))
cef_genus <- genus_tax_agg_metadata_shared
cef_genus$clindamycin <- NULL
cef_genus$streptomycin <- NULL
clinda_genus <- genus_tax_agg_metadata_shared
clinda_genus$cefoperazone <- NULL
clinda_genus$streptomycin <- NULL
strep_genus <- genus_tax_agg_metadata_shared
strep_genus$clindamycin <- NULL
strep_genus$cefoperazone <- NULL
rm(genus_tax_agg_shared, genus_tax_agg_metadata_shared)

# Using previously defined lines, find outliers to y = x
strep_630_outliers <- subset(strep_annotated, strep_annotated$strep_630_metaT_reads > strep_annotated$strep_mock_metaT_reads + 2)
strep_mock_outliers <- subset(strep_annotated, strep_annotated$strep_mock_metaT_reads > strep_annotated$strep_630_metaT_reads + 2)
cef_630_outliers <- subset(cef_annotated, cef_annotated$cef_630_metaT_reads > cef_annotated$cef_mock_metaT_reads + 2)
cef_mock_outliers <- subset(cef_annotated, cef_annotated$cef_mock_metaT_reads > cef_annotated$cef_630_metaT_reads + 2)
clinda_630_outliers <- clinda_annotated[(clinda_annotated$clinda_630_metaT_reads > clinda_annotated$clinda_mock_metaT_reads + 2),]
clinda_mock_outliers <- clinda_annotated[(clinda_annotated$clinda_mock_metaT_reads > clinda_annotated$clinda_630_metaT_reads + 2),]

# Remove outliers from the rest of the genes
strep_annotated <- strep_annotated[!row.names(strep_annotated) %in% row.names(strep_630_outliers), ]
strep_annotated <- strep_annotated[!row.names(strep_annotated) %in% row.names(strep_mock_outliers), ]
cef_annotated <- cef_annotated[!row.names(cef_annotated) %in% row.names(cef_630_outliers), ]
cef_annotated <- cef_annotated[!row.names(cef_annotated) %in% row.names(cef_mock_outliers), ]
clinda_annotated <- clinda_annotated[!row.names(clinda_annotated) %in% row.names(clinda_630_outliers), ]
clinda_annotated <- clinda_annotated[!row.names(clinda_annotated) %in% row.names(clinda_mock_outliers), ]

# Calculate mean distance of outliers to x=y
dists_630 <- c()
for (i in 1:nrow(strep_630_outliers)){
  dists_630[i] <- dist_xy(c(strep_630_outliers$strep_630_metaT_reads[i], strep_630_outliers$strep_mock_metaT_reads[i]))
}
dists_mock <- c()
for (i in 1:nrow(strep_mock_outliers)){
  dists_mock[i] <- dist_xy(c(strep_mock_outliers$strep_630_metaT_reads[i], strep_mock_outliers$strep_mock_metaT_reads[i]))
}
round(mean(c(dists_630, dists_mock)), 3)
as.numeric(length(c(dists_630, dists_mock)))
dists_630 <- c()
for (i in 1:nrow(cef_630_outliers)){
  dists_630[i] <- dist_xy(c(cef_630_outliers$cef_630_metaT_reads[i], cef_630_outliers$cef_mock_metaT_reads[i]))
}
dists_mock <- c()
for (i in 1:nrow(cef_mock_outliers)){
  dists_mock[i] <- dist_xy(c(cef_mock_outliers$cef_630_metaT_reads[i], cef_mock_outliers$cef_mock_metaT_reads[i]))
}
round(mean(c(dists_630, dists_mock)), 3)
as.numeric(length(c(dists_630, dists_mock)))
dists_630 <- c()
for (i in 1:nrow(clinda_630_outliers)){
  dists_630[i] <- dist_xy(c(clinda_630_outliers$clinda_630_metaT_reads[i], clinda_630_outliers$clinda_mock_metaT_reads[i]))
}
dists_mock <- c()
for (i in 1:nrow(clinda_mock_outliers)){
  dists_mock[i] <- dist_xy(c(clinda_mock_outliers$clinda_630_metaT_reads[i], clinda_mock_outliers$clinda_mock_metaT_reads[i]))
}
round(mean(c(dists_630, dists_mock)), 3)
as.numeric(length(c(dists_630, dists_mock)))
rm(dists_630, dists_mock)

#-------------------------------------------------------------------------------------------------------------------------#

# Break down outliers into taxonomic groups

# Drop levels
strep_630_outliers[] <- lapply(strep_630_outliers, as.character)
strep_mock_outliers[] <- lapply(strep_mock_outliers, as.character)
cef_630_outliers[] <- lapply(cef_630_outliers, as.character)
cef_mock_outliers[] <- lapply(cef_mock_outliers, as.character)
clinda_630_outliers[] <- lapply(clinda_630_outliers, as.character)
clinda_mock_outliers[] <- lapply(clinda_mock_outliers, as.character)
strep_annotated[] <- lapply(strep_annotated, as.character)
cef_annotated[] <- lapply(cef_annotated, as.character)
clinda_annotated[] <- lapply(clinda_annotated, as.character)

# Save KEGG ID names
strep_630_outliers$kegg_id <- rownames(strep_630_outliers)
strep_mock_outliers$kegg_id <- rownames(strep_mock_outliers)
cef_630_outliers$kegg_id <- rownames(cef_630_outliers)
cef_mock_outliers$kegg_id <- rownames(cef_mock_outliers)
clinda_630_outliers$kegg_id <- rownames(clinda_630_outliers)
clinda_mock_outliers$kegg_id <- rownames(clinda_mock_outliers)
strep_annotated$kegg_id <- rownames(strep_annotated)
cef_annotated$kegg_id <- rownames(cef_annotated)
clinda_annotated$kegg_id <- rownames(clinda_annotated)

# Merge with KEGG taxonomy
strep_630_outliers <- merge(x=strep_630_outliers, y=kegg_tax, by.x='organism', by.y='org_code', all.x=TRUE)
strep_630_outliers$org_code <- NULL
strep_630_outliers$strep_630_metaT_reads <- as.numeric(strep_630_outliers$strep_630_metaT_reads)
strep_630_outliers$strep_mock_metaT_reads <- as.numeric(strep_630_outliers$strep_mock_metaT_reads)
strep_mock_outliers <- merge(x=strep_mock_outliers, y=kegg_tax, by.x='organism', by.y='org_code', all.x=TRUE)
strep_mock_outliers$org_code <- NULL
strep_mock_outliers$strep_630_metaT_reads <- as.numeric(strep_mock_outliers$strep_630_metaT_reads)
strep_mock_outliers$strep_mock_metaT_reads <- as.numeric(strep_mock_outliers$strep_mock_metaT_reads)
cef_630_outliers <- merge(x=cef_630_outliers, y=kegg_tax, by.x='organism', by.y='org_code', all.x=TRUE)
cef_630_outliers$org_code <- NULL
cef_630_outliers$cef_630_metaT_reads <- as.numeric(cef_630_outliers$cef_630_metaT_reads)
cef_630_outliers$cef_mock_metaT_reads <- as.numeric(cef_630_outliers$cef_mock_metaT_reads)
cef_mock_outliers <- merge(x=cef_mock_outliers, y=kegg_tax, by.x='organism', by.y='org_code', all.x=TRUE)
cef_mock_outliers$org_code <- NULL
cef_mock_outliers$cef_630_metaT_reads <- as.numeric(cef_mock_outliers$cef_630_metaT_reads)
cef_mock_outliers$cef_mock_metaT_reads <- as.numeric(cef_mock_outliers$cef_mock_metaT_reads)
clinda_630_outliers <- merge(x=clinda_630_outliers, y=kegg_tax, by.x='organism', by.y='org_code', all.x=TRUE)
clinda_630_outliers$org_code <- NULL
clinda_630_outliers$clinda_630_metaT_reads <- as.numeric(clinda_630_outliers$clinda_630_metaT_reads)
clinda_630_outliers$clinda_mock_metaT_reads <- as.numeric(clinda_630_outliers$clinda_mock_metaT_reads)
clinda_mock_outliers <- merge(x=clinda_mock_outliers, y=kegg_tax, by.x='organism', by.y='org_code', all.x=TRUE)
clinda_mock_outliers$org_code <- NULL
clinda_mock_outliers$clinda_630_metaT_reads <- as.numeric(clinda_mock_outliers$clinda_630_metaT_reads)
clinda_mock_outliers$clinda_mock_metaT_reads <- as.numeric(clinda_mock_outliers$clinda_mock_metaT_reads)
strep_annotated <- merge(x=strep_annotated, y=kegg_tax, by.x='organism', by.y='org_code', all.x=TRUE)
strep_annotated$org_code <- NULL
strep_annotated$strep_630_metaT_reads <- as.numeric(strep_annotated$strep_630_metaT_reads)
strep_annotated$strep_mock_metaT_reads <- as.numeric(strep_annotated$strep_mock_metaT_reads)
cef_annotated <- merge(x=cef_annotated, y=kegg_tax, by.x='organism', by.y='org_code', all.x=TRUE)
cef_annotated$org_code <- NULL
cef_annotated$cef_630_metaT_reads <- as.numeric(cef_annotated$cef_630_metaT_reads)
cef_annotated$cef_mock_metaT_reads <- as.numeric(cef_annotated$cef_mock_metaT_reads)
clinda_annotated <- merge(x=clinda_annotated, y=kegg_tax, by.x='organism', by.y='org_code', all.x=TRUE)
clinda_annotated$org_code <- NULL
clinda_annotated$clinda_630_metaT_reads <- as.numeric(clinda_annotated$clinda_630_metaT_reads)
clinda_annotated$clinda_mock_metaT_reads <- as.numeric(clinda_annotated$clinda_mock_metaT_reads)
rm(kegg_tax)

#-------------------------------------------------------------------------------------------------------------------------#

# Get points ready for plotting

# Define colors based on genus
strep_630_outliers <- merge(x=strep_630_outliers, y=tax_colors, by.x='genus', by.y='taxonomy', all.x=TRUE)
strep_630_outliers$color[is.na(strep_630_outliers$color)] <- 'white'
strep_mock_outliers <- merge(x=strep_mock_outliers, y=tax_colors, by.x='genus', by.y='taxonomy', all.x=TRUE)
strep_mock_outliers$color[is.na(strep_mock_outliers$color)] <- 'white'
cef_630_outliers <- merge(x=cef_630_outliers, y=tax_colors, by.x='genus', by.y='taxonomy', all.x=TRUE)
cef_630_outliers$color[is.na(cef_630_outliers$color)] <- 'white'
cef_mock_outliers <- merge(x=cef_mock_outliers, y=tax_colors, by.x='genus', by.y='taxonomy', all.x=TRUE)
cef_mock_outliers$color[is.na(cef_mock_outliers$color)] <- 'white'
clinda_630_outliers <- merge(x=clinda_630_outliers, y=tax_colors, by.x='genus', by.y='taxonomy', all.x=TRUE)
clinda_630_outliers$color[is.na(clinda_630_outliers$color)] <- 'white'
clinda_mock_outliers <- merge(x=clinda_mock_outliers, y=tax_colors, by.x='genus', by.y='taxonomy', all.x=TRUE)
clinda_mock_outliers$color[is.na(clinda_mock_outliers$color)] <- 'white'

# Find difference in expression for outliers
strep_outliers <- rbind(strep_630_outliers, strep_mock_outliers)
strep_difference <- abs(((2^strep_outliers$strep_mock_metaT_reads) - 1) - ((2^strep_outliers$strep_630_metaT_reads) - 1))
strep_difference <- as.data.frame(cbind(strep_outliers$genus, strep_difference))
colnames(strep_difference) <- c('genus', 'difference')
strep_gene_count <- as.data.frame(table(strep_difference$genus))
colnames(strep_gene_count) <- c('genus', 'count')
strep_difference <- aggregate(. ~ genus, data=strep_difference, FUN=sum)
strep_difference <- merge(strep_difference, strep_gene_count, by='genus')
strep_difference$difference <- strep_difference$difference / strep_difference$count
strep_difference$count <- NULL
rm(strep_gene_count)
strep_difference$difference <- log2(as.numeric(as.character(strep_difference$difference)) + 1)
strep_genus_diff <- merge(strep_genus, strep_difference, by.x='row.names', by.y='genus')
rownames(strep_genus_diff) <- strep_genus_diff$Row.names
strep_genus_diff$Row.names <- NULL
colnames(strep_genus_diff) <- c('relAbund','transcriptChange')
rm(strep_outliers, strep_difference, strep_genus)
cef_outliers <- rbind(cef_630_outliers, cef_mock_outliers)
cef_difference <- abs(((2^cef_outliers$cef_mock_metaT_reads) - 1) - ((2^cef_outliers$cef_630_metaT_reads) - 1))
cef_difference <- as.data.frame(cbind(cef_outliers$genus, cef_difference))
colnames(cef_difference) <- c('genus', 'difference')
cef_gene_count <- as.data.frame(table(cef_difference$genus))
colnames(cef_gene_count) <- c('genus', 'count')
cef_difference <- aggregate(. ~ genus, data=cef_difference, FUN=sum)
cef_difference <- merge(cef_difference, cef_gene_count, by='genus')
cef_difference$difference <- cef_difference$difference / cef_difference$count
cef_difference$count <- NULL
rm(cef_gene_count)
cef_difference$difference <- log2(as.numeric(as.character(cef_difference$difference)) + 1)
cef_genus_diff <- merge(cef_genus, cef_difference, by.x='row.names', by.y='genus')
rownames(cef_genus_diff) <- cef_genus_diff$Row.names
cef_genus_diff$Row.names <- NULL
colnames(cef_genus_diff) <- c('relAbund','transcriptChange')
rm(cef_outliers, cef_difference, cef_genus)
clinda_outliers <- rbind(clinda_630_outliers, clinda_mock_outliers)
clinda_difference <- abs(((2^clinda_outliers$clinda_mock_metaT_reads) - 1) - ((2^clinda_outliers$clinda_630_metaT_reads) - 1))
clinda_difference <- as.data.frame(cbind(clinda_outliers$genus, clinda_difference))
colnames(clinda_difference) <- c('genus', 'difference')
clinda_gene_count <- as.data.frame(table(clinda_difference$genus))
colnames(clinda_gene_count) <- c('genus', 'count')
clinda_difference <- aggregate(. ~ genus, data=clinda_difference, FUN=sum)
clinda_difference <- merge(clinda_difference, clinda_gene_count, by='genus')
clinda_difference$difference <- clinda_difference$difference / clinda_difference$count
clinda_difference$count <- NULL
rm(clinda_gene_count)
clinda_difference$difference <- log2(as.numeric(as.character(clinda_difference$difference)) + 1)
clinda_genus_diff <- merge(clinda_genus, clinda_difference, by.x='row.names', by.y='genus')
rownames(clinda_genus_diff) <- clinda_genus_diff$Row.names
clinda_genus_diff$Row.names <- NULL
colnames(clinda_genus_diff) <- c('relAbund','transcriptChange')
rm(clinda_outliers, clinda_difference, clinda_genus)

# Add color for genera
cef_genus_diff <- merge(cef_genus_diff, tax_colors, by.x='row.names', by.y='taxonomy', all.x=TRUE)
rownames(cef_genus_diff) <- cef_genus_diff$Row.names
cef_genus_diff$Row.names <- NULL
cef_genus_diff[is.na(cef_genus_diff)] <- 'white'
clinda_genus_diff <- merge(clinda_genus_diff, tax_colors, by.x='row.names', by.y='taxonomy', all.x=TRUE)
rownames(clinda_genus_diff) <- clinda_genus_diff$Row.names
clinda_genus_diff$Row.names <- NULL
clinda_genus_diff[is.na(clinda_genus_diff)] <- 'white'
strep_genus_diff <- merge(strep_genus_diff, tax_colors, by.x='row.names', by.y='taxonomy', all.x=TRUE)
rownames(strep_genus_diff) <- strep_genus_diff$Row.names
strep_genus_diff$Row.names <- NULL
strep_genus_diff[is.na(strep_genus_diff)] <- 'white'
rm(tax_colors)

# Convert relative abundance to categorical variable of ranges
cef_genus_diff_01 <- subset(cef_genus_diff, cef_genus_diff$relAbund < 0.1)
cef_genus_diff_01_1 <- subset(cef_genus_diff, cef_genus_diff$relAbund >= 0.1 & cef_genus_diff$relAbund < 1)
cef_genus_diff_1_10 <- subset(cef_genus_diff, cef_genus_diff$relAbund >= 1 & cef_genus_diff$relAbund <= 10)
cef_genus_diff_10_100 <- subset(cef_genus_diff, cef_genus_diff$relAbund > 10 & cef_genus_diff$relAbund <= 100)
rm(cef_genus_diff)
clinda_genus_diff_01 <- subset(clinda_genus_diff, clinda_genus_diff$relAbund < 0.1)
clinda_genus_diff_01_1 <- subset(clinda_genus_diff, clinda_genus_diff$relAbund >= 0.1 & clinda_genus_diff$relAbund < 1)
clinda_genus_diff_1_10 <- subset(clinda_genus_diff, clinda_genus_diff$relAbund >= 1 & clinda_genus_diff$relAbund <= 10)
clinda_genus_diff_10_100 <- subset(clinda_genus_diff, clinda_genus_diff$relAbund > 10 & clinda_genus_diff$relAbund <= 100)
rm(clinda_genus_diff)
strep_genus_diff_01 <- subset(strep_genus_diff, strep_genus_diff$relAbund < 0.1)
strep_genus_diff_01_1 <- subset(strep_genus_diff, strep_genus_diff$relAbund >= 0.1 & strep_genus_diff$relAbund < 1)
strep_genus_diff_1_10 <- subset(strep_genus_diff, strep_genus_diff$relAbund >= 1 & strep_genus_diff$relAbund <= 10)
strep_genus_diff_10_100 <- subset(strep_genus_diff, strep_genus_diff$relAbund > 10 & strep_genus_diff$relAbund <= 100)
rm(strep_genus_diff)

# Subset out other bacteria group
strep_630_outliers_other <- subset(strep_630_outliers, color == 'white')
strep_630_outliers <- subset(strep_630_outliers, color != 'white')
strep_mock_outliers_other <- subset(strep_mock_outliers, color == 'white')
strep_mock_outliers <- subset(strep_mock_outliers, color != 'white')
cef_630_outliers_other <- subset(cef_630_outliers, color == 'white')
cef_630_outliers <- subset(cef_630_outliers, color != 'white')
cef_mock_outliers_other <- subset(cef_mock_outliers, color == 'white')
cef_mock_outliers <- subset(cef_mock_outliers, color != 'white')
clinda_630_outliers_other <- subset(clinda_630_outliers, color == 'white')
clinda_630_outliers <- subset(clinda_630_outliers, color != 'white')
clinda_mock_outliers_other <- subset(clinda_mock_outliers, color == 'white')
clinda_mock_outliers <- subset(clinda_mock_outliers, color != 'white')

# Make sure Archeael points are visible
strep_630_archeae <- subset(strep_630_outliers, color == '#FF8000')
strep_mock_archeae <- subset(strep_mock_outliers, color == '#FF8000')
cef_630_archeae <- subset(cef_630_outliers, color == '#FF8000')
cef_mock_archeae <- subset(cef_mock_outliers, color == '#FF8000')
clinda_630_archeae <- subset(clinda_630_outliers, color == '#FF8000')
clinda_mock_archeae <- subset(clinda_mock_outliers, color == '#FF8000')

# Make sure Actinobacteria points are visible
strep_630_actino <- rbind(subset(strep_630_outliers, color == '#009900'), 
                          subset(strep_630_outliers, color == '#006600'), 
                          subset(strep_630_outliers, color == '#33FF33'))
strep_mock_actino <- rbind(subset(strep_mock_outliers, color == '#009900'), 
                           subset(strep_mock_outliers, color == '#006600'), 
                           subset(strep_mock_outliers, color == '#33FF33'))
cef_630_actino <- rbind(subset(cef_630_outliers, color == '#009900'), 
                        subset(cef_630_outliers, color == '#006600'), 
                        subset(cef_630_outliers, color == '#33FF33'))
cef_mock_actino <- rbind(subset(cef_mock_outliers, color == '#009900'), 
                         subset(cef_mock_outliers, color == '#006600'), 
                         subset(cef_mock_outliers, color == '#33FF33'))
clinda_630_actino <- rbind(subset(clinda_630_outliers, color == '#009900'), 
                           subset(clinda_630_outliers, color == '#006600'), 
                           subset(clinda_630_outliers, color == '#33FF33'))
clinda_mock_actino <- rbind(subset(clinda_mock_outliers, color == '#009900'), 
                            subset(clinda_mock_outliers, color == '#006600'), 
                            subset(clinda_mock_outliers, color == '#33FF33'))
clinda_mock_ecoli <- subset(clinda_mock_outliers, color == '#CCCC00')

#-------------------------------------------------------------------------------------------------------------------------#

# Perform some final calculations

# Spearman correlation of 16S and metatransctiptomic change
strep_16S_r <- as.character(round(cor.test(strep_difference$transcriptDiff, strep_difference$relAbund, method='spearman', exact=FALSE)$estimate, digits=5))
strep_16S_p <- as.character(round(cor.test(strep_difference$transcriptDiff, strep_difference$relAbund, method='spearman', exact=FALSE)$p.value, digits=5))
cef_16S_r <- as.character(round(cor.test(cef_difference$transcriptDiff, cef_difference$relAbund, method='spearman', exact=FALSE)$estimate, digits=5))
cef_16S_p <- as.character(round(cor.test(cef_difference$transcriptDiff, cef_difference$relAbund, method='spearman', exact=FALSE)$p.value, digits=5))
clinda_16S_r <- as.character(round(cor.test(clinda_difference$transcriptDiff, clinda_difference$relAbund, method='spearman', exact=FALSE)$estimate, digits=5))
clinda_16S_p <- as.character(round(cor.test(clinda_difference$transcriptDiff, clinda_difference$relAbund, method='spearman', exact=FALSE)$p.value, digits=5))

#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
tiff(filename=plot_file, width=10, height=10, units='in', 
     res=300, pointsize=12, compression='none')
layout(matrix(c(1,2,
                3,4), 
              nrow=2, ncol=2, byrow = TRUE))
par(mar=c(5, 5, 1, 1), mgp=c(3,0.7,0))

#-------------------#

# Streptomycin
plot(0, type='n', xlim=c(0,12), ylim=c(0,12), pch=20, xaxt='n', yaxt='n', xlab='', ylab='')
filledrectangle(wx=20, wy=2.8, col='gray80', mid=c(6,6), angle=45)
box(lwd=2)
points(x=strep_annotated$strep_mock_metaT_reads, y=strep_annotated$strep_630_metaT_reads, pch=20, cex=1.3, col='gray40')
segments(-2, -2, 14, 14, lty=2)
axis(1, at=seq(0,12,2), label=seq(0,12,2), cex=1.3)
axis(2, at=seq(0,12,2), label=seq(0,12,2), cex=1.3, las=1)
minor.ticks.axis(1, 10, mn=0, mx=12)
minor.ticks.axis(2, 10, mn=0, mx=12)
mtext(expression(paste('Normalized cDNA Abundance (',log[2],')')), side=1, padj=2.25, cex=0.9)
mtext('Mock-Infected', side=1, padj=3.7, font=2)
mtext(expression(paste('Normalized cDNA Abundance (',log[2],')')), side=2, padj=-1.75, cex=0.9)
mtext(expression(bolditalic('C. difficile')~bold('-Infected')), side=2, padj=-3.5, font=2)
legend('topleft', c('Streptomycin-pretreated', as.expression(bquote(paste(italic('R'),' = ',.(strep_corr))))), 
       bty='n', cex=1.3, text.col=c(strep_col,'black'))
mtext('A', side=2, line=2, las=2, adj=1, padj=-9, cex=1.7, font=2)

points(x=strep_630_outliers_other$strep_mock_metaT_reads, y=strep_630_outliers_other$strep_630_metaT_reads, cex=1.7, pch=21, col='gray10', lwd=1.5, bg=strep_630_outliers_other$color)
points(x=strep_mock_outliers_other$strep_mock_metaT_reads, y=strep_mock_outliers_other$strep_630_metaT_reads, cex=1.7, pch=21, col='gray10', lwd=1.5, bg=strep_mock_outliers_other$color)
points(x=strep_630_outliers$strep_mock_metaT_reads, y=strep_630_outliers$strep_630_metaT_reads, cex=1.7, pch=21, bg=strep_630_outliers$color, col='gray10')
points(x=strep_mock_outliers$strep_mock_metaT_reads, y=strep_mock_outliers$strep_630_metaT_reads, cex=1.7, pch=21, bg=strep_mock_outliers$color, col='gray10')
points(x=strep_630_actino$strep_mock_metaT_reads, y=strep_630_actino$strep_630_metaT_reads, cex=1.7, pch=21, bg=strep_630_actino$color, col='gray10')
points(x=strep_mock_actino$strep_mock_metaT_reads, y=strep_mock_actino$strep_630_metaT_reads, cex=1.7, pch=21, bg=strep_mock_actino$color, col='gray10')
points(x=strep_630_archeae$strep_mock_metaT_reads, y=strep_630_archeae$strep_630_metaT_reads, cex=1.7, pch=21, bg=strep_630_archeae$color, col='gray10')
points(x=strep_mock_archeae$strep_mock_metaT_reads, y=strep_mock_archeae$strep_630_metaT_reads, cex=1.7, pch=21, bg=strep_mock_archeae$color, col='gray10')

#-------------------#

# Cefoperazone
plot(0, type='n', xlim=c(0,12), ylim=c(0,12), pch=20, xaxt='n', yaxt='n', xlab='', ylab='')
filledrectangle(wx=20, wy=2.8, col='gray80', mid=c(6,6), angle=45)
box(lwd=2)
points(x=cef_annotated$cef_mock_metaT_reads, y=cef_annotated$cef_630_metaT_reads, pch=20, cex=1.3, col='gray40')
segments(-2, -2, 14, 14, lty=2)
axis(1, at=seq(0,12,2), label=seq(0,12,2), cex=1.3)
axis(2, at=seq(0,12,2), label=seq(0,12,2), cex=1.3, las=1)
minor.ticks.axis(1, 10, mn=0, mx=12)
minor.ticks.axis(2, 10, mn=0, mx=12)
mtext(expression(paste('Normalized cDNA Abundance (',log[2],')')), side=1, padj=2.25, cex=0.9)
mtext('Mock-Infected', side=1, padj=3.7, font=2)
mtext(expression(paste('Normalized cDNA Abundance (',log[2],')')), side=2, padj=-1.75, cex=0.9)
mtext(expression(bolditalic('C. difficile')~bold('-Infected')), side=2, padj=-3.5, font=2)
legend('topleft', c('Cefoperazone-pretreated', as.expression(bquote(paste(italic('R'),' = ',.(cef_corr))))), 
       bty='n', cex=1.3, text.col=c(cef_col,'black'))
mtext('B', side=2, line=2, las=2, adj=1, padj=-9, cex=1.7, font=2)

points(x=cef_630_outliers_other$cef_mock_metaT_reads, y=cef_630_outliers_other$cef_630_metaT_reads, cex=1.7, pch=21, col='gray10', lwd=1.5, bg=cef_630_outliers_other$color)
points(x=cef_mock_outliers_other$cef_mock_metaT_reads, y=cef_mock_outliers_other$cef_630_metaT_reads, cex=1.7, pch=21, col='gray10', lwd=1.5, bg=cef_mock_outliers_other$color)
points(x=cef_630_outliers$cef_mock_metaT_reads, y=cef_630_outliers$cef_630_metaT_reads, cex=1.7, pch=21, bg=cef_630_outliers$color, col='gray10')
points(x=cef_mock_outliers$cef_mock_metaT_reads, y=cef_mock_outliers$cef_630_metaT_reads, cex=1.7, pch=21, bg=cef_mock_outliers$color, col='gray10')
points(x=cef_630_actino$cef_mock_metaT_reads, y=cef_630_actino$cef_630_metaT_reads, cex=1.7, pch=21, bg=cef_630_actino$color, col='gray10')
points(x=cef_mock_actino$cef_mock_metaT_reads, y=cef_mock_actino$cef_630_metaT_reads, cex=1.7, pch=21, bg=cef_mock_actino$color, col='gray10')
points(x=cef_630_archeae$cef_mock_metaT_reads, y=cef_630_archeae$cef_630_metaT_reads, cex=1.7, pch=21, bg=cef_630_archeae$color, col='gray10')
points(x=cef_mock_archeae$cef_mock_metaT_reads, y=cef_mock_archeae$cef_630_metaT_reads, cex=1.7, pch=21, bg=cef_mock_archeae$color, col='gray10')

#-------------------#

# Clindamycin
plot(0, type='n', xlim=c(0,12), ylim=c(0,12), pch=20, xaxt='n', yaxt='n', xlab='', ylab='')
filledrectangle(wx=20, wy=2.8, col='gray80', mid=c(6,6), angle=45)
box(lwd=2)
points(x=clinda_annotated$clinda_mock_metaT_reads, y=clinda_annotated$clinda_630_metaT_reads, pch=20, cex=1.3, col='gray40')
segments(-2, -2, 14, 14, lty=2)
axis(1, at=seq(0,12,2), label=seq(0,12,2), cex=1.3)
axis(2, at=seq(0,12,2), label=seq(0,12,2), cex=1.3, las=1)
minor.ticks.axis(1, 10, mn=0, mx=12)
minor.ticks.axis(2, 10, mn=0, mx=12)
mtext(expression(paste('Normalized cDNA Abundance (',log[2],')')), side=1, padj=2.25, cex=0.9)
mtext('Mock-Infected', side=1, padj=3.7, font=2)
mtext(expression(paste('Normalized cDNA Abundance (',log[2],')')), side=2, padj=-1.75, cex=0.9)
mtext(expression(bolditalic('C. difficile')~bold('-Infected')), side=2, padj=-3.5, font=2)
legend('topleft', c('Clindamycin-pretreated', as.expression(bquote(paste(italic('R'),' = ',.(clinda_corr))))), 
       bty='n', cex=1.3, text.col=c(clinda_col,'black'))
mtext('C', side=2, line=2, las=2, adj=1, padj=-9, cex=1.7, font=2)

points(x=clinda_630_outliers_other$clinda_mock_metaT_reads, y=clinda_630_outliers_other$clinda_630_metaT_reads, cex=1.7, pch=21, col='gray10', lwd=1.5, bg=clinda_630_outliers_other$color)
points(x=clinda_mock_outliers_other$clinda_mock_metaT_reads, y=clinda_mock_outliers_other$clinda_630_metaT_reads, cex=1.7, pch=21, col='gray10', lwd=1.5, bg=clinda_mock_outliers_other$color)
points(x=clinda_630_outliers$clinda_mock_metaT_reads, y=clinda_630_outliers$clinda_630_metaT_reads, cex=1.7, pch=21, bg=clinda_630_outliers$color, col='gray10')
points(x=clinda_mock_outliers$clinda_mock_metaT_reads, y=clinda_mock_outliers$clinda_630_metaT_reads, cex=1.7, pch=21, bg=clinda_mock_outliers$color, col='gray10')
points(x=clinda_630_actino$clinda_mock_metaT_reads, y=clinda_630_actino$clinda_630_metaT_reads, cex=1.7, pch=21, bg=clinda_630_actino$color, col='gray10')
points(x=clinda_mock_actino$clinda_mock_metaT_reads, y=clinda_mock_actino$clinda_630_metaT_reads, cex=1.7, pch=21, bg=clinda_mock_actino$color, col='gray10')
points(x=clinda_630_archeae$clinda_mock_metaT_reads, y=clinda_630_archeae$clinda_630_metaT_reads, cex=1.7, pch=21, bg=clinda_630_archeae$color, col='gray10')
points(x=clinda_mock_archeae$clinda_mock_metaT_reads, y=clinda_mock_archeae$clinda_630_metaT_reads, cex=1.7, pch=21, bg=clinda_mock_archeae$color, col='gray10')

#-------------------#

# Taxonomic group legend
par(mar=c(0,1,3,0))
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-5,5), ylim=c(-10,10))
rect(xleft=-4.5, ybottom=10, xright=3.5, ytop=-5.5, border='black')

# Left side
text(x=c(-3,-3,1,1,1,1,0.5), y=c(9,2.4,9,5.5,2.4,0,-2.4), labels=c('Bacteroidetes', 'Firmicutes','Actinobacteria','Proteobacteria','Verrucomicrobia','Other Bacteria', 'Archeae'), cex=1.2) # Phyla
text(x=-2.5, y=c(8.3,7.6,6.9,6.2,5.5,4.8), labels=c('Allistipes','Bacteroides','Odoribacter','Parabacteroides','Prevotella','Porphymonas'), cex=0.9, font=3) # Bacteroietes
points(x=rep(-1.1,6), y=c(8.3,7.6,6.9,6.2,5.5,4.8), pch=22, cex=2.1, col='black', bg=c('#000099','#0000CC','#0000FF','#3333FF','#6666FF','#9999FF')) # blues
text(x=-2.5, y=c(1.7,1,0.3,-0.4,-1.1,-1.8,-2.5,-3.2,-3.9), labels=c('Clostridium','Enterococcus','Eubacterium','Lactobacillus','Lactococcus','Roseburia','Ruminococcus','Staphylococcus','Streptococcus'), cex=0.9, font=3) # Firmicutes
points(x=rep(-1.1,9), y=c(1.7,1,0.3,-0.4,-1.1,-1.8,-2.5,-3.2,-3.9), pch=22, cex=2.1, col='black', bg=c('#330000','#660000','#990000','#CC0000','#FF0000','#FF3333','#FF6666','#FF9999','#FFCCCC')) # reds

# Right side
text(x=1.45, y=c(8.3,7.6,6.9), labels=c('Bifidobacterium','Corynebacterium','Olsenella'), cex=0.9, font=3) # Actinobacteria
points(x=rep(2.75,3), y=c(8.3,7.6,6.9), pch=22, cex=2.2, col='black', bg=c('#006600','#009900','#33FF33')) # greens
text(x=1.5, y=4.8, labels='Escherichia', cex=0.9, font=3) # Proteobacteria
points(x=2.75, y=4.8, pch=22, cex=2.2, col='black', bg='#CCCC00') # yellow
text(x=1.5, y=1.7, labels='Akkermansia', cex=0.9, font=3) # Verrucomicrobia
points(x=2.75, y=1.7, pch=22, cex=2.2, col='black', bg='#990099') # purple
text(x=1.5, y=-0.85, labels='<=0.01% Each', cex=0.9) # Other bacteria
points(x=2.75, y=-0.85, pch=22, cex=2.2, col='black', bg='white') # Other Bacteria - white
text(x=1.3, y=-3.1, labels='Methanobrevibacter', cex=0.9, font=3) # Archeae
points(x=2.75, y=-3.1, pch=22, cex=2.2, col='black', bg='#FF8000') # orange

dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)}
setwd(starting_dir)
rm(list=ls())
gc()
