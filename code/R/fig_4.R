
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
plot_file <- 'results/figures/figure_4.pdf'

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

# Collate 16S abundances at genus-level 
genus_shared <- t(genus_shared)
genus_tax_shared <- merge(genus_tax, genus_shared, by='row.names')
genus_tax_shared$Row.names <- NULL
genus_tax_agg_shared <- aggregate(. ~ Genus, data=genus_tax_shared, FUN=sum)
rownames(genus_tax_agg_shared) <- genus_tax_agg_shared$Genus
genus_tax_agg_shared$Genus <- NULL
genus_tax_agg_shared <- t(genus_tax_agg_shared)
rm(genus_tax, genus_shared, genus_tax_shared)

# Combine with metadata, subset to correct groups
genus_tax_agg_metadata_shared <- merge(metadata, genus_tax_agg_shared, by='row.names')
genus_tax_agg_metadata_shared$Row.names <- NULL
genus_tax_agg_metadata_shared$cage <- NULL
genus_tax_agg_metadata_shared$mouse <- NULL
genus_tax_agg_metadata_shared$gender <- NULL
genus_tax_agg_metadata_shared$type <- NULL
genus_tax_agg_metadata_shared$susceptibility <- NULL
genus_tax_agg_metadata_shared$infection <- NULL
genus_tax_agg_metadata_shared <- subset(genus_tax_agg_metadata_shared, genus_tax_agg_metadata_shared$abx %in% c('cefoperazone','clindamycin','streptomycin'))

# Compare richness of lowest abundance genera
strep_shared <- subset(genus_tax_agg_metadata_shared, abx == 'streptomycin')
strep_shared$abx <- NULL
strep_shared <- (strep_shared / rowSums(strep_shared)) * 100
strep_otu_medians <- apply(strep_shared, 2, median)
strep_shared <- strep_shared[,which(strep_otu_medians <= 1)]
strep_otu_medians <- apply(strep_shared, 2, median)
strep_shared <- strep_shared[,which(strep_otu_medians > 0)]
strep_minority_genera <- as.numeric(ncol(strep_shared))
strep_shared <- as.data.frame(t(strep_shared))
cef_shared <- subset(genus_tax_agg_metadata_shared, abx == 'cefoperazone')
cef_shared$abx <- NULL
cef_shared <- (cef_shared / rowSums(cef_shared)) * 100
cef_otu_medians <- apply(cef_shared, 2, median)
cef_shared <- cef_shared[,which(cef_otu_medians <= 1)]
cef_otu_medians <- apply(cef_shared, 2, median)
cef_shared <- cef_shared[,which(cef_otu_medians > 0)]
cef_minority_genera <- as.numeric(ncol(cef_shared))
cef_shared <- as.data.frame(t(cef_shared))
clinda_shared <- subset(genus_tax_agg_metadata_shared, abx == 'clindamycin')
clinda_shared$abx <- NULL
clinda_shared <- (clinda_shared / rowSums(clinda_shared)) * 100
clinda_otu_medians <- apply(clinda_shared, 2, median)
clinda_shared <- clinda_shared[,which(clinda_otu_medians <= 1)]
clinda_otu_medians <- apply(clinda_shared, 2, median)
clinda_shared <- clinda_shared[,which(clinda_otu_medians > 0)]
clinda_minority_genera <- as.numeric(ncol(clinda_shared))
clinda_shared <- as.data.frame(t(clinda_shared))
rm(strep_otu_medians, cef_otu_medians, clinda_otu_medians)

# Combine shared tables and and compare jaccard beta diversity
strep_cef_shared <- merge(strep_shared, cef_shared, by='row.names', all=TRUE)
rownames(strep_cef_shared) <- strep_cef_shared$Row.names
strep_cef_shared$Row.names <- NULL
strep_cef_shared[is.na(strep_cef_shared)] <- 0
strep_cef_shared <- as.data.frame(t(strep_cef_shared))
strep_cef_dist <- vegdist(strep_cef_shared, method='jaccard')
strep_cef_shared$abx <- c(rep('strep',ncol(strep_shared)), rep('cef', ncol(cef_shared)))
strep_cef_jaccard_pval <- as.character(round(adonis(strep_cef_dist ~ strep_cef_shared$abx, strep_cef_shared[,1:ncol(strep_cef_shared)-1], perm=999)$aov.tab[[6]][1], 3))
rm(strep_cef_shared, strep_cef_dist)
strep_clinda_shared <- merge(strep_shared, clinda_shared, by='row.names', all=TRUE)
rownames(strep_clinda_shared) <- strep_clinda_shared$Row.names
strep_clinda_shared$Row.names <- NULL
strep_clinda_shared[is.na(strep_clinda_shared)] <- 0
strep_clinda_shared <- as.data.frame(t(strep_clinda_shared))
strep_clinda_dist <- vegdist(strep_clinda_shared, method='jaccard')
strep_clinda_shared$abx <- c(rep('strep',ncol(strep_shared)), rep('clinda', ncol(clinda_shared)))
strep_clinda_jaccard_pval <- as.character(round(adonis(strep_clinda_dist ~ strep_clinda_shared$abx, strep_clinda_shared[,1:ncol(strep_clinda_shared)-1], perm=999)$aov.tab[[6]][1], 3))
rm(strep_clinda_shared, strep_clinda_dist)
cef_clinda_shared <- merge(cef_shared, clinda_shared, by='row.names', all=TRUE)
rownames(cef_clinda_shared) <- cef_clinda_shared$Row.names
cef_clinda_shared$Row.names <- NULL
cef_clinda_shared[is.na(cef_clinda_shared)] <- 0
cef_clinda_shared <- as.data.frame(t(cef_clinda_shared))
cef_clinda_dist <- vegdist(cef_clinda_shared, method='jaccard')
cef_clinda_shared$abx <- c(rep('cef',ncol(cef_shared)), rep('clinda', ncol(clinda_shared)))
cef_clinda_jaccard_pval <- as.character(round(adonis(cef_clinda_dist ~ cef_clinda_shared$abx, cef_clinda_shared[,1:ncol(cef_clinda_shared)-1], perm=999)$aov.tab[[6]][1], 3))
rm(cef_clinda_shared, cef_clinda_dist)

# Collate within treatment
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
strep_extra <- subset(strep_difference, genus %in% c("Cellulosilyticum",'Butyrivibrio','Dorea',"Lactococcus","Alistipes","Odoribacter","Staphylococcus","Streptococcus","Roseburia","Ruminococcus","Porphyromonas","Faecalibacterium"))
strep_extra$streptomycin <- rep(0.001, nrow(strep_extra))
strep_extra <- strep_extra[,c(1,3,2)]
colnames(strep_extra) <- c('Row.names','streptomycin','difference')
strep_genus_diff <- as.data.frame(rbind(strep_genus_diff, strep_extra))
strep_genus_diff <- strep_genus_diff[match(unique(strep_genus_diff$Row.names), strep_genus_diff$Row.names),]
rownames(strep_genus_diff) <- strep_genus_diff$Row.names
strep_genus_diff$Row.names <- NULL
colnames(strep_genus_diff) <- c('relAbund','transcriptChange')
rm(strep_outliers, strep_difference, strep_genus, strep_extra)
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
cef_extra <- subset(cef_difference, genus %in% c('Escherichia', "Porphyromonas","Ruminococcus","Clostridium","Lactococcus","Odoribacter","Prevotella"))
cef_extra$ceftomycin <- rep(0.001, nrow(cef_extra))
cef_extra <- cef_extra[,c(1,3,2)]
colnames(cef_extra) <- c('Row.names','ceftomycin','difference')
cef_genus_diff <- as.data.frame(rbind(cef_genus_diff, cef_extra))
cef_genus_diff <- cef_genus_diff[match(unique(cef_genus_diff$Row.names), cef_genus_diff$Row.names),]
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
clinda_extra <- subset(clinda_difference, genus %in% c("Cellulosilyticum","Unclassified_Clostridiales","Clostridium","Lactococcus",'Butyrivibrio','Dorea','Roseburia','Bacteroides',"Prevotella","Ruminococcus","Porphyromonas","Faecalibacterium"))
clinda_extra$clindatomycin <- rep(0.001, nrow(clinda_extra))
clinda_extra <- clinda_extra[,c(1,3,2)]
colnames(clinda_extra) <- c('Row.names','clindatomycin','difference')
clinda_genus_diff <- as.data.frame(rbind(clinda_genus_diff, clinda_extra))
clinda_genus_diff <- clinda_genus_diff[match(unique(clinda_genus_diff$Row.names), clinda_genus_diff$Row.names),]
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
cef_genus_diff_01 <- subset(cef_genus_diff, cef_genus_diff$relAbund <= 0.1)
cef_genus_diff_01_1 <- subset(cef_genus_diff, cef_genus_diff$relAbund > 0.1 & cef_genus_diff$relAbund < 1)
cef_genus_diff_1_10 <- subset(cef_genus_diff, cef_genus_diff$relAbund >= 1 & cef_genus_diff$relAbund <= 10)
cef_genus_diff_10_100 <- subset(cef_genus_diff, cef_genus_diff$relAbund > 10 & cef_genus_diff$relAbund <= 100)
rm(cef_genus_diff)
clinda_genus_diff_01 <- subset(clinda_genus_diff, clinda_genus_diff$relAbund <= 0.1)
clinda_genus_diff_01_1 <- subset(clinda_genus_diff, clinda_genus_diff$relAbund > 0.1 & clinda_genus_diff$relAbund < 1)
clinda_genus_diff_1_10 <- subset(clinda_genus_diff, clinda_genus_diff$relAbund >= 1 & clinda_genus_diff$relAbund <= 10)
clinda_genus_diff_10_100 <- subset(clinda_genus_diff, clinda_genus_diff$relAbund > 10 & clinda_genus_diff$relAbund <= 100)
rm(clinda_genus_diff)
strep_genus_diff_01 <- subset(strep_genus_diff, strep_genus_diff$relAbund <= 0.1)
strep_genus_diff_01_1 <- subset(strep_genus_diff, strep_genus_diff$relAbund > 0.1 & strep_genus_diff$relAbund < 1)
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

# Plot the figure
pdf(file=plot_file, width=13, height=7)

# Transcript abundance by species abundance
par(mar=c(5, 5, 2, 2), yaxs='i')
plot(1, type='n', ylim=c(0,15), xlim=c(0,8),
     ylab='', xlab='', xaxt='n', yaxt='n', las=1)
axis(1, at=c(1,3,5,7), label=c('<=0.1%','>0.1% and <=1%','>1% and <=10%','>10% and <=100%'), 
     tick=FALSE, cex.axis=1.4, font=2)
axis(2, at=c(0,5,10,15), cex.axis=1.5, las=1)
mtext('Genus-level Relative Abundance (16S)', side=1, padj=3, cex=1.4)
mtext(expression(paste('Absolute Difference in Metatranscriptome (',log[2],')')), side=2, padj=-2, cex=1.4)
minor.ticks.axis(2, 20, mn=0, mx=25)
abline(v=c(2,4,6), lwd=2)
abline(v=c(0.5,1,1.5), lwd=2, lty=5, col=c(strep_col,cef_col,clinda_col))
abline(v=c(2.5,3,3.5), lwd=2, lty=5, col=c(strep_col,cef_col,clinda_col))
abline(v=c(4.5,5,5.5), lwd=2, lty=5, col=c(strep_col,cef_col,clinda_col))
abline(v=c(6.5,7,7.5), lwd=2, lty=5, col=c(strep_col,cef_col,clinda_col))
box(lwd=2)

# <0.1%
strep_01_n <- as.numeric(nrow(strep_genus_diff_01))
text(x=0.35, y=0.5, strep_01_n, cex=1.4)
other <- subset(strep_genus_diff_01, strep_genus_diff_01$color == 'white')
strep_genus_diff_01 <- subset(strep_genus_diff_01, strep_genus_diff_01$color != 'white')
stripchart(at=0.5, other$transcriptChange, pch=21, bg=other$color, 
           method='jitter', jitter=0.09, cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)
stripchart(at=0.5, rev(strep_genus_diff_01$transcriptChange), pch=21, bg=rev(strep_genus_diff_01$color), 
           method='jitter', jitter=0.09, cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)
temp <- subset(strep_genus_diff_01, rownames(strep_genus_diff_01) %in% c('Lactococcus','Roseburia','Prevotella','Akkermansia','Bifidobacterium') )
stripchart(at=0.5, temp$transcriptChange, pch=21, bg=temp$color, 
           method='jitter', jitter=0.09, cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)

cef_01_n <- as.numeric(nrow(cef_genus_diff_01))
text(x=0.85, y=0.5, cef_01_n, cex=1.4)
other <- subset(cef_genus_diff_01, cef_genus_diff_01$color == 'white')
cef_genus_diff_01 <- subset(cef_genus_diff_01, cef_genus_diff_01$color != 'white')
stripchart(at=1, other$transcriptChange, pch=21, bg=other$color, 
           method='jitter', jitter=0.09, cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)
stripchart(at=1, cef_genus_diff_01$transcriptChange, pch=21, bg=cef_genus_diff_01$color, 
           method='jitter', jitter=0.09, cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)
temp <- subset(cef_genus_diff_01, rownames(cef_genus_diff_01) %in% c('Allistipes','Odoribacter','Parabacteroides') )
stripchart(at=1, temp$transcriptChange, pch=21, bg=temp$color, 
           method='jitter', jitter=0.09, cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)

clinda_01_n <- as.numeric(nrow(clinda_genus_diff_01))
text(x=1.35, y=0.5, clinda_01_n, cex=1.4)
other <- subset(clinda_genus_diff_01, clinda_genus_diff_01$color == 'white')
clinda_genus_diff_01 <- subset(clinda_genus_diff_01, clinda_genus_diff_01$color != 'white')
stripchart(at=1.5, other$transcriptChange, pch=21, bg=other$color, 
           method='jitter', jitter=0.09, cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)
stripchart(at=1.5, clinda_genus_diff_01$transcriptChange, pch=21, bg=clinda_genus_diff_01$color, 
           method='jitter', jitter=0.09, cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)
temp <- subset(clinda_genus_diff_01, rownames(clinda_genus_diff_01) %in% c('Bifidobacterium','Parabacteroides') )
stripchart(at=1.5, temp$transcriptChange, pch=21, bg=temp$color, 
           method='jitter', jitter=0.09, cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)

# >0.1% and <1%
stripchart(at=2.5, strep_genus_diff_01_1$transcriptChange, pch=21, bg=strep_genus_diff_01_1$color, 
           cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)
stripchart(at=3, cef_genus_diff_01_1$transcriptChange, pch=21, bg=cef_genus_diff_01_1$color, 
           cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)
stripchart(at=3.5, clinda_genus_diff_01_1$transcriptChange, pch=21, bg=clinda_genus_diff_01_1$color, 
           cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)

# >1% and <10%
stripchart(at=4.5, strep_genus_diff_1_10$transcriptChange, pch=21, bg=strep_genus_diff_1_10$color, 
           cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)
stripchart(at=5, cef_genus_diff_1_10$transcriptChange, pch=21, bg=cef_genus_diff_1_10$color, 
           cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)
stripchart(at=5.5, clinda_genus_diff_1_10$transcriptChange, pch=21, bg=clinda_genus_diff_1_10$color, 
           cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)

# >10% and <100%
stripchart(at=6.5, strep_genus_diff_10_100$transcriptChange, pch=21, bg=strep_genus_diff_10_100$color, 
           cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)
stripchart(at=7, cef_genus_diff_10_100$transcriptChange, pch=21, bg=cef_genus_diff_10_100$color, 
           cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)
stripchart(at=7.5, clinda_genus_diff_10_100$transcriptChange, pch=21, bg=clinda_genus_diff_10_100$color, 
           cex=2.5, lwd=1.5, vertical=TRUE, add=TRUE)

par(xpd=TRUE)
legend(x=0.25, y=16.5, legend=c('Streptomycin-pretreated','Cefoperazone-pretreated','Clindamycin-pretreated'), 
       bty='n', col=c(strep_col,cef_col,clinda_col), pch=c(1,16), cex=1.5, pt.cex=0, lty=5, lwd=2.5, ncol=3)

dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)}
setwd(starting_dir)
rm(list=ls())
gc()
