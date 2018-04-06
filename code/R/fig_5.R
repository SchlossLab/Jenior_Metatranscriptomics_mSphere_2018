
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

# Output plot
plot_file <- 'results/figures/figure_5.pdf'

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
kegg_tax$domain <- NULL
kegg_tax$group <- NULL
kegg_tax$species <- NULL

# 16S data
genus_shared <- read.delim(genus_shared, sep='\t', header=TRUE, row.names=2)
genus_shared$label <- NULL
genus_shared$numOtus <- NULL
genus_tax <- read.delim(genus_tax, sep='\t', header=TRUE)
rownames(genus_tax) <- genus_tax$OTU
genus_tax$OTU <- NULL

#-------------------------------------------------------------------------------------------------------------------------#

# Format 16S
# Collate 16S abundances within genera in each sample
genus_shared <- t(genus_shared)
genus_tax_shared <- clean_merge(genus_tax, genus_shared)
genus_tax_agg_shared <- aggregate(. ~ Genus, data=genus_tax_shared, FUN=sum)
rownames(genus_tax_agg_shared) <- genus_tax_agg_shared$Genus
genus_tax_agg_shared$Genus <- NULL
genus_tax_agg_shared <- as.data.frame(t(genus_tax_agg_shared))
rm(genus_tax, genus_shared, genus_tax_shared)
# Combine with metadata and collate within antibiotic treatment
genus_tax_agg_shared <- clean_merge(metadata, genus_tax_agg_shared)
genus_tax_agg_shared <- subset(genus_tax_agg_shared, abx %in% c('cefoperazone','clindamycin','streptomycin'))
genus_tax_agg_shared <- subset(genus_tax_agg_shared, infection == 'mock')
genus_tax_agg_shared$cage <- NULL
genus_tax_agg_shared$mouse <- NULL
genus_tax_agg_shared$gender <- NULL
genus_tax_agg_shared$type <- NULL
genus_tax_agg_shared$susceptibility <- NULL
genus_tax_agg_shared$infection <- NULL
abx_genus_tax_agg_shared <- aggregate(. ~ abx, data=genus_tax_agg_shared, FUN=median)
rownames(abx_genus_tax_agg_shared) <- abx_genus_tax_agg_shared$abx
abx_genus_tax_agg_shared$abx <- NULL
rm(genus_tax_agg_shared, metadata)
# Filter to minority genera in each group (<1%)
abx_genus_tax_agg_shared <- abx_genus_tax_agg_shared[, which(colSums(abx_genus_tax_agg_shared) > 0)]
abx_genus_tax_agg_shared <- (abx_genus_tax_agg_shared / rowSums(abx_genus_tax_agg_shared)) * 100 # Convert to relative abundance
abx_genus_tax_agg_shared[abx_genus_tax_agg_shared > 1] <- 0 # 1 percent cutoff
minority_genera <- abx_genus_tax_agg_shared[, which(colSums(abx_genus_tax_agg_shared) != 0)]
rm(abx_genus_tax_agg_shared)
# Subset to each antibiotic treatment
minority_genera <- as.data.frame(t(minority_genera))
strep_minority_genera <- rownames(minority_genera[which(minority_genera$streptomycin > 0),])
cef_minority_genera <- rownames(minority_genera[which(minority_genera$cefoperazone > 0),])
clinda_minority_genera <- rownames(minority_genera[which(minority_genera$clindamycin > 0),])
rm(minority_genera)

# Format transcription
# Remove those genes without an organism annotation
cef_annotated <- cef_normalized_reads[!is.na(cef_normalized_reads$organism),]
clinda_annotated <- clinda_normalized_reads[!is.na(clinda_normalized_reads$organism),]
strep_annotated <- strep_normalized_reads[!is.na(strep_normalized_reads$organism),]
rm(cef_normalized_reads, clinda_normalized_reads, strep_normalized_reads)
# Remove any C. difficile genes with transcription only in infected
cdiff_omit <- c('cdf','pdc','cdc','cdl','pdf')
cef_annotated <- subset(cef_annotated, !(cef_annotated$organism %in% cdiff_omit))
clinda_annotated <- subset(clinda_annotated, !(clinda_annotated$organism %in% cdiff_omit))
strep_annotated <- subset(strep_annotated, !(strep_annotated$organism %in% cdiff_omit))
rm(cdiff_omit)
# Remove non-microbial genes
mamm_omit <- c('fab','cfa','ggo','hgl','hsa','mcc','mdo','pon','aml',
               'ptr','rno','shr','ssc','aml','bta','cge','ecb',
               'pps','fca','mmu','oaa','gga','ola','acs','aga')
strep_annotated <- subset(strep_annotated, !(strep_annotated$organism %in% mamm_omit))
cef_annotated <- subset(cef_annotated, !(cef_annotated$organism %in% mamm_omit))
clinda_annotated <- subset(clinda_annotated, !(clinda_annotated$organism %in% mamm_omit))
rm(mamm_omit)
# Merge with genus level taxonomic classifications
strep_annotated$gene <- rownames(strep_annotated)
strep_annotated <- merge(kegg_tax, strep_annotated, by.x='org_code', by.y='organism', all.y=TRUE)
cef_annotated$gene <- rownames(cef_annotated)
cef_annotated <- merge(kegg_tax, cef_annotated, by.x='org_code', by.y='organism', all.y=TRUE)
clinda_annotated$gene <- rownames(clinda_annotated)
clinda_annotated <- merge(kegg_tax, clinda_annotated, by.x='org_code', by.y='organism', all.y=TRUE)
rm(kegg_tax)

#-------------------------------------------------------------------------------------------------------------------------#

# Tanscriptomic analysis of minority taxa
# Subset by minority genera from 16S analysis
strep_minority <- subset(strep_annotated, genus %in% strep_minority_genera)
cef_minority <- subset(cef_annotated, genus %in% cef_minority_genera)
clinda_minority <- subset(clinda_annotated, genus %in% clinda_minority_genera)
rm(strep_annotated, cef_annotated, clinda_annotated,
   strep_minority_genera, cef_minority_genera, clinda_minority_genera)
# Remove unused columns
strep_minority$org_code <- NULL
strep_minority$genus <- NULL
strep_minority$description <- NULL
strep_minority$ko <- NULL
strep_minority$gene <- NULL
cef_minority$org_code <- NULL
cef_minority$genus <- NULL
cef_minority$description <- NULL
cef_minority$ko <- NULL
cef_minority$gene <- NULL
clinda_minority$org_code <- NULL
clinda_minority$genus <- NULL
clinda_minority$description <- NULL
clinda_minority$ko <- NULL
clinda_minority$gene <- NULL
# Remove genes with no pathway annotation
strep_minority <- strep_minority[complete.cases(strep_minority), ]
cef_minority <- cef_minority[complete.cases(cef_minority), ]
clinda_minority <- clinda_minority[complete.cases(clinda_minority), ]
# Reformat pathway names
strep_minority$pathways <- vapply(strsplit(as.character(strep_minority$pathways),':'), `[`, 1, FUN.VALUE=character(1))
strep_minority$pathways <- gsub('_', ' ', strep_minority$pathways)
cef_minority$pathways <- vapply(strsplit(as.character(cef_minority$pathways),':'), `[`, 1, FUN.VALUE=character(1))
cef_minority$pathways <- gsub('_', ' ', cef_minority$pathways)
clinda_minority$pathways <- vapply(strsplit(as.character(clinda_minority$pathways),':'), `[`, 1, FUN.VALUE=character(1))
clinda_minority$pathways <- gsub('_', ' ', clinda_minority$pathways)
# Aggregate by pathway
strep_minority <- aggregate(. ~ pathways, data=strep_minority, FUN=sum)
rownames(strep_minority) <- strep_minority$pathways
strep_minority$pathways <- NULL
cef_minority <- aggregate(. ~ pathways, data=cef_minority, FUN=sum)
rownames(cef_minority) <- cef_minority$pathways
cef_minority$pathways <- NULL
clinda_minority <- aggregate(. ~ pathways, data=clinda_minority, FUN=sum)
rownames(clinda_minority) <- clinda_minority$pathways
clinda_minority$pathways <- NULL
# Remove groups with no difference and rank differences
strep_minority$abs_diff <- abs(strep_minority[,1] - strep_minority[,2])
strep_minority <- subset(strep_minority, abs_diff != 0)
strep_minority <- strep_minority[order(-strep_minority$abs_diff),]
cef_minority$abs_diff <- abs(cef_minority[,1] - cef_minority[,2])
cef_minority <- subset(cef_minority, abs_diff != 0)
cef_minority <- cef_minority[order(-cef_minority$abs_diff),]
clinda_minority$abs_diff <- abs(clinda_minority[,1] - clinda_minority[,2])
clinda_minority <- subset(clinda_minority, abs_diff != 0)
clinda_minority <- clinda_minority[order(-clinda_minority$abs_diff),]
# Remove useless pathways
strep_minority <- subset(strep_minority, !rownames(strep_minority) %in% c('Ribosome','Metabolic pathways'))
cef_minority <- subset(cef_minority, !rownames(cef_minority) %in% c('Ribosome','Metabolic pathways'))
clinda_minority <- subset(clinda_minority, !rownames(clinda_minority) %in% c('Ribosome','Metabolic pathways'))

# Calculate order of magnitude transcriptional shift
strep_cef_diff <- round(log2(abs(sum(strep_minority$abs_diff) - sum(cef_minority$abs_diff))), 2)
strep_clinda_diff <- round(log2(abs(sum(strep_minority$abs_diff) - sum(clinda_minority$abs_diff))), 2)
cef_clinda_diff <- round(log2(abs(sum(cef_minority$abs_diff) - sum(clinda_minority$abs_diff))), 2)
clearedVcolonized_diff <- round(log2(abs(mean(c(sum(strep_minority$abs_diff), sum(cef_minority$abs_diff))) - sum(clinda_minority$abs_diff))), 2)
# Remove absolute difference column
strep_minority$abs_diff <- NULL
cef_minority$abs_diff <- NULL
clinda_minority$abs_diff <- NULL
# Log2 transform
strep_minority <- log2(strep_minority + 1)
cef_minority <- log2(cef_minority + 1)
clinda_minority <- log2(clinda_minority + 1)
# Subset to top hits at level of minimum between groups
path_level <- as.numeric(min(nrow(strep_minority), nrow(cef_minority), nrow(clinda_minority)))
strep_minority <- strep_minority[1:path_level,]
cef_minority <- cef_minority[1:path_level,]
clinda_minority <- clinda_minority[1:path_level,]
rm(path_level)
# Calculate overlap of most altered pathways
strep_cef_overlap_top <- intersect(rownames(strep_minority), rownames(cef_minority))
strep_clinda_overlap_top <- intersect(rownames(strep_minority), rownames(clinda_minority))
cef_clinda_overlap_top <- intersect(rownames(cef_minority), rownames(clinda_minority))
clearedVcolonized_distinct_top <- setdiff(intersect(rownames(strep_minority), rownames(cef_minority)), rownames(clinda_minority))
clearedVcolonized_overlap_top <- intersect(intersect(rownames(strep_minority), rownames(cef_minority)), rownames(clinda_minority))

#-------------------------------------------------------------------------------------------------------------------------#

# Get data ready for plotting
#strep_minority$pathway <- rownames(strep_minority)
#colnames(strep_minority) <- c('infected','mock','pathway',)
#cef_minority$pathway <- rownames(cef_minority)
#colnames(cef_minority) <- c('infected','mock','pathway')
#clinda_minority$pathway <- rownames(clinda_minority)
#colnames(clinda_minority) <- c('infected','mock','pathway')
#minority_pathways <- as.data.frame(rbind(strep_minority, cef_minority, clinda_minority))
#minority_pathways$colors <- c(rep(strep_col,nrow(strep_minority)), rep(cef_col,nrow(cef_minority)), rep(clinda_col,nrow(clinda_minority)))
#rm(strep_minority, cef_minority, clinda_minority)

#-------------------------------------------------------------------------------------------------------------------------#

# Generate figure
pdf(filename=plot_file, width=15, height=5)
layout(matrix(c(1,
                2,
                3), 
              nrow=3, ncol=1, byrow = TRUE))
par(mar=c(4, 4, 1, 1), mgp=c(3,0.7,0))

#-------------------#


# Streptomycin




# Cefoperazone


abline(h=, lty=1) # separates cleared and colonized


# Clindamycin





dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
#for (dep in deps){
#  pkg <- paste('package:', dep, sep='')
#  detach(pkg, character.only = TRUE)}
#setwd(starting_dir)
#rm(list=ls())
#gc()
