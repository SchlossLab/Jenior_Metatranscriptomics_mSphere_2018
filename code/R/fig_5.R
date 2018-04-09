
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
noabx_normalized_reads <- 'data/read_mapping/noabx_normalized_metaT.tsv'
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
noabx_normalized_reads <- read.delim(noabx_normalized_reads, sep='\t', header=TRUE, row.names=6)
cef_normalized_reads <- read.delim(cef_normalized_reads, sep='\t', header=TRUE, row.names=7)
cef_normalized_reads$cef_630_metaT_reads <- NULL
clinda_normalized_reads <- read.delim(clinda_normalized_reads, sep='\t', header=TRUE, row.names=7)
clinda_normalized_reads$clinda_630_metaT_reads <- NULL
strep_normalized_reads <- read.delim(strep_normalized_reads, sep='\t', header=TRUE, row.names=7)
strep_normalized_reads$strep_630_metaT_reads <- NULL

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
# Rarefy to the same levels
sub_level <- min(sum(round(noabx_normalized_reads$noabx_mock_metaT_reads)),
                 sum(round(strep_normalized_reads$strep_mock_metaT_reads)),
                 sum(round(cef_normalized_reads$cef_mock_metaT_reads)),
                 sum(round(clinda_normalized_reads$clinda_mock_metaT_reads)))
noabx_normalized_reads$noabx_mock_metaT_reads <- as.vector(rrarefy(round(noabx_normalized_reads$noabx_mock_metaT_reads), sample=sub_level))
strep_normalized_reads$strep_mock_metaT_reads <- as.vector(rrarefy(round(strep_normalized_reads$strep_mock_metaT_reads), sample=sub_level))
cef_normalized_reads$cef_mock_metaT_reads <- as.vector(rrarefy(round(cef_normalized_reads$cef_mock_metaT_reads), sample=sub_level))
clinda_normalized_reads$clinda_mock_metaT_reads <- as.vector(rrarefy(round(clinda_normalized_reads$clinda_mock_metaT_reads), sample=sub_level))
rm(sub_level)
# Remove those genes without an organism annotation
noabx_annotated <- noabx_normalized_reads[!is.na(noabx_normalized_reads$organism),]
cef_annotated <- cef_normalized_reads[!is.na(cef_normalized_reads$organism),]
clinda_annotated <- clinda_normalized_reads[!is.na(clinda_normalized_reads$organism),]
strep_annotated <- strep_normalized_reads[!is.na(strep_normalized_reads$organism),]
rm(cef_normalized_reads, clinda_normalized_reads, strep_normalized_reads, noabx_normalized_reads)
# Remove any C. difficile genes with transcription only in infected
cdiff_omit <- c('cdf','pdc','cdc','cdl','pdf')
noabx_annotated <- subset(noabx_annotated, !(noabx_annotated$organism %in% cdiff_omit))
cef_annotated <- subset(cef_annotated, !(cef_annotated$organism %in% cdiff_omit))
clinda_annotated <- subset(clinda_annotated, !(clinda_annotated$organism %in% cdiff_omit))
strep_annotated <- subset(strep_annotated, !(strep_annotated$organism %in% cdiff_omit))
rm(cdiff_omit)
# Remove non-microbial genes
mamm_omit <- c('fab','cfa','ggo','hgl','hsa','mcc','mdo','pon','aml',
               'ptr','rno','shr','ssc','aml','bta','cge','ecb',
               'pps','fca','mmu','oaa','gga','ola','acs','aga')
noabx_annotated <- subset(noabx_annotated, !(noabx_annotated$organism %in% mamm_omit))
strep_annotated <- subset(strep_annotated, !(strep_annotated$organism %in% mamm_omit))
cef_annotated <- subset(cef_annotated, !(cef_annotated$organism %in% mamm_omit))
clinda_annotated <- subset(clinda_annotated, !(clinda_annotated$organism %in% mamm_omit))
rm(mamm_omit)
# Merge with genus level taxonomic classifications
noabx_annotated$gene <- rownames(noabx_annotated)
noabx_annotated <- merge(kegg_tax, noabx_annotated, by.x='org_code', by.y='organism', all.y=TRUE)
strep_annotated$gene <- rownames(strep_annotated)
strep_annotated <- merge(kegg_tax, strep_annotated, by.x='org_code', by.y='organism', all.y=TRUE)
cef_annotated$gene <- rownames(cef_annotated)
cef_annotated <- merge(kegg_tax, cef_annotated, by.x='org_code', by.y='organism', all.y=TRUE)
clinda_annotated$gene <- rownames(clinda_annotated)
clinda_annotated <- merge(kegg_tax, clinda_annotated, by.x='org_code', by.y='organism', all.y=TRUE)
rm(kegg_tax)
# Reformat no antibiotics control
noabx_annotated$org_code <- NULL
noabx_annotated$gene <- NULL
noabx_annotated$ko <- NULL
noabx_annotated <- noabx_annotated[complete.cases(noabx_annotated), ]
# Subset into minority grouos from 16s
noabx_strep_minority <- subset(noabx_annotated, genus %in% strep_minority_genera)
noabx_cef_minority <- subset(noabx_annotated, genus %in% cef_minority_genera)
noabx_clinda_minority <- subset(noabx_annotated, genus %in% clinda_minority_genera)
rm(noabx_annotated)
# Aggregate by genes
noabx_strep_minority$genus <- NULL
noabx_strep_minority$pathways <- NULL
noabx_strep_minority <- aggregate(. ~ description, data=noabx_strep_minority, FUN=sum)
noabx_cef_minority$genus <- NULL
noabx_cef_minority$pathways <- NULL
noabx_cef_minority <- aggregate(. ~ description, data=noabx_cef_minority, FUN=sum)
noabx_clinda_minority$genus <- NULL
noabx_clinda_minority$pathways <- NULL
noabx_clinda_minority <- aggregate(. ~ description, data=noabx_clinda_minority, FUN=sum)

#-------------------------------------------------------------------------------------------------------------------------#

# Tanscriptomic analysis of minority taxa
# Subset by minority genera from 16S analysis
strep_minority <- subset(strep_annotated, genus %in% strep_minority_genera)
cef_minority <- subset(cef_annotated, genus %in% cef_minority_genera)
clinda_minority <- subset(clinda_annotated, genus %in% clinda_minority_genera)
rm(strep_annotated, cef_annotated, clinda_annotated,
   strep_minority_genera, cef_minority_genera, clinda_minority_genera)
# Remove unused columns
strep_pathways <- strep_minority[,c(5,7)]
strep_pathways$pathways <- vapply(strsplit(as.character(strep_pathways$pathways),':'), `[`, 1, FUN.VALUE=character(1))
strep_pathways$pathways <- gsub('_', ' ', strep_pathways$pathways)
strep_minority$org_code <- NULL
strep_minority$ko <- NULL
strep_minority$genus <- NULL
strep_minority$pathways <- NULL
strep_minority$gene <- NULL
cef_pathways <- cef_minority[,c(5,7)]
cef_pathways$pathways <- vapply(strsplit(as.character(cef_pathways$pathways),':'), `[`, 1, FUN.VALUE=character(1))
cef_pathways$pathways <- gsub('_', ' ', cef_pathways$pathways)
cef_minority$org_code <- NULL
cef_minority$ko <- NULL
cef_minority$genus <- NULL
cef_minority$pathways <- NULL
cef_minority$gene <- NULL
clinda_pathways <- clinda_minority[,c(5,7)]
clinda_pathways$pathways <- vapply(strsplit(as.character(clinda_pathways$pathways),':'), `[`, 1, FUN.VALUE=character(1))
clinda_pathways$pathways <- gsub('_', ' ', clinda_pathways$pathways)
clinda_minority$org_code <- NULL
clinda_minority$ko <- NULL
clinda_minority$genus <- NULL
clinda_minority$pathways <- NULL
clinda_minority$gene <- NULL
# Agregate each treatment by gene description
strep_minority <- aggregate(. ~ description, data=strep_minority, FUN=sum)
strep_minority <- strep_minority[complete.cases(strep_minority), ]
cef_minority <- aggregate(. ~ description, data=cef_minority, FUN=sum)
cef_minority <- cef_minority[complete.cases(cef_minority), ]
clinda_minority <- aggregate(. ~ description, data=clinda_minority, FUN=sum)
clinda_minority <- clinda_minority[complete.cases(clinda_minority), ]
# Merge each treatment with no antibiotic controls
strep_minority <- merge(strep_minority, noabx_strep_minority, by='description')
cef_minority <- merge(cef_minority, noabx_cef_minority, by='description')
clinda_minority <- merge(clinda_minority, noabx_clinda_minority, by='description')
rm(noabx_strep_minority, noabx_cef_minority, noabx_clinda_minority)
# Reassociate pathway information
strep_minority <- merge(strep_minority, strep_pathways, by='description')
strep_minority <- strep_minority[complete.cases(strep_minority), ]
strep_minority <- subset(strep_minority, !strep_minority$pathways %in% c('Ribosome','Metabolic pathways'))
cef_minority <- merge(cef_minority, cef_pathways, by='description')
cef_minority <- cef_minority[complete.cases(cef_minority), ]
cef_minority <- subset(cef_minority, !cef_minority$pathways %in% c('Ribosome','Metabolic pathways'))
clinda_minority <- merge(clinda_minority, clinda_pathways, by='description')
clinda_minority <- clinda_minority[complete.cases(clinda_minority), ]
clinda_minority <- subset(clinda_minority, !clinda_minority$pathways %in% c('Ribosome','Metabolic pathways'))
rm(strep_pathways, cef_pathways, clinda_pathways)
# Aggregate by pathway
strep_minority$description <- NULL
strep_minority <- aggregate(. ~ pathways, data=strep_minority, FUN=sum)
rownames(strep_minority) <- strep_minority$pathways
strep_minority$pathways <- NULL
colnames(strep_minority) <- c('mock', 'resistant')
cef_minority$description <- NULL
cef_minority <- aggregate(. ~ pathways, data=cef_minority, FUN=sum)
rownames(cef_minority) <- cef_minority$pathways
cef_minority$pathways <- NULL
colnames(cef_minority) <- c('mock', 'resistant')
clinda_minority$description <- NULL
clinda_minority <- aggregate(. ~ pathways, data=clinda_minority, FUN=sum)
rownames(clinda_minority) <- clinda_minority$pathways
clinda_minority$pathways <- NULL
colnames(clinda_minority) <- c('mock', 'resistant')
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
# Subset to top hits at level of minimum between groups
path_level <- as.numeric(min(nrow(strep_minority), nrow(cef_minority), nrow(clinda_minority)))
strep_minority <- strep_minority[1:path_level,]
cef_minority <- cef_minority[1:path_level,]
clinda_minority <- clinda_minority[1:path_level,]
rm(path_level)
# Calculate magnitube of transcriptional change in minority groups within each treatment
strep_difference <- sum(strep_minority$abs_diff)
cef_difference <- sum(cef_minority$abs_diff)
clinda_difference <- sum(clinda_minority$abs_diff)
# Remove absolute difference column
strep_minority$abs_diff <- NULL
cef_minority$abs_diff <- NULL
clinda_minority$abs_diff <- NULL
# Log2 transform
#strep_minority <- log2(strep_minority + 1)
#cef_minority <- log2(cef_minority + 1)
#clinda_minority <- log2(clinda_minority + 1)

#-------------------------------------------------------------------------------------------------------------------------#

# Get data ready for plotting
minority_pathways <- as.data.frame(rbind(clinda_minority, cef_minority, strep_minority))
minority_pathways <- as.matrix(t(minority_pathways))
treatment_colors <- c(rep(c(clinda_col,noabx_col), nrow(clinda_minority)), 
                      rep(c(cef_col,noabx_col), nrow(cef_minority)),
                      rep(c(strep_col,noabx_col), nrow(strep_minority)))
pathway_names <- c("Phosphotransferase system (PTS)","Starch and sucrose metabolism","Two-component system","Pentose phosphate pathway",
                   "Fructose and mannose metabolism","RNA degradation","Amino sugar and nucleotide sugar metabolism","Carbon metabolism",
                   "Histidine metabolism","ABC transporters","Bacterial secretion system","Glycolysis / Gluconeogenesis",
                   "Phosphotransferase system (PTS)","Ascorbate and aldarate metabolism","Flagellar assembly","Glyoxylate and dicarboxylate metabolism",
                   "Pyrimidine metabolism","Biosynthesis of amino acids","Peptidoglycan biosynthesis","Glycerophospholipid metabolism",
                   "Protein export","Amino sugar and nucleotide sugar metabolism","Glycerolipid metabolism","Purine metabolism",
                   "RNA degradation","Thiamine metabolism","Bacterial secretion system")

#-------------------------------------------------------------------------------------------------------------------------#

# Generate figure
pdf(file=plot_file, width=12, height=20)
par(mar=c(4,25,1,2), mgp=c(2.5, 0.75, 0), las=1, xaxs='i', yaxs='i')
barplot(minority_pathways, yaxt='n', xlim=c(0,1200), ylim=c(0,82), beside=TRUE, horiz=TRUE, 
        xlab='', ylab='', col=treatment_colors, cex.axis=1.2)
axis(side=2, labels=pathway_names, at=seq(2,80,3), tick=FALSE, cex.axis=1.3)
mtext(expression(paste('Metagenome-normalized cDNA Reads (',log[2],')')), side=1, padj=2.5, cex=1.3)
text(x=1050, y=c(28.2,26.8), labels=c('Colonized', 'Clear'), font=2, cex=1.2)
legend('topright', legend=c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin'), 
       pt.bg=c(noabx_col, strep_col,cef_col,clinda_col), pch=22, pt.cex=2.8, cex=1.6)
abline(h=27.5, lwd=2) # separates cleared and colonized
box()
dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)}
setwd(starting_dir)
rm(list=ls())
gc()
