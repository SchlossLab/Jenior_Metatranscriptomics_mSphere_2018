
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
plot_file <- 'results/figures/figure_6.pdf'
legend_file <- 'results/figures/figure_6_legend.pdf'

#-------------------------------------------------------------------------------------------------------------------------#

# Read in data

# Metadata
metadata <- read.delim(metadata, sep='\t', header=TRUE, row.names=1)

# Normalized Metatranscriptomes
noabx_normalized_reads <- read.delim(noabx_normalized_reads, sep='\t', header=TRUE, row.names=6)
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
genus_tax_agg_shared$cage <- NULL
genus_tax_agg_shared$mouse <- NULL
genus_tax_agg_shared$gender <- NULL
genus_tax_agg_shared$type <- NULL
genus_tax_agg_shared$susceptibility <- NULL
genus_tax_agg_shared$infection <- NULL
abx_genus_tax_agg_shared <- aggregate(. ~ abx, data=genus_tax_agg_shared, FUN=max)
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
                 sum(round(strep_normalized_reads$strep_630_metaT_reads)),
                 sum(round(strep_normalized_reads$strep_mock_metaT_reads)),
                 sum(round(cef_normalized_reads$cef_630_metaT_reads)),
                 sum(round(cef_normalized_reads$cef_mock_metaT_reads)),
                 sum(round(clinda_normalized_reads$clinda_630_metaT_reads)),
                 sum(round(clinda_normalized_reads$clinda_mock_metaT_reads)))
noabx_normalized_reads$noabx_mock_metaT_reads <- as.vector(rrarefy(round(noabx_normalized_reads$noabx_mock_metaT_reads), sample=sub_level))
strep_normalized_reads$strep_630_metaT_reads <- as.vector(rrarefy(round(strep_normalized_reads$strep_630_metaT_reads), sample=sub_level))
strep_normalized_reads$strep_mock_metaT_reads <- as.vector(rrarefy(round(strep_normalized_reads$strep_mock_metaT_reads), sample=sub_level))
cef_normalized_reads$cef_630_metaT_reads <- as.vector(rrarefy(round(cef_normalized_reads$cef_630_metaT_reads), sample=sub_level))
cef_normalized_reads$cef_mock_metaT_reads <- as.vector(rrarefy(round(cef_normalized_reads$cef_mock_metaT_reads), sample=sub_level))
clinda_normalized_reads$clinda_630_metaT_reads <- as.vector(rrarefy(round(clinda_normalized_reads$clinda_630_metaT_reads), sample=sub_level))
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
noabx_annotated <- subset(noabx_annotated, noabx_annotated$pathways != 'Naphthalene_degradation:Xenobiotics_biodegradation_and_metabolism')
strep_annotated <- subset(strep_annotated, !(strep_annotated$organism %in% mamm_omit))
strep_annotated <- subset(strep_annotated, strep_annotated$pathways != 'Naphthalene_degradation:Xenobiotics_biodegradation_and_metabolism')
cef_annotated <- subset(cef_annotated, !(cef_annotated$organism %in% mamm_omit))
cef_annotated <- subset(cef_annotated, cef_annotated$pathways != 'Naphthalene_degradation:Xenobiotics_biodegradation_and_metabolism')
clinda_annotated <- subset(clinda_annotated, !(clinda_annotated$organism %in% mamm_omit))
clinda_annotated <- subset(clinda_annotated, clinda_annotated$pathways != 'Naphthalene_degradation:Xenobiotics_biodegradation_and_metabolism')
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
# Format untreated control
noabx_annotated$pathways <- vapply(strsplit(as.character(noabx_annotated$pathways),':'), `[`, 1, FUN.VALUE=character(1))
noabx_annotated$pathways <- gsub('_', ' ', noabx_annotated$pathways)
noabx_annotated$org_code <- NULL
noabx_annotated$gene <- NULL
noabx_annotated$ko <- NULL

#-------------------------------------------------------------------------------------------------------------------------#

# Find transcription levels of minority taxa
# Curate minority taxa that are at higher taxonomic level than genus
strep_minority_genera <- c(strep_minority_genera,"Cellulosilyticum",'Butyrivibrio','Dorea',"Lactococcus","Alistipes","Odoribacter","Staphylococcus","Streptococcus","Roseburia","Ruminococcus","Porphyromonas","Faecalibacterium")
strep_minority_genera <- unique(strep_minority_genera)
cef_minority_genera <- gsub('Shigella', 'Escherichia', cef_minority_genera)
cef_minority_genera <- c(cef_minority_genera, "Porphyromonas","Ruminococcus","Clostridium","Lactococcus","Odoribacter","Prevotella")
cef_minority_genera <- unique(cef_minority_genera)
clinda_minority_genera <- c(clinda_minority_genera,"Cellulosilyticum","Unclassified_Clostridiales","Clostridium","Lactococcus",'Butyrivibrio','Dorea','Roseburia','Bacteroides',"Prevotella","Ruminococcus","Porphyromonas","Faecalibacterium")
clinda_minority_genera <- unique(clinda_minority_genera)
# Subset untreated metatranscriptome to each treatment minority by 16S
noabx_strep_minority <- subset(noabx_annotated, genus %in% strep_minority_genera)
noabx_strep_minority$genus <- NULL
noabx_strep_minority <- noabx_strep_minority[complete.cases(noabx_strep_minority), ]
noabx_cef_minority <- subset(noabx_annotated, genus %in% cef_minority_genera)
noabx_cef_minority$genus <- NULL
noabx_cef_minority <- noabx_cef_minority[complete.cases(noabx_cef_minority), ]
noabx_clinda_minority <- subset(noabx_annotated, genus %in% clinda_minority_genera)
noabx_clinda_minority$genus <- NULL
noabx_clinda_minority <- noabx_clinda_minority[complete.cases(noabx_clinda_minority), ]
rm(noabx_annotated)
# Subset by minority genera from 16S analysis
strep_minority <- subset(strep_annotated, genus %in% strep_minority_genera)
strep_minority <- strep_minority[complete.cases(strep_minority), ]
cef_minority <- subset(cef_annotated, genus %in% cef_minority_genera)
cef_minority <- cef_minority[complete.cases(cef_minority), ]
clinda_minority <- subset(clinda_annotated, genus %in% clinda_minority_genera)
clinda_minority <- clinda_minority[complete.cases(clinda_minority), ]
rm(strep_annotated, cef_annotated, clinda_annotated,
   strep_minority_genera, cef_minority_genera, clinda_minority_genera)
# Reformat columns
strep_minority$pathways <- vapply(strsplit(as.character(strep_minority$pathways),':'), `[`, 1, FUN.VALUE=character(1))
strep_minority$pathways <- gsub('_', ' ', strep_minority$pathways)
strep_minority$org_code <- NULL
strep_minority$ko <- NULL
strep_minority$genus <- NULL
strep_minority$gene <- NULL
cef_minority$pathways <- vapply(strsplit(as.character(cef_minority$pathways),':'), `[`, 1, FUN.VALUE=character(1))
cef_minority$pathways <- gsub('_', ' ', cef_minority$pathways)
cef_minority$org_code <- NULL
cef_minority$ko <- NULL
cef_minority$genus <- NULL
cef_minority$gene <- NULL
clinda_minority$pathways <- vapply(strsplit(as.character(clinda_minority$pathways),':'), `[`, 1, FUN.VALUE=character(1))
clinda_minority$pathways <- gsub('_', ' ', clinda_minority$pathways)
clinda_minority$org_code <- NULL
clinda_minority$ko <- NULL
clinda_minority$genus <- NULL
clinda_minority$gene <- NULL
# Remove uninformative pathways and genes
noabx_strep_minority <- subset(noabx_strep_minority, !noabx_strep_minority$pathways %in% c('Ribosome','Metabolic pathways'))
noabx_cef_minority <- subset(noabx_cef_minority, !noabx_cef_minority$pathways %in% c('Ribosome','Metabolic pathways'))
noabx_clinda_minority <- subset(noabx_clinda_minority, !noabx_clinda_minority$pathways %in% c('Ribosome','Metabolic pathways'))
strep_minority <- subset(strep_minority, !strep_minority$pathways %in% c('Ribosome','Metabolic pathways'))
cef_minority <- subset(cef_minority, !cef_minority$pathways %in% c('Ribosome','Metabolic pathways'))
clinda_minority <- subset(clinda_minority, !clinda_minority$pathways %in% c('Ribosome','Metabolic pathways'))
noabx_strep_minority <- subset(noabx_strep_minority, !noabx_strep_minority$description %in% c('hypothetical_protein','putative'))
noabx_cef_minority <- subset(noabx_cef_minority, !noabx_cef_minority$description %in% c('hypothetical_protein','putative'))
noabx_clinda_minority <- subset(noabx_clinda_minority, !noabx_clinda_minority$description %in% c('hypothetical_protein','putative'))
strep_minority <- subset(strep_minority, !strep_minority$description %in% c('hypothetical_protein','putative'))
cef_minority <- subset(cef_minority, !cef_minority$description %in% c('hypothetical_protein','putative'))
clinda_minority <- subset(clinda_minority, !clinda_minority$description %in% c('hypothetical_protein','putative'))
# Aggregate each treatment by pathways or genes
# Pathways
noabx_strep_minority_pathways <- aggregate(. ~ pathways, data=noabx_strep_minority, FUN=sum)
noabx_strep_minority_pathways$description <- NULL
noabx_cef_minority_pathways <- aggregate(. ~ pathways, data=noabx_cef_minority, FUN=sum)
noabx_cef_minority_pathways$description <- NULL
noabx_clinda_minority_pathways <- aggregate(. ~ pathways, data=noabx_clinda_minority, FUN=sum)
noabx_clinda_minority_pathways$description <- NULL
strep_minority_pathways <- aggregate(. ~ pathways, data=strep_minority, FUN=sum)
strep_minority_pathways$description <- NULL
cef_minority_pathways <- aggregate(. ~ pathways, data=cef_minority, FUN=sum)
cef_minority_pathways$description <- NULL
clinda_minority_pathways <- aggregate(. ~ pathways, data=clinda_minority, FUN=sum)
clinda_minority_pathways$description <- NULL
# Genes
noabx_strep_minority$pathways <- NULL
noabx_strep_minority_genes <- aggregate(. ~ description, data=noabx_strep_minority, FUN=sum)
noabx_cef_minority$pathways <- NULL
noabx_cef_minority_genes <- aggregate(. ~ description, data=noabx_cef_minority, FUN=sum)
noabx_clinda_minority$pathways <- NULL
noabx_clinda_minority_genes <- aggregate(. ~ description, data=noabx_clinda_minority, FUN=sum)
strep_minority$pathways <- NULL
strep_minority_genes <- aggregate(. ~ description, data=strep_minority, FUN=sum)
cef_minority$pathways <- NULL
cef_minority_genes <- aggregate(. ~ description, data=cef_minority, FUN=sum)
clinda_minority$pathways <- NULL
clinda_minority_genes <- aggregate(. ~ description, data=clinda_minority, FUN=sum)
rm(clinda_minority, cef_minority, strep_minority,
   noabx_clinda_minority, noabx_cef_minority, noabx_strep_minority)
# Merge treatment and control minority groups
strep_minority_pathways <- merge(strep_minority_pathways, noabx_strep_minority_pathways, by='pathways')
cef_minority_pathways <- merge(cef_minority_pathways, noabx_cef_minority_pathways, by='pathways')
clinda_minority_pathways <- merge(clinda_minority_pathways, noabx_clinda_minority_pathways, by='pathways')
strep_minority_genes <- merge(strep_minority_genes, noabx_strep_minority_genes, by='description')
cef_minority_genes <- merge(cef_minority_genes, noabx_cef_minority_genes, by='description')
clinda_minority_genes <- merge(clinda_minority_genes, noabx_clinda_minority_genes, by='description')
rm(noabx_strep_minority_pathways, noabx_cef_minority_pathways, noabx_clinda_minority_pathways,
   noabx_strep_minority_genes, noabx_cef_minority_genes, noabx_clinda_minority_genes)
# Set row names
rownames(strep_minority_pathways) <- strep_minority_pathways$pathways
strep_minority_pathways$pathways <- NULL
colnames(strep_minority_pathways) <- c('infected','mock','resistant')
rownames(cef_minority_pathways) <- cef_minority_pathways$pathways
cef_minority_pathways$pathways <- NULL
colnames(cef_minority_pathways) <- c('infected','mock','resistant')
rownames(clinda_minority_pathways) <- clinda_minority_pathways$pathways
clinda_minority_pathways$pathways <- NULL
colnames(clinda_minority_pathways) <- c('infected','mock','resistant')
rownames(strep_minority_genes) <- strep_minority_genes$description
strep_minority_genes$description <- NULL
colnames(strep_minority_genes) <- c('infected','mock','resistant')
rownames(cef_minority_genes) <- cef_minority_genes$description
cef_minority_genes$description <- NULL
colnames(cef_minority_genes) <- c('infected','mock','resistant')
rownames(clinda_minority_genes) <- clinda_minority_genes$description
clinda_minority_genes$description <- NULL
colnames(clinda_minority_genes) <- c('infected','mock','resistant')

#-------------------------------------------------------------------------------------------------------------------------#

# Tanscriptomic analysis of minority taxa
# Find correct overlaps by clearance
colonized_pathways <- intersect(rownames(strep_minority_pathways), rownames(cef_minority_pathways))
colonized_genes <- intersect(rownames(strep_minority_genes), rownames(cef_minority_genes))
cleared_pathways <- rownames(clinda_minority_pathways)
cleared_genes <- rownames(clinda_minority_genes)
diff_pathways <- setdiff(colonized_pathways, cleared_pathways)
diff_pathways <- intersect(diff_pathways, colonized_pathways)
diff_genes <- setdiff(colonized_genes, cleared_genes)
diff_genes <- intersect(diff_genes, colonized_genes)
rm(colonized_pathways, cleared_pathways, colonized_genes, cleared_genes)
# Subset treatments
strep_minority_pathways <- subset(strep_minority_pathways, rownames(strep_minority_pathways) %in% diff_pathways)
cef_minority_pathways <- subset(cef_minority_pathways, rownames(cef_minority_pathways) %in% diff_pathways)
clinda_minority_pathways <- subset(clinda_minority_pathways, !rownames(clinda_minority_pathways) %in% diff_pathways)
strep_minority_genes <- subset(strep_minority_genes, rownames(strep_minority_genes) %in% diff_genes)
cef_minority_genes <- subset(cef_minority_genes, rownames(cef_minority_genes) %in% diff_genes)
clinda_minority_genes <- subset(clinda_minority_genes, !rownames(clinda_minority_genes) %in% diff_genes)
rm(diff_pathways, diff_genes)
rm(clinda_minority_pathways, cef_minority_pathways, strep_minority_pathways)
# Find degree of change caused by infection
strep_minority_genes$abs_diff <- abs(strep_minority_genes$infected - strep_minority_genes$mock)
strep_minority_genes <- subset(strep_minority_genes, abs_diff != 0)
cef_minority_genes$abs_diff <- abs(cef_minority_genes$infected - cef_minority_genes$mock)
cef_minority_genes <- subset(cef_minority_genes, abs_diff != 0)
clinda_minority_genes$abs_diff <- abs(clinda_minority_genes$infected - clinda_minority_genes$mock)
clinda_minority_genes <- subset(clinda_minority_genes, abs_diff != 0)

# Remove genes with very little transcriptional shift
sub_level <- round(sum(strep_minority_genes$abs_diff) * 0.02)
strep_minority_genes <- subset(strep_minority_genes, strep_minority_genes$abs_diff >= sub_level)
sub_level <- round(sum(cef_minority_genes$abs_diff) * 0.02)
cef_minority_genes <- subset(cef_minority_genes, cef_minority_genes$abs_diff >= sub_level)
sub_level <- round(sum(clinda_minority_genes$abs_diff) * 0.02)
clinda_minority_genes <- subset(clinda_minority_genes, clinda_minority_genes$abs_diff >= sub_level)
rm(sub_level)

# Sort tables by overall expression across genes/pathways
#strep_minority_pathways <- strep_minority_pathways[order(rowSums(strep_minority_pathways)),]
#cef_minority_pathways <- cef_minority_pathways[order(rowSums(cef_minority_pathways)),]
#clinda_minority_pathways <- clinda_minority_pathways[order(rowSums(clinda_minority_pathways)),]
strep_minority_genes <- strep_minority_genes[order(strep_minority_genes$abs_diff),]
strep_minority_genes$abs_diff <- NULL
cef_minority_genes <- cef_minority_genes[order(cef_minority_genes$abs_diff),]
cef_minority_genes$abs_diff <- NULL
clinda_minority_genes <- clinda_minority_genes[order(clinda_minority_genes$abs_diff),]
clinda_minority_genes$abs_diff <- NULL

# Transform data
#strep_minority_pathways <- log2(strep_minority_pathways + 1)
#cef_minority_pathways <- log2(cef_minority_pathways + 1)
#clinda_minority_pathways <- log2(clinda_minority_pathways + 1)
strep_minority_genes <- log2(strep_minority_genes + 1)
cef_minority_genes <- log2(cef_minority_genes + 1)
clinda_minority_genes <- log2(clinda_minority_genes + 1)

# Reformat gene names
strep_minority_genes$product <- rownames(strep_minority_genes)
strep_minority_genes$product <- capitalize(strep_minority_genes$product)
strep_minority_genes$product <- gsub('_', ' ', strep_minority_genes$product)
strep_minority_genes$product <- c("Molecular chaperone DnaK [A]","Glycerol kinase [B]",
                                  "UDP-galactose 4-epimerase [C,D]","Galactokinase [C,D]",
                                  "PTS transporter subunit IIC [E]","Phospho-beta-glucosidase [E,F]")
strep_minority_genes <- aggregate(. ~ product, data=strep_minority_genes, FUN=sum)
rownames(strep_minority_genes) <- strep_minority_genes$product
strep_minority_genes$product <- NULL
cef_minority_genes$product <- rownames(cef_minority_genes)
cef_minority_genes$product <- capitalize(cef_minority_genes$product)
cef_minority_genes$product <- gsub('_', ' ', cef_minority_genes$product)
cef_minority_genes$product <- c("Glucose-6-P isomerase [D,E,F]","Glycogen phosphorylase [E]",
                                "Uridylyltransferase [C,D]","Phosphoglycerate kinase [F]",
                                "Molecular chaperone DnaK [A]","ssDNA-binding protein [G]",
                                "PTS transporter subunit [E]","Predicted transcriptional regulator",
                                "Pyruvate kinase [F]","Galactokinase [C,D]","Galactokinase [C,D]",
                                "Glycoside hydrolase [D]")
cef_minority_genes <- subset(cef_minority_genes, product != 'Predicted transcriptional regulator')
cef_minority_genes <- aggregate(. ~ product, data=cef_minority_genes, FUN=sum)
rownames(cef_minority_genes) <- cef_minority_genes$product
cef_minority_genes$product <- NULL
clinda_minority_genes$product <- rownames(clinda_minority_genes)
clinda_minority_genes$product <- capitalize(clinda_minority_genes$product)
clinda_minority_genes$product <- gsub('_', ' ', clinda_minority_genes$product)
clinda_minority_genes$product <- c("Formate acetyltransferase [H]","Arginase [I]",
                                   "DNA-binding response regulator [G]","1-Phosphofructokinase [J]",
                                   "Glycosyl hydrolase [D]","Iron ABC transporter permease",
                                   "Sensor histidine kinase","Deoxyribose-phosphate aldolase [K]",
                                   "PTS transporter subunit IIA [E]","Glycerol-3-P dehydrogenase [B]",
                                   "Lactose/cellobiose IIC component [C,E]","PTS transporter subunit IID [E]",
                                   "6-phospho-beta-glucosidase [E,F]","Putative nucleotidyltransferase",
                                   "Glycosyl hydrolase [D]")
clinda_minority_genes <- subset(clinda_minority_genes, product != 'Putative nucleotidyltransferase')
clinda_minority_genes <- aggregate(. ~ product, data=clinda_minority_genes, FUN=sum)
rownames(clinda_minority_genes) <- clinda_minority_genes$product
clinda_minority_genes$product <- NULL

#-------------------------------------------------------------------------------------------------------------------------#

# Get data ready for plotting
# Pathways
#minority_pathways <- as.data.frame(rbind(clinda_minority_pathways, cef_minority_pathways, strep_minority_pathways))
#minority_pathways <- as.matrix(t(minority_pathways))
#treatment_colors <- c(rep(c(clinda_col,clinda_col,'gray70'), nrow(clinda_minority_pathways)), 
#                      rep(c(cef_col,cef_col,'gray70'), nrow(cef_minority_pathways)),
#                      rep(c(strep_col,strep_col,'gray70'), nrow(strep_minority_pathways)))
#pathway_names <- c("Base excision repair","Glycerophospholipid metabolism","Pentose phosphate pathway","Fructose and mannose metabolism",
#                   "Starch and sucrose metabolism","Glycolysis / Gluconeogenesis","Microbial metabolism in diverse environments","Two-component system",
#                   "Phosphotransferase system (PTS)","ABC transporters","Tyrosine metabolism","Cysteine and methionine metabolism",
#                   "Homologous recombination","Pyrimidine metabolism*","Purine metabolism*","Oxidative phosphorylation*","Aminoacyl-tRNA biosynthesis*",
#                   "Propanoate metabolism","Nucleotide excision repair","Chloroalkane and chloroalkene degradation","Citrate cycle (TCA cycle)",
#                   "Peptidoglycan biosynthesis","Glycerolipid metabolism","Methane metabolism","Bacterial secretion system","Purine metabolism*",
#                   "Pyrimidine metabolism*","Oxidative phosphorylation*","Aminoacyl-tRNA biosynthesis*")
#rm(clinda_minority_pathways, cef_minority_pathways, strep_minority_pathways)

# Genes
minority_genes <- as.data.frame(rbind(clinda_minority_genes, cef_minority_genes, strep_minority_genes))
minority_genes <- as.matrix(t(minority_genes))
treatment_colors <- c(rep(c(clinda_col,clinda_col,'gray70'), nrow(clinda_minority_genes)), 
                      rep(c(cef_col,cef_col,'gray70'), nrow(cef_minority_genes)),
                      rep(c(strep_col,strep_col,'gray70'), nrow(strep_minority_genes)))
gene_names <- c("1-Phosphofructokinase [J]","6-phospho-beta-glucosidase [E,F]",
                "Arginase [I]","Deoxyribose-phosphate aldolase [K]",
                "DNA-binding response regulator [G]","Formate acetyltransferase [H]",
                "Glycerol-3-P dehydrogenase [B]","Glycosyl hydrolase [D]",
                "Iron ABC transporter permease","Lactose/cellobiose IIC component [C,E]",
                "PTS transporter subunit IIA [E]","PTS transporter subunit IID [E]",
                "Sensor histidine kinase","Galactokinase [C,D]**",
                "Glucose-6-P isomerase [D,E,F]","Glycogen phosphorylase [E]",
                "Glycoside hydrolase [D]","Molecular chaperone DnaK [A]**",
                "Phosphoglycerate kinase [F]","PTS transporter subunit [E]",
                "Pyruvate kinase [F]","ssDNA-binding protein [G]",
                "Uridylyltransferase [C,D]","Galactokinase [C,D]**",
                "Glycerol kinase [B]","Molecular chaperone DnaK [A]**",
                "Phospho-beta-glucosidase [E,F]","PTS transporter subunit IIC [E]",
                "UDP-galactose 4-epimerase [C,D]")
rm(clinda_minority_genes, cef_minority_genes, strep_minority_genes)

pathways <- c('A. RNA degradation','B. Glycerolipid metabolism','C. Galactose metabolism',
              'D. Amino sugar and nucleotide sugar metabolism','E. Starch and sucrose metabolism',
              'F. Glycolysis / Gluconeogenesis','G. DNA replication','H. Butanoate / Propanoate metabolism',
              'I. Arginine and proline metabolism','J. Fructose and mannose metabolism','K. Pentose phosphate pathway')

#-------------------------------------------------------------------------------------------------------------------------#

# Generate figure
pdf(file=plot_file, width=15, height=20)
par(mar=c(5,26,1,2), mgp=c(2.5, 0.75, 0), las=1, xaxs='i', yaxs='i')
barplot(minority_genes, xaxt='n', yaxt='n', xlim=c(0,10), ylim=c(0,88), beside=TRUE, horiz=TRUE, width=0.75,
        xlab='', ylab='', col=treatment_colors, cex.axis=1.2, angle=20, density=c(NA,20,NA))
axis(1, at=seq(0,10,2), label=seq(0,10,2), cex.axis=1.6)
minor.ticks.axis(1, 10, mn=0, mx=10)
axis(side=2, labels=gene_names, at=seq(2,86,3), tick=FALSE, cex.axis=1.7)
mtext(expression(paste('Normalized cDNA Reads (',log[2],')')), side=1, padj=2.5, cex=1.8)
abline(h=seq(3.4,84.4,3), lty=2)
legend('topright', legend=c('No Antibiotics','Streptomycin-pretreated','Cefoperazone-pretreated','Clindamycin-pretreated'), 
       pt.bg=c('gray70',strep_col,cef_col,clinda_col), pch=22, pt.cex=3.3, cex=2, box.lwd=1, box.col="black", bg="white")
legend(x=6.6, y=78.25, legend=c('630-infected','Mock-infected'), angle=20, density=c(20,NA), 
       pt.cex=3.3, cex=2, box.lwd=1, box.col="black", bg="white")
box()
dev.off()

# Pathway legend
pdf(file=legend_file, width=40, height=10)
par(mar=c(0,0,0,0))
plot(0, type='n', xlim=c(-25,25), ylim=c(-6,6), pch=20, xaxt='n', yaxt='n', xlab='', ylab='')
legend('center', legend=pathways, pt.cex=0, cex=3, ncol=3, box.lwd=2)
dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)}
setwd(starting_dir)
rm(list=ls())
gc()
