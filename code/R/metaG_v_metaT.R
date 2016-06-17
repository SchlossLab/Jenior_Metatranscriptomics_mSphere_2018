deps <- c('vegan');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Define input file names
metagenome_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/gene_mapping/metagenome/Cefoperazone.DNA_reads2metaG.all.pool.norm.remove.annotated.txt'
metatranscriptome1_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/gene_mapping/metatranscriptome/remove_cdf/cefoperazone_630.RNA_reads2metaG.cdf.all.pool.norm.remove.annotated.txt'
metatranscriptome2_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/data/gene_mapping/metatranscriptome/remove_cdf/cefoperazone_mock.RNA_reads2metaG.cdf.all.pool.norm.remove.annotated.txt'

# Load in data
metagenome <- read.delim(metagenome_file, sep='\t', header=TRUE, row.names=4)
metatranscriptome1 <- read.delim(metatranscriptome1_file, sep='\t', header=TRUE, row.names=4)
metatranscriptome2 <- read.delim(metatranscriptome2_file, sep='\t', header=TRUE, row.names=4)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data for merging
metagenome$gene_annotation <- NULL
metagenome$evalue <- NULL
metagenome$KEGG_ortholog <- NULL
metagenome$pathway_annotation <- NULL
metatranscriptome1$gene_annotation <- NULL
metatranscriptome1$evalue <- NULL
metatranscriptome1$KEGG_ortholog <- NULL
metatranscriptome1$pathway_annotation <- NULL

# Merge tables
combined_mapping <- merge(metagenome, metatranscriptome1, by='row.names')
rownames(combined_mapping) <- combined_mapping$Row.names
combined_mapping$Row.names <- NULL
combined_mapping <- merge(combined_mapping, metatranscriptome2, by='row.names')
rownames(combined_mapping) <- combined_mapping$Row.names
combined_mapping$Row.names <- NULL
colnames(combined_mapping) <- c('metagenome', 'infected_metatranscriptome', 'mock_metatranscriptome', 'gene_annotation', 'evalue', 'KEGG_ortholog', 'pathway_annotation')
combined_mapping$gene_annotation <- NULL
combined_mapping$evalue <- NULL
combined_mapping$KEGG_ortholog <- NULL

# Rarefy mappings to be equal within sequencing type
read_totals <- colSums(combined_mapping[,c(1:3)])
metaG_size <- round(read_totals[1]*0.9)
metaT_size <- round(min(read_totals[2:3])*0.9)
combined_mapping$metagenome <- t(rrarefy(combined_mapping$metagenome, sample=metaG_size))
combined_mapping$infected_metatranscriptome <- t(rrarefy(combined_mapping$infected_metatranscriptome, sample=metaT_size))
combined_mapping$mock_metatranscriptome <- t(rrarefy(combined_mapping$mock_metatranscriptome, sample=metaT_size))

# Eliminate genes with no metagenomic mappings
combined_mapping <- subset(combined_mapping, metagenome != 0)

# Log10 transform mappings
combined_mapping[combined_mapping == 0] <- 1
combined_mapping$metagenome <- log10(combined_mapping$metagenome)
combined_mapping$infected_metatranscriptome <- log10(combined_mapping$infected_metatranscriptome)
combined_mapping$mock_metatranscriptome <- log10(combined_mapping$mock_metatranscriptome)

#-------------------------------------------------------------------------------------------------------------------------#

# Plotting Metatranscriptome VS Metagenome

# Subset data to color specific groups 
# Finer scale  -  need to add sugar alcohols and stickland fermentation
amino_sugar <- subset(combined_mapping, grepl('*Amino_sugar*', combined_mapping$pathway_annotation))[,c(1:3)]
alanine_aspartate_glutamate <- subset(combined_mapping, grepl('*Alanine,_aspartate_and_glutamate*', combined_mapping$pathway_annotation))[,c(1:3)]
glycolysis_gluconeogenesis <- subset(combined_mapping, grepl('*Glycolysis_/_Gluconeogenesis*', combined_mapping$pathway_annotation))[,c(1:3)]
galactose <- subset(combined_mapping, grepl('*Galactose*', combined_mapping$pathway_annotation))[,c(1:3)]
cysteine_methionine <- subset(combined_mapping, grepl('*Cysteine_and_methionine*', combined_mapping$pathway_annotation))[,c(1:3)]
fructose_mannose <- subset(combined_mapping, grepl('*Fructose_and_mannose*', combined_mapping$pathway_annotation))[,c(1:3)]
aline_leucine_isoleucine <- subset(combined_mapping, grepl('*Valine,_leucine_and_isoleucine*', combined_mapping$pathway_annotation))[,c(1:3)]
glycine_serine_threonine <- subset(combined_mapping, grepl('*Glycine,_serine_and_threonine*', combined_mapping$pathway_annotation))[,c(1:3)]
starch_sucrose <- subset(combined_mapping, grepl('*Starch_and_sucrose*', combined_mapping$pathway_annotation))[,c(1:3)]
butanoate <- subset(combined_mapping, grepl('*Butanoate*', combined_mapping$pathway_annotation))[,c(1:3)]
pentose_glucuronate <- subset(combined_mapping, grepl('*Pentose_and_glucuronate*', combined_mapping$pathway_annotation))[,c(1:3)]
phenylalanine_tyrosine_tryptophan <- subset(combined_mapping, grepl('*Phenylalanine,_tyrosine_and_tryptophan*', combined_mapping$pathway_annotation))[,c(1:3)]
tca_cycle <- subset(combined_mapping, grepl('*TCA_cycle*', combined_mapping$pathway_annotation))[,c(1:3)]
lipopolysaccharide <- subset(combined_mapping, grepl('*Lipopolysaccharide*', combined_mapping$pathway_annotation))[,c(1:3)]
oxidative_phosphorylation <- subset(combined_mapping, grepl('*Oxidative_phosphorylation*', combined_mapping$pathway_annotation))[,c(1:3)]
glyoxylate_dicarboxylate <- subset(combined_mapping, grepl('*Glyoxylate_and_dicarboxylate*', combined_mapping$pathway_annotation))[,c(1:3)]
streptomycin_biosynthesis <- subset(combined_mapping, grepl('*Streptomycin_biosynthesis*', combined_mapping$pathway_annotation))[,c(1:3)]
polyketide_sugar_biosynthesis <- subset(combined_mapping, grepl('*Polyketide_sugar_unit_biosynthesis*', combined_mapping$pathway_annotation))[,c(1:3)]
lysine <- subset(combined_mapping, grepl('*Lysine*', combined_mapping$pathway_annotation))[,c(1:3)]
arginine_proline <- subset(combined_mapping, grepl('*Arginine_and_proline*', combined_mapping$pathway_annotation))[,c(1:3)]
pyruvate <- subset(combined_mapping, grepl('*Pyruvate*', combined_mapping$pathway_annotation))[,c(1:3)]
peptidoglycan <- subset(combined_mapping, grepl('*Peptidoglycan*', combined_mapping$pathway_annotation))[,c(1:3)]
glycan_degradation <- subset(combined_mapping, grepl('*glycan_degradation*', combined_mapping$pathway_annotation))[,c(1:3)]
secondary_bile <- subset(combined_mapping, grepl('*Secondary_bile*', combined_mapping$pathway_annotation))[,c(1:3)]

# More general scale
energy_metabolism <- subset(combined_mapping, grepl('*Energy_metabolism*', combined_mapping$pathway_annotation))[,c(1:3)]
carbohydrate_metabolism <- subset(combined_mapping, grepl('*Carbohydrate_metabolism*', combined_mapping$pathway_annotation))[,c(1:3)]
amino_acid_metabolism <- subset(combined_mapping, grepl('*Amino_acid_metabolism*', combined_mapping$pathway_annotation))[,c(1:3)]
nucleotide_metabolism <- subset(combined_mapping, grepl('*Nucleotide_metabolism*', combined_mapping$pathway_annotation))[,c(1:3)]
lipid_metabolism <- subset(combined_mapping, grepl('*Lipid_metabolism*', combined_mapping$pathway_annotation))[,c(1:3)]
cofactors_and_vitamins <- subset(combined_mapping, grepl('*cofactors_and_vitamins*', combined_mapping$pathway_annotation))[,c(1:3)]
diverse_environments <- subset(combined_mapping, grepl('*diverse_environments*', combined_mapping$pathway_annotation))[,c(1:3)]
sulfur_metabolism <- subset(combined_mapping, grepl('*Sulfur_metabolism*', combined_mapping$pathway_annotation))[,c(1:3)]
secondary_metabolites <- subset(combined_mapping, grepl('*secondary_metabolites*', combined_mapping$pathway_annotation))[,c(1:3)]

#-------------------------------------------------------------------------------------------------------------------------#

# Define which pathway to plot and the ouput file name
pathway <- secondary_bile
point_color <- 'darkorange2'
fig_title <- 'Secondary Bile Acid Synthesis'
plot_file <- '~/Desktop/Repositories/Jenior_Transcriptomics_2015/results/figures/cefoperazone.stretomycin.metaTvsmetaG.pdf'


# Plot the most interesting pathway changes

#pdf(file=plot_file, width=25, height=14)
layout(matrix(1:2, nrow=1, byrow=T))
par(mar=c(4, 4, 3, 1), mgp=c(2.4,0.7,0))
plot(x=combined_mapping$metagenome, y=combined_mapping$mock_metatranscriptome, 
     xlim=c(0,6), ylim=c(0,6), pch=20, col='gray25', xaxt='n', yaxt='n',
     xlab='Normalized Metagenomic Read Abundance', ylab='Normalized Metatranscriptomic Read Abundance')
segments(-2, -2, 8, 8, lwd=2)
legend('topleft', legend='Mock Infected', bty='n') 
axis_labels <- parse(text=paste(rep(10,7), '^', seq(0,6,1), sep=''))
axis(side=1, at=c(0:6), axis_labels, tick=TRUE)
axis(side=2, at=c(0:6), axis_labels, tick=TRUE, las=1)
points(pathway[,c(1,3)], cex=1.5, pch=21, bg=point_color, col='black')
par(mar=c(4, 4, 3, 1) + 0.1, mgp=c(2.3,0.7,0))
plot(x=combined_mapping$metagenome, y=combined_mapping$infected_metatranscriptome, 
     xlim=c(0,6), ylim=c(0,6), pch=20, col='gray25', xaxt='n', yaxt='n',
     xlab='Normalized Metagenomic Read Abundance', ylab='Normalized Metatranscriptomic Read Abundance')
segments(-2, -2, 8, 8, lwd=2)
plot_label <- substitute(paste(italic('C. diffile'), ' 630 Infected'))
legend('topleft', legend=plot_label, bty='n') 
axis_labels <- parse(text=paste(rep(10,7), '^', seq(0,6,1), sep=''))
axis(side=1, at=c(0:6), axis_labels, tick=TRUE)
axis(side=2, at=c(0:6), axis_labels, tick=TRUE, las=1)
points(pathway[,c(1,2)], cex=1.5, pch=21, bg=point_color, col='black')
mtext(fig_title, side = 3, line = -2, outer = TRUE, cex=2)
#dev.off()

