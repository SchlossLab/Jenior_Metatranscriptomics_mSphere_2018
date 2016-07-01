
deps <- c('vegan', 'wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Define input file names
cef_pangenome_file <- '/media/mjenior/Data/mapping/metagenome/Cefoperazone.DNA_reads2metaG.all.pool.norm.remove.annotated.txt'
clinda_pangenome_file <- '/media/mjenior/Data/mapping/metagenome/Cefoperazone.DNA_reads2metaG.all.pool.norm.remove.annotated.txt'
strep_pangenome_file <- '/media/mjenior/Data/mapping/metagenome/Cefoperazone.DNA_reads2metaG.all.pool.norm.remove.annotated.txt'
conv_pangenome_file <- '/media/mjenior/Data/mapping/metagenome/Cefoperazone.DNA_reads2metaG.all.pool.norm.remove.annotated.txt'

cef_630_metatranscriptome_file <- '/media/mjenior/Data/mapping/metatranscriptome/remove_cdf/cefoperazone_630.RNA_reads2metaG.cdf.all.pool.norm.remove.annotated.txt'
cef_mock_metatranscriptome_file <- '/media/mjenior/Data/mapping/metatranscriptome/remove_cdf/cefoperazone_mock.RNA_reads2metaG.cdf.all.pool.norm.remove.annotated.txt'
clinda_630_metatranscriptome_file <- '/media/mjenior/Data/mapping/metatranscriptome/remove_cdf/cefoperazone_630.RNA_reads2metaG.cdf.all.pool.norm.remove.annotated.txt'
clinda_mock_metatranscriptome_file <- '/media/mjenior/Data/mapping/metatranscriptome/remove_cdf/cefoperazone_mock.RNA_reads2metaG.cdf.all.pool.norm.remove.annotated.txt'
strep_630_metatranscriptome_file <- '/media/mjenior/Data/mapping/metatranscriptome/remove_cdf/cefoperazone_630.RNA_reads2metaG.cdf.all.pool.norm.remove.annotated.txt'
strep_mock_metatranscriptome_file <- '/media/mjenior/Data/mapping/metatranscriptome/remove_cdf/cefoperazone_mock.RNA_reads2metaG.cdf.all.pool.norm.remove.annotated.txt'
conv_metatranscriptome_file <- '/media/mjenior/Data/mapping/metatranscriptome/remove_cdf/cefoperazone_630.RNA_reads2metaG.cdf.all.pool.norm.remove.annotated.txt'

# Load in data
cef_pangenome <- read.delim(metagenome_file, sep='\t', header=TRUE, row.names=4)
clinda_pangenome <- read.delim(metagenome_file, sep='\t', header=TRUE, row.names=4)
strep_pangenome <- read.delim(metagenome_file, sep='\t', header=TRUE, row.names=4)
conv_pangenome <- read.delim(metagenome_file, sep='\t', header=TRUE, row.names=4)
rm(cef_pangenome, clinda_pangenome, strep_pangenome, conv_pangenome)

cef_630_metatranscriptome <- read.delim(cef_630_metatranscriptome_file, sep='\t', header=TRUE, row.names=4)
cef_mock_metatranscriptome <- read.delim(cef_mock_metatranscriptome_file, sep='\t', header=TRUE, row.names=4)
rm(cef_630_metatranscriptome_file, cef_mock_metatranscriptome_file)
clinda_630_metatranscriptome <- read.delim(clinda_630_metatranscriptome_file, sep='\t', header=TRUE, row.names=4)
clinda_mock_metatranscriptome <- read.delim(clinda_mock_metatranscriptome_file, sep='\t', header=TRUE, row.names=4)
rm(clinda_630_metatranscriptome_file, clinda_mock_metatranscriptome_file)
strep_630_metatranscriptome <- read.delim(strep_630_metatranscriptome_file, sep='\t', header=TRUE, row.names=4)
strep_mock_metatranscriptome <- read.delim(strep_mock_metatranscriptome_file, sep='\t', header=TRUE, row.names=4)
rm(strep_630_metatranscriptome_file, strep_mock_metatranscriptome_file)
conv_metatranscriptome <- read.delim(conv_metatranscriptome_file, sep='\t', header=TRUE, row.names=4)
rm(conv_metatranscriptome_file)

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
metaG_size <- round(read_totals[1] * 0.9) # 8309039
metaT_size <- round(min(read_totals[2:3]) * 0.9) # 3375147
combined_mapping$metagenome <- t(rrarefy(combined_mapping$metagenome, sample=metaG_size))
combined_mapping$infected_metatranscriptome <- t(rrarefy(combined_mapping$infected_metatranscriptome, sample=metaT_size))
combined_mapping$mock_metatranscriptome <- t(rrarefy(combined_mapping$mock_metatranscriptome, sample=metaT_size))

# Eliminate genes with no metagenomic mappings
combined_mapping <- subset(combined_mapping, metagenome != 0)

# Normalize to metagenomic coverage
combined_mapping <- subset(combined_mapping, infected_metatranscriptome != 0 & mock_metatranscriptome != 0)
combined_mapping$infected_metatranscriptome <- combined_mapping$infected_metatranscriptome / combined_mapping$metagenome
combined_mapping$mock_metatranscriptome <- combined_mapping$mock_metatranscriptome / combined_mapping$metagenome
combined_mapping$metagenome <- NULL

# Filter lowly abundant categories from normalized data, then log10 transform the data
combined_mapping <- subset(combined_mapping, infected_metatranscriptome > 1 & mock_metatranscriptome > 1)
combined_mapping$infected_metatranscriptome <- log10(combined_mapping$infected_metatranscriptome)
combined_mapping$mock_metatranscriptome <- log10(combined_mapping$mock_metatranscriptome)
norm_max <- max(combined_mapping[,c(1,2)])

#-------------------------------------------------------------------------------------------------------------------------#

# Plotting Metatranscriptome VS Metatranscriptome

# Subset data to color specific groups 
# Highest resolution - Known C. difficile 630 carbon sources
glycolysis <- subset(combined_mapping, grepl('*Glycolysis_/_Gluconeogenesis*', combined_mapping$pathway_annotation))
amino_sugar <- subset(combined_mapping, grepl('*Amino_sugar*', combined_mapping$pathway_annotation))
galactose <- subset(combined_mapping, grepl('*Galactose*', combined_mapping$pathway_annotation))
fructose_mannose <- subset(combined_mapping, grepl('*Fructose_and_mannose*', combined_mapping$pathway_annotation))
starch_sucrose <- subset(combined_mapping, grepl('*Starch_and_sucrose*', combined_mapping$pathway_annotation))
glycan_degradation <- subset(combined_mapping, grepl('*glycan_degradation*', combined_mapping$pathway_annotation))

pentose <- subset(combined_mapping, grepl('*Pentose*', combined_mapping$pathway_annotation))
cdf_carbohydates <- rbind(glycolysis, amino_sugar, galactose, fructose_mannose, starch_sucrose, glycan_degradation)

amino_sugar <- rbind(amino_sugar, glycan_degradation)
sugars <- rbind(glycolysis, fructose_mannose, starch_sucrose)

proline <- subset(combined_mapping, grepl('*roline*', combined_mapping$pathway_annotation))

cysteine_methionine <- subset(combined_mapping, grepl('*Cysteine_and_methionine*', combined_mapping$pathway_annotation))
valine_leucine_isoleucine <- subset(combined_mapping, grepl('*Valine,_leucine_and_isoleucine*', combined_mapping$pathway_annotation))
glycine_serine_threonine <- subset(combined_mapping, grepl('*Glycine,_serine_and_threonine*', combined_mapping$pathway_annotation))
alanine_aspartate_glutamate <- subset(combined_mapping, grepl('*Alanine,_aspartate_and_glutamate*', combined_mapping$pathway_annotation))
cdf_amino_acids <- rbind(cysteine_methionine, valine_leucine_isoleucine, glycine_serine_threonine, alanine_aspartate_glutamate)

cdf_carbon_sources <- rbind(cdf_carbohydates, cdf_amino_acids)


butanoate <- subset(combined_mapping, grepl('*Butanoate*', combined_mapping$pathway_annotation))[,c(1:2)]
pentose_glucuronate <- subset(combined_mapping, grepl('*Pentose_and_glucuronate*', combined_mapping$pathway_annotation))[,c(1:2)]
phenylalanine_tyrosine_tryptophan <- subset(combined_mapping, grepl('*Phenylalanine,_tyrosine_and_tryptophan*', combined_mapping$pathway_annotation))[,c(1:2)]
tca_cycle <- subset(combined_mapping, grepl('*TCA_cycle*', combined_mapping$pathway_annotation))[,c(1:2)]
lipopolysaccharide <- subset(combined_mapping, grepl('*Lipopolysaccharide*', combined_mapping$pathway_annotation))[,c(1:2)]
oxidative_phosphorylation <- subset(combined_mapping, grepl('*Oxidative_phosphorylation*', combined_mapping$pathway_annotation))[,c(1:2)]
glyoxylate_dicarboxylate <- subset(combined_mapping, grepl('*Glyoxylate_and_dicarboxylate*', combined_mapping$pathway_annotation))[,c(1:2)]
streptomycin_biosynthesis <- subset(combined_mapping, grepl('*Streptomycin_biosynthesis*', combined_mapping$pathway_annotation))[,c(1:2)]
polyketide_sugar_biosynthesis <- subset(combined_mapping, grepl('*Polyketide_sugar_unit_biosynthesis*', combined_mapping$pathway_annotation))[,c(1:2)]
lysine <- subset(combined_mapping, grepl('*Lysine*', combined_mapping$pathway_annotation))[,c(1:2)]
arginine_proline <- subset(combined_mapping, grepl('*Arginine_and_proline*', combined_mapping$pathway_annotation))[,c(1:2)]
pyruvate <- subset(combined_mapping, grepl('*Pyruvate*', combined_mapping$pathway_annotation))[,c(1:2)]
peptidoglycan <- subset(combined_mapping, grepl('*Peptidoglycan*', combined_mapping$pathway_annotation))[,c(1:2)]
glycan_degradation <- subset(combined_mapping, grepl('*glycan_degradation*', combined_mapping$pathway_annotation))[,c(1:2)]
secondary_bile <- subset(combined_mapping, grepl('*Secondary_bile*', combined_mapping$pathway_annotation))[,c(1:2)]
# More general scale
energy_metabolism <- subset(combined_mapping, grepl('*Energy_metabolism*', combined_mapping$pathway_annotation))[,c(1:2)]
carbohydrate_metabolism <- subset(combined_mapping, grepl('*Carbohydrate_metabolism*', combined_mapping$pathway_annotation))[,c(1:2)]
amino_acid_metabolism <- subset(combined_mapping, grepl('*Amino_acid_metabolism*', combined_mapping$pathway_annotation))[,c(1:2)]
nucleotide_metabolism <- subset(combined_mapping, grepl('*Nucleotide_metabolism*', combined_mapping$pathway_annotation))[,c(1:2)]
lipid_metabolism <- subset(combined_mapping, grepl('*Lipid_metabolism*', combined_mapping$pathway_annotation))[,c(1:2)]
cofactors_and_vitamins <- subset(combined_mapping, grepl('*cofactors_and_vitamins*', combined_mapping$pathway_annotation))[,c(1:2)]
diverse_environments <- subset(combined_mapping, grepl('*diverse_environments*', combined_mapping$pathway_annotation))[,c(1:2)]
sulfur_metabolism <- subset(combined_mapping, grepl('*Sulfur_metabolism*', combined_mapping$pathway_annotation))[,c(1:2)]
secondary_metabolites <- subset(combined_mapping, grepl('*secondary_metabolites*', combined_mapping$pathway_annotation))[,c(1:2)]

#-------------------------------------------------------------------------------------------------------------------------#

# Define which pathway to plot and the ouput file name
pathway <- diverse_environments
point_color <- 'forestgreen'
#pathway_name <- 'Carbohydrate Metabolism'
plot_file <- '~/Desktop/test3.pdf'

# Plot it!
pdf(file=plot_file, width=10, height=9)
par(mar=c(4, 4, 1, 1), mgp=c(2.4,0.7,0))
plot(x=combined_mapping$infected_metatranscriptome, y=combined_mapping$mock_metatranscriptome, 
     xlim=c(0,4), ylim=c(0,4), pch=20, col='black', xaxt='n', yaxt='n', cex.label=1.5,
     xlab='Normalized cDNA Read Abundance (Infected)', ylab='Normalized cDNA Read Abundance (Mock)')
#legend('topleft', legend=pathway_name, bty='n', cex=1.5) 
segments(-2, -2, 8, 8, lwd=2, lty=2)
axis_labels <- parse(text=paste(rep(10,4), '^', seq(0,4,1), sep=''))
axis(side=1, at=c(0:4), axis_labels, tick=TRUE)
axis(side=2, at=c(0:4), axis_labels, tick=TRUE, las=1)
#points(pathway[,c(1,2)], cex=1.5, pch=21, bg=point_color, col='black')
 

points(sugars[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[1], col='black')
points(amino_sugar[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[2], col='black')
points(galactose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[3], col='black')
points(pentose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[4], col='black')
points(valine_leucine_isoleucine[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling2")[4], col='black')
points(glycine_serine_threonine[,c(1,2)], cex=1.5, pch=21, bg='green3', col='black')
legend('bottomright', legend=c('Glucose/Fructose/Sucrose/Mannose', 'Amino sugars', 'Galactose', 'Pentose', 'Valine/Leucine/Isoleucine', 'Glycine/Serine/Threonine'), 
       pt.bg=c(wes_palette("Darjeeling")[1], wes_palette("Darjeeling")[2], wes_palette("Darjeeling")[3], wes_palette("Darjeeling")[4], wes_palette("Darjeeling2")[4], 'green3'), 
       pch=21, cex=1.2, pt.cex=1.5, bty='n')
dev.off()







# Subset cef treatment to zoom in on mock upregulation
mock_enriched <- subset(combined_mapping, infected_metatranscriptome < 1.2 & mock_metatranscriptome > 1.8)
par(mar=c(4, 4, 1, 1), mgp=c(2.4,0.7,0))
plot(x=mock_enriched$infected_metatranscriptome, y=mock_enriched$mock_metatranscriptome, 
     xlim=c(0,1.2), ylim=c(1.8,3), pch=20, col='gray25', xaxt='n', yaxt='n',
     xlab='Normalized cDNA Read Abundance (Infected)', ylab='Normalized cDNA Read Abundance (Mock)')
#legend('topleft', legend=pathway_name, bty='n', cex=1.5) 
abline(a=0, b=1, lwd=2)

xaxis_labels <- parse(text=paste(rep(10,4), '^', seq(0,1.2,0.2), sep=''))
yaxis_labels <- parse(text=paste(rep(10,4), '^', seq(1.8,3,0.2), sep=''))
axis(side=1, at=seq(0,1.2,0.2), xaxis_labels, tick=TRUE)
axis(side=2, at=seq(1.8,3,0.2), yaxis_labels, tick=TRUE, las=1)
points(pathway[,c(1,2)], cex=1.5, pch=21, bg=point_color, col='black')

abline(h=1.3, lwd=2)
abline(v=1.8, lwd=2)


# Subset cef treatment to zoom in on mock upregulation
mock_enriched <- subset(combined_mapping, infected_metatranscriptome < 1.2 & mock_metatranscriptome > 1.8)
par(mar=c(4, 4, 1, 1), mgp=c(2.4,0.7,0))
plot(x=mock_enriched$infected_metatranscriptome, y=mock_enriched$mock_metatranscriptome, 
     xlim=c(0,1.2), ylim=c(1.8,3), pch=20, col='gray25', xaxt='n', yaxt='n',
     xlab='Normalized cDNA Read Abundance (Infected)', ylab='Normalized cDNA Read Abundance (Mock)')
#legend('topleft', legend=pathway_name, bty='n', cex=1.5) 
abline(a=0, b=1, lwd=2)

xaxis_labels <- parse(text=paste(rep(10,4), '^', seq(0,1.2,0.2), sep=''))
yaxis_labels <- parse(text=paste(rep(10,4), '^', seq(1.8,3,0.2), sep=''))
axis(side=1, at=seq(0,1.2,0.2), xaxis_labels, tick=TRUE)
axis(side=2, at=seq(1.8,3,0.2), yaxis_labels, tick=TRUE, las=1)
points(pathway[,c(1,2)], cex=1.5, pch=21, bg=point_color, col='black')




