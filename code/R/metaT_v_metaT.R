
deps <- c('vegan', 'wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Define input file names
cef_metagenome_file <- '/media/mjenior/Jenior\ HD/data/mapping/pangenome/metagenome/Cefoperazone.DNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_metagenome_file <- '/media/mjenior/Jenior\ HD/data/mapping/pangenome/metagenome/Cefoperazone.DNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_metagenome_file <- '/media/mjenior/Jenior\ HD/data/mapping/pangenome/metagenome/Cefoperazone.DNA_reads2pangenome.all.norm.remove.annotated.txt'
conv_metagenome_file <- '/media/mjenior/Jenior\ HD/data/mapping/pangenome/metagenome/Cefoperazone.DNA_reads2pangenome.all.norm.remove.annotated.txt'

cef_630_metatranscriptome_file <- '/media/mjenior/Jenior\ HD/data/mapping/pangenome/metatranscriptome/cefoperazone_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'
cef_mock_metatranscriptome_file <- '/media/mjenior/Jenior\ HD/data/mapping/pangenome/metatranscriptome/cefoperazone_mock.RNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_630_metatranscriptome_file <- '/media/mjenior/Jenior\ HD/data/mapping/pangenome/metatranscriptome/cefoperazone_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_mock_metatranscriptome_file <- '/media/mjenior/Jenior\ HD/data/mapping/pangenome/metatranscriptome/cefoperazone_mock.RNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_630_metatranscriptome_file <- '/media/mjenior/Jenior\ HD/data/mapping/pangenome/metatranscriptome/cefoperazone_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_mock_metatranscriptome_file <- '/media/mjenior/Jenior\ HD/data/mapping/pangenome/metatranscriptome/cefoperazone_mock.RNA_reads2pangenome.all.norm.remove.annotated.txt'
conv_metatranscriptome_file <- '/media/mjenior/Jenior\ HD/data/mapping/pangenome/metatranscriptome/cefoperazone_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'

# Load in data
# Metagenomes
cef_metagenome <- read.delim(cef_metagenome_file, sep='\t', header=FALSE, row.names=1)
colnames(cef_metagenome) <- c('reads', 'ko', 'gene', 'pathway')
clinda_metagenome <- read.delim(clinda_metagenome_file, sep='\t', header=FALSE, row.names=1)
colnames(clinda_metagenome) <- c('reads', 'ko', 'gene', 'pathway')
strep_metagenome <- read.delim(strep_metagenome_file, sep='\t', header=FALSE, row.names=1)
colnames(strep_metagenome) <- c('reads', 'ko', 'gene', 'pathway')
conv_metagenome <- read.delim(conv_metagenome_file, sep='\t', header=FALSE, row.names=1)
colnames(conv_metagenome) <- c('reads', 'ko', 'gene', 'pathway')
rm(cef_metagenome_file, clinda_metagenome_file, strep_metagenome_file, conv_metagenome_file)
# Metatranscriptomes
cef_630_metatranscriptome <- read.delim(cef_630_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(cef_630_metatranscriptome) <- c('reads', 'ko', 'gene', 'pathway')
cef_mock_metatranscriptome <- read.delim(cef_mock_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(cef_mock_metatranscriptome) <- c('reads', 'ko', 'gene', 'pathway')
rm(cef_630_metatranscriptome_file, cef_mock_metatranscriptome_file)
clinda_630_metatranscriptome <- read.delim(clinda_630_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(clinda_630_metatranscriptome) <- c('reads', 'ko', 'gene', 'pathway')
clinda_mock_metatranscriptome <- read.delim(clinda_mock_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(clinda_mock_metatranscriptome) <- c('reads', 'ko', 'gene', 'pathway')
rm(clinda_630_metatranscriptome_file, clinda_mock_metatranscriptome_file)
strep_630_metatranscriptome <- read.delim(strep_630_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(strep_630_metatranscriptome) <- c('reads', 'ko', 'gene', 'pathway')
strep_mock_metatranscriptome <- read.delim(strep_mock_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(strep_mock_metatranscriptome) <- c('reads', 'ko', 'gene', 'pathway')
rm(strep_630_metatranscriptome_file, strep_mock_metatranscriptome_file)
conv_metatranscriptome <- read.delim(conv_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(conv_metatranscriptome) <- c('reads', 'ko', 'gene', 'pathway')
rm(conv_metatranscriptome_file)

#-------------------------------------------------------------------------------------------------------------------------#

# Format metagenomic data for merging
cef_metagenome$ko <- NULL
cef_metagenome$gene <- NULL
cef_metagenome$pathway <- NULL
clinda_metagenome$ko <- NULL
clinda_metagenome$gene <- NULL
clinda_metagenome$pathway <- NULL
strep_metagenome$ko <- NULL
strep_metagenome$gene <- NULL
strep_metagenome$pathway <- NULL

# Merge metagenome tables
all_metagenome <- merge(cef_metagenome, clinda_metagenome, by='row.names')
rownames(all_metagenome) <- all_metagenome$Row.names
all_metagenome$Row.names <- NULL
all_metagenome <- merge(all_metagenome, strep_metagenome, by='row.names')
rownames(all_metagenome) <- all_metagenome$Row.names
all_metagenome$Row.names <- NULL
all_metagenome <- merge(all_metagenome, conv_metagenome, by='row.names')
rownames(all_metagenome) <- all_metagenome$Row.names
all_metagenome$Row.names <- NULL
colnames(all_metagenome) <- c('cefoperazone', 'clindamycin', 'streptomycin', 'conventional', 'ko', 'gene', 'pathway')
rm(cef_metagenome, clinda_metagenome, strep_metagenome, conv_metagenome)

#-------------------------------------------------------------------------------------------------------------------------#

# Format metatranscriptomic data for merging
cef_630_metatranscriptome$ko <- NULL
cef_630_metatranscriptome$gene <- NULL
cef_630_metatranscriptome$pathway <- NULL
cef_mock_metatranscriptome$ko <- NULL
cef_mock_metatranscriptome$gene <- NULL
cef_mock_metatranscriptome$pathway <- NULL
clinda_630_metatranscriptome$ko <- NULL
clinda_630_metatranscriptome$gene <- NULL
clinda_630_metatranscriptome$pathway <- NULL
clinda_mock_metatranscriptome$ko <- NULL
clinda_mock_metatranscriptome$gene <- NULL
clinda_mock_metatranscriptome$pathway <- NULL
strep_630_metatranscriptome$ko <- NULL
strep_630_metatranscriptome$gene <- NULL
strep_630_metatranscriptome$pathway <- NULL
strep_mock_metatranscriptome$ko <- NULL
strep_mock_metatranscriptome$gene <- NULL
strep_mock_metatranscriptome$pathway <- NULL

# Merge metatranscriptome tables
all_metatranscriptome <- merge(cef_630_metatranscriptome, cef_mock_metatranscriptome, by='row.names')
rownames(all_metatranscriptome) <- all_metatranscriptome$Row.names
all_metatranscriptome$Row.names <- NULL
all_metatranscriptome <- merge(all_metatranscriptome, clinda_630_metatranscriptome, by='row.names')
rownames(all_metatranscriptome) <- all_metatranscriptome$Row.names
all_metatranscriptome$Row.names <- NULL
all_metatranscriptome <- merge(all_metatranscriptome, clinda_mock_metatranscriptome, by='row.names')
rownames(all_metatranscriptome) <- all_metatranscriptome$Row.names
all_metatranscriptome$Row.names <- NULL
all_metatranscriptome <- merge(all_metatranscriptome, strep_630_metatranscriptome, by='row.names')
rownames(all_metatranscriptome) <- all_metatranscriptome$Row.names
all_metatranscriptome$Row.names <- NULL
all_metatranscriptome <- merge(all_metatranscriptome, strep_mock_metatranscriptome, by='row.names')
rownames(all_metatranscriptome) <- all_metatranscriptome$Row.names
all_metatranscriptome$Row.names <- NULL
all_metatranscriptome <- merge(all_metatranscriptome, conv_metatranscriptome, by='row.names')
rownames(all_metatranscriptome) <- all_metatranscriptome$Row.names
all_metatranscriptome$Row.names <- NULL
colnames(all_metatranscriptome) <- c('cefoperazone_630', 'cefoperazone_mock', 'clindamycin_630', 'clindamycin_mock', 'streptomycin_630', 'streptomycin_mock','conventional','ko', 'gene', 'pathway')
rm(cef_630_metatranscriptome, cef_mock_metatranscriptome, clinda_630_metatranscriptome, 
   clinda_mock_metatranscriptome, strep_630_metatranscriptome, strep_mock_metatranscriptome, conv_metatranscriptome)

#-------------------------------------------------------------------------------------------------------------------------#

# Determine subsample sizes
metagenome_totals <- colSums(all_metagenome[,c(1:4)])
metatranscriptome_totals <- colSums(all_metatranscriptome[,c(1:7)])
metaG_size <- round(min(metagenome_totals) * 0.9) # 
metaT_size <- round(min(metatranscriptome_totals) * 0.9) # 
rm(metagenome_totals, metatranscriptome_totals)

# Rarefy metagenomic data
all_metagenome$cefoperazone <- t(rrarefy(all_metagenome$cefoperazone, sample=metaG_size))
all_metagenome$clindamycin <- t(rrarefy(all_metagenome$clindamycin, sample=metaG_size))
all_metagenome$streptomycin <- t(rrarefy(all_metagenome$streptomycin, sample=metaG_size))
all_metagenome$conventional <- t(rrarefy(all_metagenome$conventional, sample=metaG_size))
rm(metaG_size)

# Eliminate genes with no metagenomic mappings
all_metagenome$rowSums <- rowSums(all_metagenome[,c(1:4)])
all_metagenome <- subset(all_metagenome, all_metagenome$rowSums != 0)
all_metagenome$rowSums <- NULL
   
# Rarefy metatranscriptomic data
all_metatranscriptome$cefoperazone_630 <- t(rrarefy(all_metatranscriptome$cefoperazone_630, sample=metaT_size))
all_metatranscriptome$cefoperazone_mock <- t(rrarefy(all_metatranscriptome$cefoperazone_mock, sample=metaT_size))
all_metatranscriptome$clindamycin_630 <- t(rrarefy(all_metatranscriptome$clindamycin_630, sample=metaT_size))
all_metatranscriptome$clindamycin_mock <- t(rrarefy(all_metatranscriptome$clindamycin_mock, sample=metaT_size))
all_metatranscriptome$streptomycin_630 <- t(rrarefy(all_metatranscriptome$streptomycin_630, sample=metaT_size))
all_metatranscriptome$streptomycin_mock <- t(rrarefy(all_metatranscriptome$streptomycin_mock, sample=metaT_size))
all_metatranscriptome$conventional <- t(rrarefy(all_metatranscriptome$conventional, sample=metaT_size))
rm(metaT_size)

# Normalize to metagenomic coverage
all_metatranscriptome$cefoperazone_630 <- all_metatranscriptome$cefoperazone_630 / all_metagenome$cefoperazone
all_metatranscriptome$cefoperazone_mock <- all_metatranscriptome$cefoperazone_mock / all_metagenome$cefoperazone
all_metatranscriptome$clindamycin_630 <- all_metatranscriptome$clindamycin_630 / all_metagenome$clindamycin
all_metatranscriptome$clindamycin_mock <- all_metatranscriptome$clindamycin_mock / all_metagenome$clindamycin
all_metatranscriptome$streptomycin_630 <- all_metatranscriptome$streptomycin_630 / all_metagenome$streptomycin
all_metatranscriptome$streptomycin_mock <- all_metatranscriptome$streptomycin_mock / all_metagenome$streptomycin
all_metatranscriptome$conventional <- all_metatranscriptome$conventional / all_metagenome$conventional
rm(all_metagenome)

# Separate treatment groups
cef_metatranscriptome <- all_metatranscriptome[,c(1,2,8:10)]
clinda_metatranscriptome <- all_metatranscriptome[,c(3,4,8:10)]
strep_metatranscriptome <- all_metatranscriptome[,c(5,6,8:10)]
conv_metatranscriptome <- all_metatranscriptome[,c(7:10)]
rm(all_metatranscriptome)

# Filter lowly abundant categories from normalized data, then log10 transform the data
cef_metatranscriptome$cefoperazone_630[cef_metatranscriptome$cefoperazone_630 == 0] <- 1
cef_metatranscriptome$cefoperazone_6mock[cef_metatranscriptome$cefoperazone_mock == 0] <- 1
cef_metatranscriptome[,c(1,2)] <- log10(cef_metatranscriptome[,c(1,2)])
clinda_metatranscriptome$clindamycin_630[clinda_metatranscriptome$clindamycin_630 == 0] <- 1
clinda_metatranscriptome$clindamycin_mock[clinda_metatranscriptome$clindamycin_mock == 0] <- 1
clinda_metatranscriptome[,c(1,2)] <- log10(clinda_metatranscriptome[,c(1,2)])
strep_metatranscriptome$streptomycin_630[strep_metatranscriptome$streptomycin_630 == 0] <- 1
strep_metatranscriptome$streptomycin_mock[strep_metatranscriptome$streptomycin_mock == 0] <- 1
strep_metatranscriptome[,c(1,2)] <- log10(strep_metatranscriptome[,c(1,2)])
conv_metatranscriptome$conventional[conv_metatranscriptome$conventional == 0] <- 1
conv_metatranscriptome[,1] <- log10(conv_metatranscriptome[,1])

#-------------------------------------------------------------------------------------------------------------------------#






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




