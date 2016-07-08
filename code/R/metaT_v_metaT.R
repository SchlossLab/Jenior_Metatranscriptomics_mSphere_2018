
deps <- c('vegan', 'wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

# Define input file names
cef_metagenome_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/mapping/pangenome/Cefoperazone.DNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_metagenome_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/mapping/pangenome/Cefoperazone.DNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_metagenome_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/mapping/pangenome/Cefoperazone.DNA_reads2pangenome.all.norm.remove.annotated.txt'
conv_metagenome_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/mapping/pangenome/Cefoperazone.DNA_reads2pangenome.all.norm.remove.annotated.txt'

cef_630_metatranscriptome_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/mapping/pangenome/cefoperazone_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'
cef_mock_metatranscriptome_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/mapping/pangenome/cefoperazone_mock.RNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_630_metatranscriptome_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/mapping/pangenome/cefoperazone_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_mock_metatranscriptome_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/mapping/pangenome/cefoperazone_mock.RNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_630_metatranscriptome_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/mapping/pangenome/cefoperazone_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_mock_metatranscriptome_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/mapping/pangenome/cefoperazone_mock.RNA_reads2pangenome.all.norm.remove.annotated.txt'
conv_metatranscriptome_file <- '/Users/mattjenior/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/mapping/pangenome/cefoperazone_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'

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
colnames(all_metatranscriptome) <- c('cefoperazone_630', 'cefoperazone_mock', 'clindamycin_630', 'clindamycin_mock', 'streptomycin_630', 'streptomycin_mock','conventional_metaT','ko', 'gene', 'pathway')
rm(cef_630_metatranscriptome, cef_mock_metatranscriptome, clinda_630_metatranscriptome, 
   clinda_mock_metatranscriptome, strep_630_metatranscriptome, strep_mock_metatranscriptome, conv_metatranscriptome)

#-------------------------------------------------------------------------------------------------------------------------#

# Determine subsample sizes
metagenome_totals <- colSums(all_metagenome[,c(1:4)])
metatranscriptome_totals <- colSums(all_metatranscriptome[,c(1:7)])
metaG_size <- round(min(metagenome_totals) * 0.9) # 
metaT_size <- round(min(metatranscriptome_totals) * 0.9) # 
optimal_size <- min(c(metaG_size, metaT_size))
rm(metagenome_totals, metatranscriptome_totals, metaG_size, metaT_size)

# Rarefy metagenomic data
all_metagenome$cefoperazone <- t(rrarefy(all_metagenome$cefoperazone, sample=optimal_size))
all_metagenome$clindamycin <- t(rrarefy(all_metagenome$clindamycin, sample=optimal_size))
all_metagenome$streptomycin <- t(rrarefy(all_metagenome$streptomycin, sample=optimal_size))
all_metagenome$conventional <- t(rrarefy(all_metagenome$conventional, sample=optimal_size))

# Eliminate genes with no metagenomic mappings
all_metagenome$cefoperazone[all_metagenome$cefoperazone == 0] <- 1
all_metagenome$clindamycin[all_metagenome$clindamycin == 0] <- 1
all_metagenome$streptomycin[all_metagenome$streptomycin == 0] <- 1
all_metagenome$conventional[all_metagenome$conventional == 0] <- 1

# Rarefy metatranscriptomic data
all_metatranscriptome$cefoperazone_630 <- t(rrarefy(all_metatranscriptome$cefoperazone_630, sample=optimal_size))
all_metatranscriptome$cefoperazone_mock <- t(rrarefy(all_metatranscriptome$cefoperazone_mock, sample=optimal_size))
all_metatranscriptome$clindamycin_630 <- t(rrarefy(all_metatranscriptome$clindamycin_630, sample=optimal_size))
all_metatranscriptome$clindamycin_mock <- t(rrarefy(all_metatranscriptome$clindamycin_mock, sample=optimal_size))
all_metatranscriptome$streptomycin_630 <- t(rrarefy(all_metatranscriptome$streptomycin_630, sample=optimal_size))
all_metatranscriptome$streptomycin_mock <- t(rrarefy(all_metatranscriptome$streptomycin_mock, sample=optimal_size))
all_metatranscriptome$conventional <- t(rrarefy(all_metatranscriptome$conventional, sample=optimal_size))
rm(optimal_size)

# Merge metagenomes and metatranscriptomes
all_metagenome$ko <- NULL
all_metagenome$gene <- NULL
all_metagenome$pathway <- NULL
full_mapping <- merge(all_metagenome, all_metatranscriptome, by='row.names')
rm(all_metagenome, all_metatranscriptome)

# Normalize to metagenomic coverage
full_mapping$cefoperazone_630 <- full_mapping$cefoperazone_630 / full_mapping$cefoperazone
full_mapping$cefoperazone_mock <- full_mapping$cefoperazone_mock / full_mapping$cefoperazone
full_mapping$clindamycin_630 <- full_mapping$clindamycin_630 / full_mapping$clindamycin
full_mapping$clindamycin_mock <- full_mapping$clindamycin_mock / full_mapping$clindamycin
full_mapping$streptomycin_630 <- full_mapping$streptomycin_630 / full_mapping$streptomycin
full_mapping$streptomycin_mock <- full_mapping$streptomycin_mock / full_mapping$streptomycin
full_mapping$conventional_metaT <- full_mapping$conventional_metaT / full_mapping$conventional

# Separate treatment groups
cef_metatranscriptome <- full_mapping[,c(1,2,13:15)]
clinda_metatranscriptome <- full_mapping[,c(3,4,13:15)]
strep_metatranscriptome <- full_mapping[,c(5,6,13:15)]
conv_metatranscriptome <- full_mapping[,c(7,8,13:15)]
rm(full_mapping)

# Filter lowly abundant categories from normalized data, then log10 transform the data
cef_metatranscriptome$cefoperazone_630[cef_metatranscriptome$cefoperazone_630 < 1] <- 1
cef_metatranscriptome$cefoperazone_6mock[cef_metatranscriptome$cefoperazone_mock < 1] <- 1
cef_metatranscriptome[,c(1,2)] <- log10(cef_metatranscriptome[,c(1,2)])
clinda_metatranscriptome$clindamycin_630[clinda_metatranscriptome$clindamycin_630 < 1] <- 1
clinda_metatranscriptome$clindamycin_mock[clinda_metatranscriptome$clindamycin_mock < 1] <- 1
clinda_metatranscriptome[,c(1,2)] <- log10(clinda_metatranscriptome[,c(1,2)])
strep_metatranscriptome$streptomycin_630[strep_metatranscriptome$streptomycin_630 < 1] <- 1
strep_metatranscriptome$streptomycin_mock[strep_metatranscriptome$streptomycin_mock < 1] <- 1
strep_metatranscriptome[,c(1,2)] <- log10(strep_metatranscriptome[,c(1,2)])
conv_metatranscriptome$conventional[conv_metatranscriptome$conventional < 1] <- 1
conv_metatranscriptome[,1] <- log10(conv_metatranscriptome[,1])

#-------------------------------------------------------------------------------------------------------------------------#

# Subset data to color specific groups 
# Known C. difficile 630 carbon sources
cef_glycolysis <- subset(cef_metatranscriptome, grepl('*Glycolysis_/_Gluconeogenesis*', cef_metatranscriptome$pathway))
cef_amino_sugar <- subset(cef_metatranscriptome, grepl('*Amino_sugar*', cef_metatranscriptome$pathway))
cef_galactose <- subset(cef_metatranscriptome, grepl('*Galactose*', cef_metatranscriptome$pathway))
cef_fructose_mannose <- subset(cef_metatranscriptome, grepl('*Fructose_and_mannose*', cef_metatranscriptome$pathway))
cef_starch_sucrose <- subset(cef_metatranscriptome, grepl('*Starch_and_sucrose*', cef_metatranscriptome$pathway))
cef_glycan_degradation <- subset(cef_metatranscriptome, grepl('*glycan_degradation*', cef_metatranscriptome$pathway))
cef_pentose <- subset(cef_metatranscriptome, grepl('*Pentose*', cef_metatranscriptome$pathway))
cef_proline <- subset(cef_metatranscriptome, grepl('*roline*', cef_metatranscriptome$pathway))
cef_glycine_serine_threonine <- subset(cef_metatranscriptome, grepl('*Glycine,_serine_and_threonine*', cef_metatranscriptome$pathway))

clinda_glycolysis <- subset(clinda_metatranscriptome, grepl('*Glycolysis_/_Gluconeogenesis*', clinda_metatranscriptome$pathway))
clinda_amino_sugar <- subset(clinda_metatranscriptome, grepl('*Amino_sugar*', clinda_metatranscriptome$pathway))
clinda_galactose <- subset(clinda_metatranscriptome, grepl('*Galactose*', clinda_metatranscriptome$pathway))
clinda_fructose_mannose <- subset(clinda_metatranscriptome, grepl('*Fructose_and_mannose*', clinda_metatranscriptome$pathway))
clinda_starch_sucrose <- subset(clinda_metatranscriptome, grepl('*Starch_and_sucrose*', clinda_metatranscriptome$pathway))
clinda_glycan_degradation <- subset(clinda_metatranscriptome, grepl('*glycan_degradation*', clinda_metatranscriptome$pathway))
clinda_pentose <- subset(clinda_metatranscriptome, grepl('*Pentose*', clinda_metatranscriptome$pathway))
clinda_proline <- subset(clinda_metatranscriptome, grepl('*roline*', clinda_metatranscriptome$pathway))
clinda_glycine_serine_threonine <- subset(clinda_metatranscriptome, grepl('*Glycine,_serine_and_threonine*', clinda_metatranscriptome$pathway))

strep_glycolysis <- subset(strep_metatranscriptome, grepl('*Glycolysis_/_Gluconeogenesis*', strep_metatranscriptome$pathway))
strep_amino_sugar <- subset(strep_metatranscriptome, grepl('*Amino_sugar*', strep_metatranscriptome$pathway))
strep_galactose <- subset(strep_metatranscriptome, grepl('*Galactose*', strep_metatranscriptome$pathway))
strep_fructose_mannose <- subset(strep_metatranscriptome, grepl('*Fructose_and_mannose*', strep_metatranscriptome$pathway))
strep_starch_sucrose <- subset(strep_metatranscriptome, grepl('*Starch_and_sucrose*', strep_metatranscriptome$pathway))
strep_glycan_degradation <- subset(strep_metatranscriptome, grepl('*glycan_degradation*', strep_metatranscriptome$pathway))
strep_pentose <- subset(strep_metatranscriptome, grepl('*Pentose*', strep_metatranscriptome$pathway))
strep_proline <- subset(strep_metatranscriptome, grepl('*roline*', strep_metatranscriptome$pathway))
strep_glycine_serine_threonine <- subset(strep_metatranscriptome, grepl('*Glycine,_serine_and_threonine*', strep_metatranscriptome$pathway))

#-------------------------------------------------------------------------------------------------------------------------#

# Define which pathway to plot and the ouput file name
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/figures/figure_3.pdf'
pdf(file=plot_file, width=14, height=14)
layout(matrix(c(1,2,
                3,4), 
              nrow=2, ncol=2, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------#

# Cefoperazone
par(mar=c(4, 4, 1, 1), mgp=c(2.4,0.7,0))
plot(x=cef_metatranscriptome$cefoperazone_630, y=cef_metatranscriptome$cefoperazone_mock, 
     xlim=c(0,4), ylim=c(0,4), pch=20, col='black', xaxt='n', yaxt='n', 
     xlab='Normalized cDNA Read Abundance (Infected)', ylab='Normalized cDNA Read Abundance (Mock)')
segments(-2, -2, 8, 8, lwd=2, lty=2)
axis_labels <- parse(text=paste(rep(10,4), '^', seqs(0,4,1), sep=''))
axis(side=1, at=c(0:4), axis_labels, tick=TRUE)
axis(side=2, at=c(0:4), axis_labels, tick=TRUE, las=1)
points(cef_glycolysis[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[1], col='black')
points(cef_amino_sugar[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[2], col='black')
points(cef_galactose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[3], col='black')
points(cef_fructose_mannose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[4], col='black')
points(cef_pentose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling2")[4], col='black')
points(cef_proline[,c(1,2)], cex=1.5, pch=21, bg='green3', col='black')
points(cef_glycine_serine_threonine[,c(1,2)], cex=1.5, pch=21, bg='red', col='black')
legend('topleft', 'Cefoperazone')

#-------------------------------------------------------------------------------------------------------------------------#

# Clindamycin
par(mar=c(4, 4, 1, 1), mgp=c(2.4,0.7,0))
plot(x=clinda_metatranscriptome$clindamycin_630, y=clinda_metatranscriptome$clindamycin_mock, 
     xlim=c(0,4), ylim=c(0,4), pch=20, col='black', xaxt='n', yaxt='n', 
     xlab='Normalized cDNA Read Abundance (Infected)', ylab='Normalized cDNA Read Abundance (Mock)')
segments(-2, -2, 8, 8, lwd=2, lty=2)
axis_labels <- parse(text=paste(rep(10,4), '^', seq(0,4,1), sep=''))
axis(side=1, at=c(0:4), axis_labels, tick=TRUE)
axis(side=2, at=c(0:4), axis_labels, tick=TRUE, las=1)
points(clinda_glycolysis[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[1], col='black')
points(clinda_amino_sugar[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[2], col='black')
points(clinda_galactose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[3], col='black')
points(clinda_fructose_mannose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[4], col='black')
points(clinda_pentose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling2")[4], col='black')
points(clinda_proline[,c(1,2)], cex=1.5, pch=21, bg='green3', col='black')
points(clinda_glycine_serine_threonine[,c(1,2)], cex=1.5, pch=21, bg='red', col='black')
legend('topleft', 'Clindamycin')

#-------------------------------------------------------------------------------------------------------------------------#

# Streptomycin
par(mar=c(4, 4, 1, 1), mgp=c(2.4,0.7,0))
plot(x=strep_metatranscriptome$streptomycin_630, y=strep_metatranscriptome$streptomycin_mock, 
     xlim=c(0,4), ylim=c(0,4), pch=20, col='black', xaxt='n', yaxt='n', 
     xlab='Normalized cDNA Read Abundance (Infected)', ylab='Normalized cDNA Read Abundance (Mock)')
segments(-2, -2, 8, 8, lwd=2, lty=2)
axis_labels <- parse(text=paste(rep(10,4), '^', seq(0,4,1), sep=''))
axis(side=1, at=c(0:4), axis_labels, tick=TRUE)
axis(side=2, at=c(0:4), axis_labels, tick=TRUE, las=1)
points(strep_glycolysis[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[1], col='black')
points(strep_amino_sugar[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[2], col='black')
points(strep_galactose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[3], col='black')
points(strep_fructose_mannose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling")[4], col='black')
points(strep_pentose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette("Darjeeling2")[4], col='black')
points(strep_proline[,c(1,2)], cex=1.5, pch=21, bg='green3', col='black')
points(strep_glycine_serine_threonine[,c(1,2)], cex=1.5, pch=21, bg='red', col='black')
legend('topleft', 'Streptomycin')

#-------------------------------------------------------------------------------------------------------------------------#

plot(1, type="n", axes=F, xlab="", ylab="")
legend('center', legend=c('Glycolysis', 'Amino sugars', 'Galactose', 'Fructose/Mannose','Pentose', 'Proline', 'Glycine/Serine/Threonine'), 
       pt.bg=c(wes_palette("Darjeeling")[1], wes_palette("Darjeeling")[2], wes_palette("Darjeeling")[3], wes_palette("Darjeeling")[4], wes_palette("Darjeeling2")[4], 'green3', 'red'), 
       pch=21, cex=1.2, pt.cex=1.5, bty='n')

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()







# Subset cef treatment to zoom in on mock upregulation
mock_enriched <- subset(combined_mapping, infected_metatranscriptome < 1.2 & mock_metatranscriptome > 1.8)
par(mar=c(4, 4, 1, 1), mgp=c(2.4,0.7,0))
plot(x=mock_enriched$infected_metatranscriptome, y=mock_enriched$mock_metatranscriptome, 
     xlim=c(0,1.2), ylim=c(1.8,3), pch=20, col='gray25', xaxt='n', yaxt='n',
     xlab='Normalized cDNA Read Abundance (Infected)', ylab='Normalized cDNA Read Abundance (Mock)')
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




