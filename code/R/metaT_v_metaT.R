
deps <- c('vegan', 'wesanderson');
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  }
  library(dep, verbose=FALSE, character.only=TRUE)
}

dist2d <- function(a, b, c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
} 


# Define input file names
cef_metagenome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Cefoperazone.DNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_metagenome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Clindamycin.DNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_metagenome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Streptomycin.DNA_reads2pangenome.all.norm.remove.annotated.txt'
conv_metagenome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Conventional.DNA_reads2pangenome.all.norm.remove.annotated.txt'

cef_630_metatranscriptome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/cefoperazone_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'
cef_mock_metatranscriptome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/cefoperazone_mock.RNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_630_metatranscriptome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/clindamycin_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_mock_metatranscriptome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/clindamycin_mock.RNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_630_metatranscriptome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/streptomycin_630.RNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_mock_metatranscriptome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/streptomycin_mock.RNA_reads2pangenome.all.norm.remove.annotated.txt'
conv_metatranscriptome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metatranscriptome/conventional.RNA_reads2pangenome.all.norm.remove.annotated.txt'




# Load in data
# Metagenomes
cef_metagenome <- read.delim(cef_metagenome_file, sep='\t', header=FALSE, row.names=1)
colnames(cef_metagenome) <- c('cefoperazone_metaG', 'ko', 'gene', 'pathway')
clinda_metagenome <- read.delim(clinda_metagenome_file, sep='\t', header=FALSE, row.names=1)
colnames(clinda_metagenome) <- c('clindamycin_metaG', 'ko', 'gene', 'pathway')
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
metaG_size <- round(min(metagenome_totals) * 0.9)
metaT_size <- round(min(metatranscriptome_totals) * 0.9)
optimal_size <- min(c(metaG_size, metaT_size))
rm(metagenome_totals, metatranscriptome_totals, metaG_size, metaT_size)

# Rarefy metagenomic data
all_metagenome$cefoperazone <- t(rrarefy(all_metagenome$cefoperazone, sample=optimal_size)) + 1
all_metagenome$clindamycin <- t(rrarefy(all_metagenome$clindamycin, sample=optimal_size)) + 1
all_metagenome$streptomycin <- t(rrarefy(all_metagenome$streptomycin, sample=optimal_size)) + 1
all_metagenome$conventional <- t(rrarefy(all_metagenome$conventional, sample=optimal_size)) + 1

# Rarefy metatranscriptomic data
all_metatranscriptome$cefoperazone_630 <- t(rrarefy(all_metatranscriptome$cefoperazone_630, sample=optimal_size)) + 1
all_metatranscriptome$cefoperazone_mock <- t(rrarefy(all_metatranscriptome$cefoperazone_mock, sample=optimal_size)) + 1
all_metatranscriptome$clindamycin_630 <- t(rrarefy(all_metatranscriptome$clindamycin_630, sample=optimal_size)) + 1
all_metatranscriptome$clindamycin_mock <- t(rrarefy(all_metatranscriptome$clindamycin_mock, sample=optimal_size)) + 1
all_metatranscriptome$streptomycin_630 <- t(rrarefy(all_metatranscriptome$streptomycin_630, sample=optimal_size)) + 1
all_metatranscriptome$streptomycin_mock <- t(rrarefy(all_metatranscriptome$streptomycin_mock, sample=optimal_size)) + 1
all_metatranscriptome$conventional_metaT <- t(rrarefy(all_metatranscriptome$conventional_metaT, sample=optimal_size)) + 1
rm(optimal_size)

# Merge metagenomes and metatranscriptomes
all_metagenome$ko <- NULL
all_metagenome$gene <- NULL
all_metagenome$pathway <- NULL
full_mapping <- merge(all_metagenome, all_metatranscriptome, by='row.names')
rownames(full_mapping) <- full_mapping$Row.names
full_mapping$Row.names <- NULL
colnames(full_mapping) <- c('cefoperazone_metaG', 'clindamycin_metaG', 'streptomycin_metaG', 'conventional_metaG', 'cefoperazone_630', 'cefoperazone_mock', 'clindamycin_630', 'clindamycin_mock', 'streptomycin_630', 'streptomycin_mock','conventional_metaT','ko', 'gene', 'pathway')
rm(all_metagenome, all_metatranscriptome)

# Normalize to metagenomic coverage
full_mapping$cefoperazone_630 <- full_mapping$cefoperazone_630 / full_mapping$cefoperazone_metaG
full_mapping$cefoperazone_mock <- full_mapping$cefoperazone_mock / full_mapping$cefoperazone_metaG
full_mapping$clindamycin_630 <- full_mapping$clindamycin_630 / full_mapping$clindamycin_metaG
full_mapping$clindamycin_mock <- full_mapping$clindamycin_mock / full_mapping$clindamycin_metaG
full_mapping$streptomycin_630 <- full_mapping$streptomycin_630 / full_mapping$streptomycin_metaG
full_mapping$streptomycin_mock <- full_mapping$streptomycin_mock / full_mapping$streptomycin_metaG
full_mapping$conventional_metaT <- full_mapping$conventional_metaT / full_mapping$conventional_metaG

# Separate treatment groups
cef_metatranscriptome <- full_mapping[,c(5,6,12:14)]
clinda_metatranscriptome <- full_mapping[,c(7,8,12:14)]
strep_metatranscriptome <- full_mapping[,c(9,10,12:14)]
rm(full_mapping)

# Log10 transform the data
cef_metatranscriptome[,c(1,2)] <- log10(cef_metatranscriptome[,c(1,2)])
clinda_metatranscriptome[,c(1,2)] <- log10(clinda_metatranscriptome[,c(1,2)])
strep_metatranscriptome[,c(1,2)] <- log10(strep_metatranscriptome[,c(1,2)])

# Remove groups with no transcription
cef_metatranscriptome <- subset(cef_metatranscriptome, cef_metatranscriptome$cefoperazone_630 != 0.0 & cef_metatranscriptome$cefoperazone_mock != 0.0)
clinda_metatranscriptome <- subset(clinda_metatranscriptome, clinda_metatranscriptome$clindamycin_630 != 0.0 & clinda_metatranscriptome$clindamycin_mock != 0.0)
strep_metatranscriptome <- subset(strep_metatranscriptome, strep_metatranscriptome$streptomycin_630 != 0.0 & strep_metatranscriptome$streptomycin_mock != 0.0)

#-------------------------------------------------------------------------------------------------------------------------#

# Calculate the distance of all points from x=y
# Will reveal which genes were most effected by c. diff colonization


dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
} 





good.dist <- sqrt((good.ord - typ.ord)^2 / 2)






#-------------------------------------------------------------------------------------------------------------------------#

# Define which pathway to plot and the ouput file name
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/figures/figure_3.pdf'
pdf(file=plot_file, width=10, height=10)
layout(matrix(c(1,2,
                3,4), 
              nrow=2, ncol=2, byrow = TRUE))

#-------------------------------------------------------------------------------------------------------------------------#

# Cefoperazone
par(mar=c(4, 4, 1, 1), mgp=c(2.4,0.7,0), xaxs='i', yaxs='i')
plot(x=cef_metatranscriptome$cefoperazone_630, y=cef_metatranscriptome$cefoperazone_mock, 
     xlim=c(0,4), ylim=c(0,4), pch=20, col='black', xaxt='n', yaxt='n', 
     xlab='Normalized cDNA Read Abundance (Infected)', ylab='Normalized cDNA Read Abundance (Mock)')
segments(-2, -2, 8, 8, lwd=2, lty=2)
axis_labels <- parse(text=paste(rep(10,4), '^', seq(0,4,1), sep=''))
axis(side=1, at=c(0:4), axis_labels, tick=TRUE)
axis(side=2, at=c(0:4), axis_labels, tick=TRUE, las=1)
points(cef_glycolysis[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('Rushmore')[5], col='black')
points(cef_amino_sugar[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('FantasticFox')[3], col='black')
points(cef_fructose_mannose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('Rushmore')[4], col='black')
points(cef_pentose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('Rushmore')[3], col='black')
points(cef_proline[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('FantasticFox')[4], col='black')
points(cef_glycine_serine_threonine[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('Zissou')[3], col='black')
legend('topleft', 'Cefoperazone', bty='n', cex=1.4)

#-------------------------------------------------------------------------------------------------------------------------#

# Clindamycin
par(mar=c(4, 4, 1, 1), mgp=c(2.4,0.7,0), xaxs='i', yaxs='i')
plot(x=clinda_metatranscriptome$clindamycin_630, y=clinda_metatranscriptome$clindamycin_mock, 
     xlim=c(0,4), ylim=c(0,4), pch=20, col='black', xaxt='n', yaxt='n', 
     xlab='Normalized cDNA Read Abundance (Infected)', ylab='Normalized cDNA Read Abundance (Mock)')
segments(-2, -2, 8, 8, lwd=2, lty=2)
axis(side=1, at=c(0:4), axis_labels, tick=TRUE)
axis(side=2, at=c(0:4), axis_labels, tick=TRUE, las=1)
points(clinda_glycolysis[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('Rushmore')[5], col='black')
points(clinda_amino_sugar[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('FantasticFox')[3], col='black')
points(clinda_fructose_mannose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('Rushmore')[4], col='black')
points(clinda_pentose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('Rushmore')[3], col='black')
points(clinda_proline[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('FantasticFox')[4], col='black')
points(clinda_glycine_serine_threonine[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('Zissou')[3], col='black')
legend('topleft', 'Clindamycin', bty='n', cex=1.4)

#-------------------------------------------------------------------------------------------------------------------------#

# Streptomycin
par(mar=c(4, 4, 1, 1), mgp=c(2.4,0.7,0), xaxs='i', yaxs='i')
plot(x=strep_metatranscriptome$streptomycin_630, y=strep_metatranscriptome$streptomycin_mock, 
     xlim=c(0,4), ylim=c(0,4), pch=20, col='black', xaxt='n', yaxt='n', 
     xlab='Normalized cDNA Read Abundance (Infected)', ylab='Normalized cDNA Read Abundance (Mock)')
segments(-2, -2, 8, 8, lwd=2, lty=2)
axis(side=1, at=c(0:4), axis_labels, tick=TRUE)
axis(side=2, at=c(0:4), axis_labels, tick=TRUE, las=1)
points(strep_glycolysis[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('Rushmore')[5], col='black')
points(strep_amino_sugar[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('FantasticFox')[3], col='black')
points(strep_fructose_mannose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('Rushmore')[4], col='black')
points(strep_pentose[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('Rushmore')[3], col='black')
points(strep_proline[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('FantasticFox')[4], col='black')
points(strep_glycine_serine_threonine[,c(1,2)], cex=1.5, pch=21, bg=wes_palette('Zissou')[3], col='black')
legend('topleft', 'Streptomycin', bty='n', cex=1.4)

#-------------------------------------------------------------------------------------------------------------------------#

plot(1, type="n", axes=F, xlab="", ylab="")
legend('center', legend=c('Glycolysis', 'Amino sugars', 'Fructose/Mannose','Pentose', 'Proline', 'Glycine/Serine/Threonine'), 
       pt.bg=c(wes_palette('Rushmore')[5], wes_palette('FantasticFox')[3], wes_palette('Rushmore')[4], wes_palette('Rushmore')[3], wes_palette('FantasticFox')[4], wes_palette('Zissou')[3]), 
       pch=21, cex=1.6, pt.cex=2.5)

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()
#rm(cef_metatranscriptome, clinda_metatranscriptome, strep_metatranscriptome, plot_file, axis_labels)
