
# Load dependencies
deps <- c('randomForest', 'vegan', 'shape')
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}
set.seed(6189)

clean_merge <- function(data_1, data_2){
  
  clean_merged <- merge(data_1, data_2, by = 'row.names')
  rownames(clean_merged) <- clean_merged$Row.names
  clean_merged$Row.names <- NULL
  
  return(clean_merged)
}

# Define input file names
# metagenomes
cef_metagenome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Cefoperazone.DNA_reads2pangenome.all.norm.remove.annotated.txt'
clinda_metagenome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Clindamycin.DNA_reads2pangenome.all.norm.remove.annotated.txt'
strep_metagenome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Streptomycin.DNA_reads2pangenome.all.norm.remove.annotated.txt'
conv_metagenome_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/read_mapping/metagenome/Conventional.DNA_reads2pangenome.all.norm.remove.annotated.txt'
# metatranscriptomes
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
colnames(cef_metagenome) <- c('cef_metaG_reads', 'ko', 'gene', 'pathway')
clinda_metagenome <- read.delim(clinda_metagenome_file, sep='\t', header=FALSE, row.names=1)
colnames(clinda_metagenome) <- c('clinda_metaG_reads', 'ko', 'gene', 'pathway')
strep_metagenome <- read.delim(strep_metagenome_file, sep='\t', header=FALSE, row.names=1)
colnames(strep_metagenome) <- c('strep_metaG_reads', 'ko', 'gene', 'pathway')
conv_metagenome <- read.delim(conv_metagenome_file, sep='\t', header=FALSE, row.names=1)
colnames(conv_metagenome) <- c('conv_metaG_reads', 'ko', 'gene', 'pathway')
rm(cef_metagenome_file, clinda_metagenome_file, strep_metagenome_file, conv_metagenome_file)
# Metatranscriptomes
cef_630_metatranscriptome <- read.delim(cef_630_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(cef_630_metatranscriptome) <- c('cef_630_metaT_reads', 'ko', 'gene', 'pathway')
cef_mock_metatranscriptome <- read.delim(cef_mock_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(cef_mock_metatranscriptome) <- c('cef_mock_metaT_reads', 'ko', 'gene', 'pathway')
rm(cef_630_metatranscriptome_file, cef_mock_metatranscriptome_file)
clinda_630_metatranscriptome <- read.delim(clinda_630_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(clinda_630_metatranscriptome) <- c('clinda_630_metaT_reads', 'ko', 'gene', 'pathway')
clinda_mock_metatranscriptome <- read.delim(clinda_mock_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(clinda_mock_metatranscriptome) <- c('clinda_mock_metaT_reads', 'ko', 'gene', 'pathway')
rm(clinda_630_metatranscriptome_file, clinda_mock_metatranscriptome_file)
strep_630_metatranscriptome <- read.delim(strep_630_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(strep_630_metatranscriptome) <- c('strep_630_metaT_reads', 'ko', 'gene', 'pathway')
strep_mock_metatranscriptome <- read.delim(strep_mock_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(strep_mock_metatranscriptome) <- c('strep_mock_metaT_reads', 'ko', 'gene', 'pathway')
rm(strep_630_metatranscriptome_file, strep_mock_metatranscriptome_file)
conv_metatranscriptome <- read.delim(conv_metatranscriptome_file, sep='\t', header=FALSE, row.names=1)
colnames(conv_metatranscriptome) <- c('conv_metaT_reads', 'ko', 'gene', 'pathway')
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
all_metagenome <- clean_merge(cef_metagenome, clinda_metagenome)
all_metagenome <- clean_merge(all_metagenome, strep_metagenome)
all_metagenome <- clean_merge(all_metagenome, conv_metagenome)
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
all_metatranscriptome <- clean_merge(cef_630_metatranscriptome, cef_mock_metatranscriptome)
all_metatranscriptome <- clean_merge(all_metatranscriptome, clinda_630_metatranscriptome)
all_metatranscriptome <- clean_merge(all_metatranscriptome, clinda_mock_metatranscriptome)
all_metatranscriptome <- clean_merge(all_metatranscriptome, strep_630_metatranscriptome)
all_metatranscriptome <- clean_merge(all_metatranscriptome, strep_mock_metatranscriptome)
all_metatranscriptome <- clean_merge(all_metatranscriptome, conv_metatranscriptome)
colnames(all_metatranscriptome) <- c('cefoperazone_630_metaT', 'cefoperazone_mock_metaT', 'clindamycin_630_metaT', 'clindamycin_mock_metaT', 
                                     'streptomycin_630_metaT', 'streptomycin_mock_metaT', 'conventional_metaT', 'ko', 'gene', 'pathway')
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
all_metagenome$cefoperazone <- t(rrarefy(all_metagenome$cefoperazone, sample=optimal_size))
all_metagenome$clindamycin <- t(rrarefy(all_metagenome$clindamycin, sample=optimal_size))
all_metagenome$streptomycin <- t(rrarefy(all_metagenome$streptomycin, sample=optimal_size))
all_metagenome$conventional <- t(rrarefy(all_metagenome$conventional, sample=optimal_size))

# Rarefy metatranscriptomic data
all_metatranscriptome$cefoperazone_630_metaT <- t(rrarefy(all_metatranscriptome$cefoperazone_630_metaT, sample=optimal_size))
all_metatranscriptome$cefoperazone_mock_metaT <- t(rrarefy(all_metatranscriptome$cefoperazone_mock_metaT, sample=optimal_size))
all_metatranscriptome$clindamycin_630_metaT <- t(rrarefy(all_metatranscriptome$clindamycin_630_metaT, sample=optimal_size))
all_metatranscriptome$clindamycin_mock_metaT <- t(rrarefy(all_metatranscriptome$clindamycin_mock_metaT, sample=optimal_size))
all_metatranscriptome$streptomycin_630_metaT <- t(rrarefy(all_metatranscriptome$streptomycin_630_metaT, sample=optimal_size))
all_metatranscriptome$streptomycin_mock_metaT <- t(rrarefy(all_metatranscriptome$streptomycin_mock_metaT, sample=optimal_size))
all_metatranscriptome$conventional_metaT <- t(rrarefy(all_metatranscriptome$conventional_metaT, sample=optimal_size))
rm(optimal_size)

# Merge metagenomes and metatranscriptomes
all_metagenome$ko <- NULL
all_metagenome$gene <- NULL
all_metagenome$pathway <- NULL
full_mapping <- merge(all_metagenome, all_metatranscriptome, by='row.names')
rownames(full_mapping) <- full_mapping$Row.names
full_mapping$Row.names <- NULL
colnames(full_mapping) <- c('cefoperazone_metaG', 'clindamycin_metaG', 'streptomycin_metaG', 'conventional_metaG', 'cefoperazone_630_metaT', 'cefoperazone_mock_metaT', 'clindamycin_630_metaT', 'clindamycin_mock_metaT', 'streptomycin_630_metaT', 'streptomycin_mock_metaT','conventional_metaT','ko', 'gene', 'pathway')
rm(all_metagenome, all_metatranscriptome)

# Normalize to metagenomic coverage
full_mapping$cefoperazone_630_metaT <- full_mapping$cefoperazone_630_metaT / full_mapping$cefoperazone_metaG
full_mapping$cefoperazone_mock_metaT <- full_mapping$cefoperazone_mock_metaT / full_mapping$cefoperazone_metaG
full_mapping$clindamycin_630_metaT <- full_mapping$clindamycin_630_metaT / full_mapping$clindamycin_metaG
full_mapping$clindamycin_mock_metaT <- full_mapping$clindamycin_mock_metaT / full_mapping$clindamycin_metaG
full_mapping$streptomycin_630_metaT <- full_mapping$streptomycin_630_metaT / full_mapping$streptomycin_metaG
full_mapping$streptomycin_mock_metaT <- full_mapping$streptomycin_mock_metaT / full_mapping$streptomycin_metaG
full_mapping$conventional_metaT <- full_mapping$conventional_metaT / full_mapping$conventional_metaG

# Separate treatment groups
cef_metatranscriptome <- full_mapping[,c(5,6,12:14)]
clinda_metatranscriptome <- full_mapping[,c(7,8,12:14)]
strep_metatranscriptome <- full_mapping[,c(9,10,12:14)]
conv_metatranscriptome <- full_mapping[,c(11,12:14)]
rm(full_mapping)

# Screen for active transcription in either condition
cef_metatranscriptome <- subset(cef_metatranscriptome, cef_metatranscriptome$cefoperazone_630 >= 1.0 | cef_metatranscriptome$cefoperazone_mock >= 1.0)
clinda_metatranscriptome <- subset(clinda_metatranscriptome, clinda_metatranscriptome$clindamycin_630 >= 1.0 | clinda_metatranscriptome$clindamycin_mock >= 1.0)
strep_metatranscriptome <- subset(strep_metatranscriptome, strep_metatranscriptome$streptomycin_630 >= 1.0 | strep_metatranscriptome$streptomycin_mock >= 1.0)
conv_metatranscriptome <- subset(conv_metatranscriptome, conv_metatranscriptome$conventional_metaT >= 1.0)

# Log2 transform the data
cef_metatranscriptome[,c(1,2)] <- log2(cef_metatranscriptome[,c(1,2)] + 1)
clinda_metatranscriptome[,c(1,2)] <- log2(clinda_metatranscriptome[,c(1,2)] + 1)
strep_metatranscriptome[,c(1,2)] <- log2(strep_metatranscriptome[,c(1,2)] + 1)
conv_metatranscriptome[,1] <- log2(conv_metatranscriptome[,1] + 1)

#-------------------------------------------------------------------------------------------------------------------------#

# Determine pathway annotations
cef_amino_sugar <- cef_metatranscriptome[grep('Amino_sugar', cef_metatranscriptome$pathway), ][,c(1,2)]
cef_fructose_mannose <- cef_metatranscriptome[grep('Fructose', cef_metatranscriptome$pathway), ][,c(1,2)]
cef_proline <- cef_metatranscriptome[grep('proline', cef_metatranscriptome$pathway), ][,c(1,2)]
cef_glycine <- cef_metatranscriptome[grep('Glycine', cef_metatranscriptome$pathway), ][,c(1,2)]
cef_galactose <- cef_metatranscriptome[grep('Galactose', cef_metatranscriptome$pathway), ][,c(1,2)]
cef_sucrose <- cef_metatranscriptome[grep('sucrose', cef_metatranscriptome$pathway), ][,c(1,2)]
cef_glycolysis <- cef_metatranscriptome[grep('Glycolysis', cef_metatranscriptome$pathway), ][,c(1,2)]

clinda_amino_sugar <- clinda_metatranscriptome[grep('Amino_sugar', clinda_metatranscriptome$pathway), ][,c(1,2)]
clinda_fructose_mannose <- clinda_metatranscriptome[grep('Fructose', clinda_metatranscriptome$pathway), ][,c(1,2)]
clinda_proline <- clinda_metatranscriptome[grep('proline', clinda_metatranscriptome$pathway), ][,c(1,2)]
clinda_glycine <- clinda_metatranscriptome[grep('Glycine', clinda_metatranscriptome$pathway), ][,c(1,2)]
clinda_galactose <- clinda_metatranscriptome[grep('Galactose', clinda_metatranscriptome$pathway), ][,c(1,2)]
clinda_sucrose <- clinda_metatranscriptome[grep('sucrose', clinda_metatranscriptome$pathway), ][,c(1,2)]
clinda_glycolysis <- clinda_metatranscriptome[grep('Glycolysis', clinda_metatranscriptome$pathway), ][,c(1,2)]

strep_amino_sugar <- strep_metatranscriptome[grep('Amino_sugar', strep_metatranscriptome$pathway), ][,c(1,2)]
strep_fructose_mannose <- strep_metatranscriptome[grep('Fructose', strep_metatranscriptome$pathway), ][,c(1,2)]
strep_proline <- strep_metatranscriptome[grep('proline', strep_metatranscriptome$pathway), ][,c(1,2)]
strep_glycine <- strep_metatranscriptome[grep('Glycine', strep_metatranscriptome$pathway), ][,c(1,2)]
strep_galactose <- strep_metatranscriptome[grep('Galactose', strep_metatranscriptome$pathway), ][,c(1,2)]
strep_sucrose <- strep_metatranscriptome[grep('sucrose', strep_metatranscriptome$pathway), ][,c(1,2)]
strep_glycolysis <- strep_metatranscriptome[grep('Glycolysis', strep_metatranscriptome$pathway), ][,c(1,2)]

conv_amino_sugar <- conv_metatranscriptome[grep('Amino_sugar', conv_metatranscriptome$pathway), ][,1]
conv_fructose_mannose <- conv_metatranscriptome[grep('Fructose', conv_metatranscriptome$pathway), ][,1]
conv_proline <- conv_metatranscriptome[grep('proline', conv_metatranscriptome$pathway), ][,1]
conv_glycine <- conv_metatranscriptome[grep('Glycine', conv_metatranscriptome$pathway), ][,1]
conv_galactose <- conv_metatranscriptome[grep('Galactose', conv_metatranscriptome$pathway), ][,1]
conv_sucrose <- conv_metatranscriptome[grep('sucrose', conv_metatranscriptome$pathway), ][,1]
conv_glycolysis <- conv_metatranscriptome[grep('Glycolysis', conv_metatranscriptome$pathway), ][,1]

#-------------------------------------------------------------------------------------------------------------------------#

# Calculate the distance of all points from x=y
# Will reveal which genes were most effected by c. diff colonization




#good.dist <- sqrt((good.ord - typ.ord)^2 / 2)



#-------------------------------------------------------------------------------------------------------------------------#

# Define which pathway to plot and the ouput file name
plot_file <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/results/figures/figure_3.pdf'
pdf(file=plot_file, width=10, height=10)
layout(matrix(c(1,1,2,2,
                3,4,5,6,
                7,7,8,8), 
              nrow=3, ncol=4, byrow = TRUE))









#-------------------#

# Streptomycin
par(mar=c(4, 4, 1, 1), mgp=c(2.4,0.7,0))
plot(x=strep_metatranscriptome$streptomycin_630_metaT, y=strep_metatranscriptome$streptomycin_mock_metaT, 
     xlim=c(0,4), ylim=c(0,4), pch=20, col='gray35', xaxt='n', yaxt='n', 
     xlab='Normalized cDNA Read Abundance (Infected)', ylab='Normalized cDNA Read Abundance (Mock)')
legend('topleft', 'Streptomycin-treated', bty='n', cex=3)
segments(-2, -2, 8, 8, lwd=2, lty=2)
axis(side=1, at=c(0:4), axis_labels, tick=TRUE)
axis(side=2, at=c(0:4), axis_labels, tick=TRUE, las=1)

points(strep_glycolysis, cex=3.5, pch=21, bg='chartreuse3', col='black')
points(strep_amino_sugar, cex=3.5, pch=21, bg='darkcyan', col='black')
points(strep_fructose_mannose, cex=3.5, pch=21, bg='firebrick3', col='black')
points(strep_proline, cex=3.5, pch=21, bg='darkgoldenrod1', col='black')
points(strep_glycine, cex=3.5, pch=21, bg='darkgoldenrod1', col='black')
points(strep_galactose, cex=3.5, pch=21, bg='darkorchid3', col='black')

#-------------------#

# Cefoperazone
par(mar=c(4, 4, 1, 1), mgp=c(2.4,0.7,0), xaxs='r', yaxs='r')
plot(x=cef_metatranscriptome$cefoperazone_630_metaT, y=cef_metatranscriptome$cefoperazone_mock_metaT, 
     xlim=c(0,4), ylim=c(0,4), pch=20, cex=1.5, col='gray35', xaxt='n', yaxt='n', 
     xlab='Normalized cDNA Read Abundance (Infected)', ylab='Normalized cDNA Read Abundance (Mock)')
segments(-2, -2, 8, 8, lwd=2, lty=2)
axis_labels <- parse(text=paste(rep(10,4), '^', seq(0,4,1), sep=''))
axis(side=1, at=c(0:4), axis_labels, tick=TRUE)
axis(side=2, at=c(0:4), axis_labels, tick=TRUE, las=1)
legend('topleft', 'Cefoperazone-treated', bty='n', cex=3)

points(cef_glycolysis, cex=3.5, pch=21, bg='chartreuse3', col='black')
points(cef_amino_sugar, cex=3.5, pch=21, bg='darkcyan', col='black')
points(cef_fructose_mannose, cex=3.5, pch=21, bg='firebrick3', col='black')
points(cef_proline, cex=3.5, pch=21, bg='darkgoldenrod1', col='black')
points(cef_glycine, cex=3.5, pch=21, bg='darkgoldenrod1', col='black')
points(cef_galactose, cex=3.5, pch=21, bg='darkorchid3', col='black')

#-------------------#

# Clindamycin
par(mar=c(4, 4, 1, 1), mgp=c(2.4,0.7,0), xaxs='r', yaxs='r')
plot(x=clinda_metatranscriptome$clindamycin_630_metaT, y=clinda_metatranscriptome$clindamycin_mock_metaT, 
     xlim=c(0,4), ylim=c(0,4), pch=20, col='gray35', xaxt='n', yaxt='n', 
     xlab='Normalized cDNA Read Abundance (Infected)', ylab='Normalized cDNA Read Abundance (Mock)')
legend('topleft', 'Clindamycin-treated', bty='n', cex=3)
segments(-2, -2, 8, 8, lwd=2, lty=2)
axis(side=1, at=c(0:4), axis_labels, tick=TRUE)
axis(side=2, at=c(0:4), axis_labels, tick=TRUE, las=1)

points(clinda_glycolysis, cex=3.5, pch=21, bg='chartreuse3', col='black')
points(clinda_amino_sugar, cex=3.5, pch=21, bg='darkcyan', col='black')
points(clinda_fructose_mannose, cex=3.5, pch=21, bg='firebrick3', col='black')
points(clinda_proline, cex=3.5, pch=21, bg='darkgoldenrod1', col='black')
points(clinda_glycine, cex=3.5, pch=21, bg='darkgoldenrod1', col='black')
points(clinda_galactose, cex=3.5, pch=21, bg='darkorchid3', col='black')

#-------------------#

plot(1, type="n", axes=F, xlab="", ylab="")
legend('center', legend=c('Glycolysis', 'Amino sugar metabolism', 'Fructose/Mannose metabolism', 'Proline/Glycine metabolism', 'Galactose metabolism'), 
       pt.bg=c('chartreuse3', 'darkcyan', 'firebrick3', 'darkgoldenrod1', 'darkorchid3'), 
       pch=21, cex=1.6, pt.cex=3, ncol=5)

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
dev.off()
#rm(cef_metatranscriptome, clinda_metatranscriptome, strep_metatranscriptome, plot_file, axis_labels)
