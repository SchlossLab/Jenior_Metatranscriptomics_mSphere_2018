
# Set up environment
rm(list=ls())
gc()

# Load in functions
starting_dir <- getwd()
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/code/R/functions.R')

# Output plot name
plot_file <- 'results/supplement/figures/figure_S1.pdf'

# Input 0.03 OTU shared file
shared_otu_file <- 'data/16S_analysis/all_treatments.0.03.unique_list.conventional.shared'

# Input Metadata
metadata_file <- 'data/metadata.tsv'

#----------------#

# Read in data

# 16S data
shared_otu <- read.delim(shared_otu_file, sep='\t', header=TRUE, row.names=2)
rm(shared_otu_file)

# Metadata
metadata <- read.delim(metadata_file, sep='\t', header=T, row.names=1)
rm(metadata_file)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data

# Metadata
metadata <- metadata[!rownames(metadata) %in% c('CefC5M2'), ] # Remove contaminated sample
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL

# 16S data
shared_otu$label <- NULL
shared_otu$numOtus <- NULL
shared_otu <- shared_otu[!rownames(shared_otu) %in% c('CefC5M2'), ] # Contaminated sample
shared_otu <- shared_otu[,!(names(shared_otu) %in% c('Otu0004','Otu0308'))] # Remove residual C. difficile OTUs
#shared_otu <- shared_otu[!rownames(shared_otu) %in% c('ConvC1M1','ConvC1M2','ConvC1M3','ConvC1M4',
#                                                       'ConvC2M1','ConvC2M2','ConvC2M3','ConvC2M4','ConvC2M5'), ] # Untreated SPF samples
sub_sample <- ceiling(min(rowSums(shared_otu)) * 0.9) # calculate rarefaction level
shared_otu <- rrarefy(shared_otu, sub_sample) # subsample shared file
shared_otu <- filter_table(shared_otu)
rm(sub_sample)

#-------------------------------------------------------------------------------------------------------------------------#

# Calculate axes
# Community structure - all
otu_nmds <- metaMDS(shared_otu, k=2, trymax=100, distance='bray')$points
otu_nmds <- clean_merge(metadata, otu_nmds)

# Subset NMDS axes to color points - all
otu_cefoperazone <- subset(otu_nmds, abx == 'cefoperazone')
otu_cefoperazone_630 <- subset(otu_cefoperazone, infection == '630')
otu_cefoperazone_mock <- subset(otu_cefoperazone, infection == 'mock')
rm(otu_cefoperazone)
otu_clindamycin <- subset(otu_nmds, abx == 'clindamycin')
otu_clindamycin_630 <- subset(otu_clindamycin, infection == '630')
otu_clindamycin_mock <- subset(otu_clindamycin, infection == 'mock')
rm(otu_clindamycin)
otu_streptomycin <- subset(otu_nmds, abx == 'streptomycin')
otu_streptomycin_630 <- subset(otu_streptomycin, infection == '630')
otu_streptomycin_mock <- subset(otu_streptomycin, infection == 'mock')
rm(otu_streptomycin)
otu_noantibiotics <- subset(otu_nmds, abx == 'none')

#----------------#

# Separate plots
# All abx
abx_otu <- shared_otu[!rownames(shared_otu) %in% c('ConvC1M1','ConvC1M2','ConvC1M3','ConvC1M4',
                                                   'ConvC2M1','ConvC2M2','ConvC2M3','ConvC2M4','ConvC2M5'), ]
# Individual abx
cef_otu <- shared_otu[rownames(shared_otu) %in% c('CefC1M1','CefC1M2','CefC1M3',
                                                  'CefC2M1','CefC2M2','CefC2M3',
                                                  'CefC3M1','CefC3M2','CefC3M3',
                                                  'CefC4M1','CefC4M2','CefC4M3',
                                                  'CefC5M1','CefC5M2','CefC5M3',
                                                  'CefC6M1','CefC6M2','CefC6M3'), ]
clinda_otu <- shared_otu[rownames(shared_otu) %in% c('ClindaC1M1','ClindaC1M2','ClindaC1M3',
                                                     'ClindaC2M1','ClindaC2M2','ClindaC2M3',
                                                     'ClindaC3M1','ClindaC3M2','ClindaC3M3',
                                                     'ClindaC4M1','ClindaC4M2','ClindaC4M3',
                                                     'ClindaC5M1','ClindaC5M2','ClindaC5M3',
                                                     'ClindaC6M1','ClindaC6M2','ClindaC6M3'), ]
strep_otu <- shared_otu[rownames(shared_otu) %in% c('StrepC1M1','StrepC1M2','StrepC1M3',
                                                    'StrepC2M1','StrepC2M2','StrepC2M3',
                                                    'StrepC3M1','StrepC3M2','StrepC3M3',
                                                    'StrepC4M1','StrepC4M2','StrepC4M3',
                                                    'StrepC5M1','StrepC5M2','StrepC5M3',
                                                    'StrepC6M1','StrepC6M2','StrepC6M3'), ]

# Calculate inv-Simpson diversity
cef_infected_diversity <- diversity(cef_otu[rownames(cef_otu) %in% c('CefC1M1','CefC1M2','CefC1M3',
                                                    'CefC2M1','CefC2M2','CefC2M3',
                                                    'CefC3M1','CefC3M2','CefC3M3'), ], 'invsimpson')
cef_mock_diversity <- diversity(cef_otu[rownames(cef_otu) %in% c('CefC4M1','CefC4M2','CefC4M3',
                                                'CefC5M1','CefC5M2','CefC5M3',
                                                'CefC6M1','CefC6M2','CefC6M3'), ], 'invsimpson')
clinda_infected_diversity <- diversity(clinda_otu[rownames(clinda_otu) %in% c('ClindaC1M1','ClindaC1M2','ClindaC1M3',
                                                          'ClindaC2M1','ClindaC2M2','ClindaC2M3',
                                                          'ClindaC3M1','ClindaC3M2','ClindaC3M3'), ], 'invsimpson')
clinda_mock_diversity <- diversity(clinda_otu[rownames(clinda_otu) %in% c('ClindaC4M1','ClindaC4M2','ClindaC4M3',
                                                      'ClindaC5M1','ClindaC5M2','ClindaC5M3',
                                                      'ClindaC6M1','ClindaC6M2','ClindaC6M3'), ], 'invsimpson')
strep_infected_diversity <- diversity(strep_otu[rownames(strep_otu) %in% c('StrepC1M1','StrepC1M2','StrepC1M3',
                                                        'StrepC2M1','StrepC2M2','StrepC2M3',
                                                        'StrepC3M1','StrepC3M2','StrepC3M3'), ], 'invsimpson')
strep_mock_diversity <- diversity(strep_otu[rownames(strep_otu) %in% c('StrepC4M1','StrepC4M2','StrepC4M3',
                                                    'StrepC5M1','StrepC5M2','StrepC5M3',
                                                    'StrepC6M1','StrepC6M2','StrepC6M3'), ], 'invsimpson')
noabx_diversity <- diversity(shared_otu[rownames(shared_otu) %in% c('ConvC1M1','ConvC1M2','ConvC1M3','ConvC1M4',
                                                    'ConvC2M1','ConvC2M2','ConvC2M3','ConvC2M4','ConvC2M5'), ], 'invsimpson')

# Test differences in diversity
# Untreated vs abx
p.adjust(c(wilcox.test(noabx_diversity, cef_infected_diversity, exact=F)$p.value,
           wilcox.test(noabx_diversity, cef_mock_diversity, exact=F)$p.value,
           wilcox.test(noabx_diversity, clinda_infected_diversity, exact=F)$p.value,
           wilcox.test(noabx_diversity, clinda_mock_diversity, exact=F)$p.value,
           wilcox.test(noabx_diversity, strep_infected_diversity, exact=F)$p.value,
           wilcox.test(noabx_diversity, strep_mock_diversity, exact=F)$p.value), method='BH')
# Strep
wilcox.test(strep_infected_diversity, strep_mock_diversity, exact=F)$p.value
# Cef
wilcox.test(cef_infected_diversity, cef_mock_diversity, exact=F)$p.value
# Clinda
wilcox.test(clinda_infected_diversity, clinda_mock_diversity, exact=F)$p.value


# Calculate axes and merge with metadata
abx_otu_nmds <- metaMDS(abx_otu, k=2, trymax=100)$points
abx_otu_nmds <- clean_merge(metadata, abx_otu_nmds)
cef_otu_nmds <- metaMDS(cef_otu, k=2, trymax=100)$points
cef_otu_nmds <- clean_merge(metadata, cef_otu_nmds)
clinda_otu_nmds <- metaMDS(clinda_otu, k=2, trymax=100)$points
clinda_otu_nmds <- clean_merge(metadata, clinda_otu_nmds)
strep_otu_nmds <- metaMDS(strep_otu, k=2, trymax=100)$points
strep_otu_nmds <- clean_merge(metadata, strep_otu_nmds)
rm(metadata)

# Calculate significant differences
otu_p <- as.character(anosim(shared_otu, otu_nmds$susceptibility, permutations=999, distance='bray')$signif)
otu_r <- as.character(round(anosim(shared_otu, otu_nmds$susceptibility, permutations=999, distance='bray')$statistic, 3))
abx_otu_p <- as.character(anosim(abx_otu, abx_otu_nmds$abx, permutations=999, distance='bray')$signif)
abx_otu_r <- as.character(round(anosim(abx_otu, abx_otu_nmds$abx, permutations=999, distance='bray')$statistic, 3))
cef_otu_p <- as.character(anosim(cef_otu, cef_otu_nmds$infection, permutations=999, distance='bray')$signif)
cef_otu_r <- as.character(round(anosim(cef_otu, cef_otu_nmds$infection, permutations=999, distance='bray')$statistic, 3))
clinda_otu_p <- as.character(anosim(clinda_otu, clinda_otu_nmds$infection, permutations=999, distance='bray')$signif)
clinda_otu_r <- as.character(round(anosim(clinda_otu, clinda_otu_nmds$infection, permutations=999, distance='bray')$statistic, 3))
strep_otu_p <- as.character(anosim(strep_otu, strep_otu_nmds$infection, permutations=999, distance='bray')$signif)
strep_otu_r <- as.character(round(anosim(strep_otu, strep_otu_nmds$infection, permutations=999, distance='bray')$statistic, 3))
rm(cef_otu, clinda_otu, strep_otu)

# Move points for easier viewing
strep_otu_nmds$MDS2 <- strep_otu_nmds$MDS2 + 0.2

# Subset to points for plot
otu_abx_cef <- subset(abx_otu_nmds, abx == 'cefoperazone')
otu_abx_cef_630 <- subset(otu_abx_cef, infection == '630')
otu_abx_cef_mock <- subset(otu_abx_cef, infection == 'mock')
rm(otu_abx_cef)
otu_abx_clinda <- subset(abx_otu_nmds, abx == 'clindamycin')
otu_abx_clinda_630 <- subset(otu_abx_clinda, infection == '630')
otu_abx_clinda_mock <- subset(otu_abx_clinda, infection == 'mock')
rm(otu_abx_clinda)
otu_abx_strep <- subset(abx_otu_nmds, abx == 'streptomycin')
otu_abx_strep_630 <- subset(otu_abx_strep, infection == '630')
otu_abx_strep_mock <- subset(otu_abx_strep, infection == 'mock')
rm(otu_abx_strep)
cef_otu_nmds_630 <- subset(cef_otu_nmds, infection == '630')
cef_otu_nmds_mock <- subset(cef_otu_nmds, infection == 'mock')
clinda_otu_nmds_630 <- subset(clinda_otu_nmds, infection == '630')
clinda_otu_nmds_mock <- subset(clinda_otu_nmds, infection == 'mock')
strep_otu_nmds_630 <- subset(strep_otu_nmds, infection == '630')
strep_otu_nmds_mock <- subset(strep_otu_nmds, infection == 'mock')

# Calculate centroids
cef_otu_centoids <- aggregate(cbind(cef_otu_nmds$MDS1,cef_otu_nmds$MDS2)~cef_otu_nmds$infection, data=cef_otu_nmds, mean)
clinda_otu_centoids <- aggregate(cbind(clinda_otu_nmds$MDS1,clinda_otu_nmds$MDS2)~clinda_otu_nmds$infection, data=clinda_otu_nmds, mean)
strep_otu_centoids <- aggregate(cbind(strep_otu_nmds$MDS1,strep_otu_nmds$MDS2)~strep_otu_nmds$infection, data=strep_otu_nmds, mean)

#-------------------------------------------------------------------------------------------------------------------------#

# Plot the figure
pdf(file=plot_file, width=12, height=12)
layout(matrix(c(1,1,1,
                2,3,4,
                5,6,7),
              nrow=3, ncol=3, byrow=TRUE))

#-------------------#

# Diversity
par(mar=c(3,4,1,1), mgp=c(2.6,1,0), las=1, yaxs='i', xaxs='r')
# No abx
stripchart(at=0.5, noabx_diversity, xlim=c(0,10), ylim=c(0,20), vertical=TRUE, cex.axis=1.4, cex.lab=1.5,
           pch=21, bg=noabx_col, method='jitter', jitter=0.15, cex=2.5, lwd=1.2, ylab='Inv. Simpson Diversity')
segments(x0=0.2, y0=median(noabx_diversity), x1=0.8, y1=median(noabx_diversity), lwd=2.5)
#Strep
stripchart(at=2.5, strep_infected_diversity, vertical=TRUE, 
           pch=21, bg=strep_col, method='jitter', jitter=0.15, cex=2.5, lwd=1.2, add=TRUE)
segments(x0=2.2, y0=median(strep_infected_diversity), x1=2.8, y1=median(strep_infected_diversity), lwd=2.5)
stripchart(at=3.5, strep_mock_diversity, vertical=TRUE, 
           pch=21, bg=strep_col, method='jitter', jitter=0.15, cex=2.5, lwd=1.2, add=TRUE)
segments(x0=3.2, y0=median(strep_mock_diversity), x1=3.8, y1=median(strep_mock_diversity), lwd=2.5)
# Cef
stripchart(at=5.5, cef_infected_diversity, vertical=TRUE, 
           pch=21, bg=cef_col, method='jitter', jitter=0.15, cex=2.5, lwd=1.2, add=TRUE)
segments(x0=5.2, y0=median(cef_infected_diversity), x1=5.8, y1=median(cef_infected_diversity), lwd=2.5)
stripchart(at=6.5, cef_mock_diversity, vertical=TRUE, 
           pch=21, bg=cef_col, method='jitter', jitter=0.15, cex=2.5, lwd=1.2, add=TRUE)
segments(x0=6.2, y0=median(cef_mock_diversity), x1=6.8, y1=median(cef_mock_diversity), lwd=2.5)
# Clinda
stripchart(at=8.5, clinda_infected_diversity, vertical=TRUE, 
           pch=21, bg=clinda_col, method='jitter', jitter=0.15, cex=2.5, lwd=1.2, add=TRUE)
segments(x0=8.2, y0=median(clinda_infected_diversity), x1=8.8, y1=median(clinda_infected_diversity), lwd=2.5)
stripchart(at=9.5, clinda_mock_diversity, vertical=TRUE, 
           pch=21, bg=clinda_col, method='jitter', jitter=0.15, cex=2.5, lwd=1.2, add=TRUE)
segments(x0=9.2, y0=median(clinda_mock_diversity), x1=9.8, y1=median(clinda_mock_diversity), lwd=2.5)
mtext('CDI:', side=1, at=-0.5, padj=1, cex=1.1)
mtext(c('-','+','-','+','-','+','-'), side=1, 
      at=c(0.5, 2.5,3.5, 5.5,6.5, 8.5,9.5), padj=1, cex=1.3)
abline(v=c(1.5,4.5,7.5), lty=5)
mtext('a', side=2, line=2, las=2, adj=2, padj=-10, cex=1.4, font=2)
text(x=c(0.5, 2.5,3.5, 5.5,6.5, 8.5,9.5), y=19.5, labels='*', col=noabx_col, cex=2.5, font=2) # Significant difference fromm untreated
segments(x0=c(5.5,8.5), y0=5, x1=c(6.5,9.5), y1=5, lwd=2.5)
text(x=c(6,9), y=6, labels='*', cex=2.5, font=2) # Significant difference within groups

#-------------------#

# All groups - NMDS
par(mar=c(4,4,1,1), las=1, mgp=c(2.5,0.75,0), xaxs='i', yaxs='i')
plot(x=otu_nmds$MDS1, y=otu_nmds$MDS2, xlim=c(-1.5,2), ylim=c(-1.5,1.5),
     xlab='NMDS axis 1', ylab='NMDS axis 2', xaxt='n', yaxt='n', pch=19, cex=0.2)
axis(side=1, at=seq(-1.5,2.0,0.5), labels=seq(-1.5,2.0,0.5))
axis(side=2, at=seq(-2.0,1.5,0.5), labels=seq(-2.0,1.5,0.5))
mtext('b', side=2, line=2, las=2, adj=1.4, padj=-9, cex=1.4, font=2)
legend('topleft', legend='All groups', pch=1, cex=1.4, pt.cex=0, bty='n')
points(x=otu_cefoperazone_630$MDS1, y=otu_cefoperazone_630$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=otu_clindamycin_630$MDS1, y=otu_clindamycin_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=otu_streptomycin_630$MDS1, y=otu_streptomycin_630$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=otu_cefoperazone_mock$MDS1, y=otu_cefoperazone_mock$MDS2, bg=cef_col, pch=24, cex=1.8, lwd=1.2)
points(x=otu_clindamycin_mock$MDS1, y=otu_clindamycin_mock$MDS2, bg=clinda_col, pch=24, cex=1.8, lwd=1.2)
points(x=otu_streptomycin_mock$MDS1, y=otu_streptomycin_mock$MDS2, bg=strep_col, pch=24, cex=1.8, lwd=1.2)
points(x=otu_noantibiotics$MDS1, y=otu_noantibiotics$MDS2, bg=noabx_col, pch=24, cex=1.8, lwd=1.2)
legend('bottomright', legend=c('Resistant vs Susceptible:',
                               as.expression(bquote(paste(italic('R'),' = ',.(otu_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(otu_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')

#-------------------#

# Abx only
plot(x=abx_otu_nmds$MDS1-0.25, y=abx_otu_nmds$MDS2, xlim=c(-1.3,1.3), ylim=c(-1.8,1.8),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('c', side=2, line=2, las=2, adj=1.4, padj=-9, cex=1.4, font=2)
legend('topleft', legend='Antibiotic pretreatments only', pch=1, cex=1.4, pt.cex=0, bty='n')
points(x=otu_abx_cef_630$MDS1-0.25, y=otu_abx_cef_630$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=otu_abx_cef_mock$MDS1-0.25, y=otu_abx_cef_mock$MDS2, bg=cef_col, pch=24, cex=2, lwd=1.2)
points(x=otu_abx_clinda_630$MDS1-0.25, y=otu_abx_clinda_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=otu_abx_clinda_mock$MDS1-0.25, y=otu_abx_clinda_mock$MDS2, bg=clinda_col, pch=24, cex=2, lwd=1.2)
points(x=otu_abx_strep_630$MDS1-0.25, y=otu_abx_strep_630$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=otu_abx_strep_mock$MDS1-0.25, y=otu_abx_strep_mock$MDS2, bg=strep_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c('Strep vs Cef vs Clinda:',
                               as.expression(bquote(paste(italic('R'),' = ',.(abx_otu_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(abx_otu_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')

#-------------------#

# Streptomycin
plot(x=strep_otu_nmds$MDS1, y=strep_otu_nmds$MDS2, xlim=c(-1.3,1.3), ylim=c(-1.8,1.8),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('d', side=2, line=2, las=2, adj=1.4, padj=-9, cex=1.4, font=2)
legend('topleft', legend='Streptomycin-pretreated', pch=1, cex=1.4, pt.cex=0, bty='n')
segments(x0=strep_otu_nmds_630$MDS1, y0=strep_otu_nmds_630$MDS2, x1=strep_otu_centoids[1,2], y1=strep_otu_centoids[1,3], col='gray30')
segments(x0=strep_otu_nmds_mock$MDS1, y0=strep_otu_nmds_mock$MDS2, x1=strep_otu_centoids[2,2], y1=strep_otu_centoids[2,3], col='gray30')
points(x=strep_otu_nmds_630$MDS1, y=strep_otu_nmds_630$MDS2, bg=strep_col, pch=21, cex=2, lwd=1.2)
points(x=strep_otu_nmds_mock$MDS1, y=strep_otu_nmds_mock$MDS2, bg=strep_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c('Mock vs Infected:',
                               as.expression(bquote(paste(italic('R'),' = ',.(strep_otu_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(strep_otu_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')

#-------------------#

# Cefoperazone
plot(x=cef_otu_nmds$MDS1, y=cef_otu_nmds$MDS2, xlim=c(-1.2,1.2), ylim=c(-0.8,0.8),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('e', side=2, line=2, las=2, adj=1.4, padj=-9, cex=1.4, font=2)
legend('topleft', legend='Cefoperazone-pretreated', pch=1, cex=1.4, pt.cex=0, bty='n')
segments(x0=cef_otu_nmds_630$MDS1, y0=cef_otu_nmds_630$MDS2, x1=cef_otu_centoids[1,2], y1=cef_otu_centoids[1,3], col='gray30')
segments(x0=cef_otu_nmds_mock$MDS1, y0=cef_otu_nmds_mock$MDS2, x1=cef_otu_centoids[2,2], y1=cef_otu_centoids[2,3], col='gray30')
points(x=cef_otu_nmds_630$MDS1, y=cef_otu_nmds_630$MDS2, bg=cef_col, pch=21, cex=2, lwd=1.2)
points(x=cef_otu_nmds_mock$MDS1, y=cef_otu_nmds_mock$MDS2, bg=cef_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c('Mock vs Infected:',
                               as.expression(bquote(paste(italic('R'),' = ',.(cef_otu_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(cef_otu_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')

#-------------------#

# Clindamycin
plot(x=clinda_otu_nmds$MDS1, y=clinda_otu_nmds$MDS2, xlim=c(-1.1,1.1), ylim=c(-1.2,1.2),
     xlab='NMDS axis 1', ylab='NMDS axis 2', pch=19, cex=0.2)
mtext('f', side=2, line=2, las=2, adj=1.4, padj=-9, cex=1.4, font=2)
legend('topleft', legend='Clindamycin-pretreated', pch=1, cex=1.4, pt.cex=0, bty='n')
segments(x0=clinda_otu_nmds_630$MDS1, y0=clinda_otu_nmds_630$MDS2, x1=clinda_otu_centoids[1,2], y1=clinda_otu_centoids[1,3], col='gray30')
segments(x0=clinda_otu_nmds_mock$MDS1, y0=clinda_otu_nmds_mock$MDS2, x1=clinda_otu_centoids[2,2], y1=clinda_otu_centoids[2,3], col='gray30')
points(x=clinda_otu_nmds_630$MDS1, y=clinda_otu_nmds_630$MDS2, bg=clinda_col, pch=21, cex=2, lwd=1.2)
points(x=clinda_otu_nmds_mock$MDS1, y=clinda_otu_nmds_mock$MDS2, bg=clinda_col, pch=24, cex=2, lwd=1.2)
legend('bottomright', legend=c('Mock vs Infected:',
                               as.expression(bquote(paste(italic('R'),' = ',.(clinda_otu_r)))),
                               as.expression(bquote(paste(italic('p'),' = ',.(clinda_otu_p))))), pch=1, cex=1.3, pt.cex=0, bty='n')

#-------------------#

# Legends
par(mar=c(1,1,1,1))
plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-5,5), ylim=c(-5,5))
legend(x=-3, y=3, legend=c('No Antibiotics','Streptomycin-pretreated','Cefoperzone-pretreated','Clindamycin-pretreated'), 
       pt.bg=c(noabx_col,strep_col,cef_col,clinda_col), pch=22, cex=1.5, pt.cex=2.7)
legend(x=-2.5, y=0, legend=c(as.expression(bquote(paste(italic('C. difficile'),'-infected'))),'Mock-infected'), 
       col='black', pch=c(16,17), cex=1.5, pt.cex=c(2.5,2.2))

dev.off()

#-------------------------------------------------------------------------------------------------------------------------#

# Clean up
for (dep in deps){
  pkg <- paste('package:', dep, sep='')
  detach(pkg, character.only = TRUE)
}
setwd(starting_dir)
rm(list=ls())
gc()
