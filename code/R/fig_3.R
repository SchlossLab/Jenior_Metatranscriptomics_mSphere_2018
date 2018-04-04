
# Set up environment
starting_dir <- getwd()
rm(list=ls())
gc()

# Load in functions
source('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/code/R/functions.R')

# Output plot name
fig3_plot <- 'results/figures/figure_3.pdf'

# Input Metabolomes
metabolome <- 'data/metabolome/scaled_intensities.log10.tsv'

# Input Metadata
metadata <- 'data/metadata.tsv'

#----------------#

# Read in data

# Metabolomes
metabolome <- read.delim(metabolome, sep='\t', header=TRUE)

# Metadata
metadata <- read.delim(metadata, sep='\t', header=T, row.names=1)

#-------------------------------------------------------------------------------------------------------------------------#

# Format data

# Metadata
metadata$type <- NULL
metadata$cage <- NULL
metadata$mouse <- NULL
metadata$gender <- NULL

# Metabolomes
metabolome$BIOCHEMICAL <- gsub('_', ' ', metabolome$BIOCHEMICAL)
rownames(metabolome) <- metabolome$BIOCHEMICAL
metabolome$BIOCHEMICAL <- NULL
metabolome$PUBCHEM <- NULL
metabolome$KEGG <- NULL
metabolome$SUB_PATHWAY <- NULL
metabolome$SUPER_PATHWAY <- NULL
metabolome <- as.data.frame(t(metabolome))
colnames(metabolome) <- make.names(colnames(metabolome))
metabolome <- data.frame(sapply(metabolome, function(x) as.numeric(as.character(x))))

#-------------------------------------------------------------------------------------------------------------------------#

# Stats
# Prep data
metabolome <- clean_merge(metadata, metabolome)
metabolome <- subset(metabolome, abx != 'germfree')
resistant_metabolome <- subset(metabolome, abx == 'none')
resistant_metabolome$susceptibility <- NULL
resistant_metabolome$abx <- NULL
resistant_metabolome$infection <- NULL
metabolome <- subset(metabolome, abx != 'none')
metabolome$susceptibility <- NULL

# Stickland metabolites
stickland_aa <- c("proline","trans.4.hydroxyproline","pro.hydroxy.pro",
                  'thioproline','N.acetylthreonine','N.acetylserine',
                  "glycine",
                  "alanine","N.acetylalanine",
                  "valine","N.acetylvaline",
                  "leucine","N.acetylleucine",
                  "isoleucine","N.acetylisoleucine",
                  "phenylalanine","N.acetylphenylalanine",
                  "N.acetyltyrosine",
                  "tyrosine","N.acetylglycine")
#stickland_keto <- c('infection',"alpha.hydroxyisocaproate","X3..4.hydroxyphenyl.lactate","phenyllactate..PLA.")
stickland_carboxy <- c("isocaproate","lactate","X5.aminovalerate","X2.3.dihydroxyisovalerate",
                       "X2..4.hydroxyphenyl.propionate","X3..3.hydroxyphenyl.propionate","X3..4.hydroxyphenyl.propionate",
                       "isovalerate","alpha.hydroxyisovalerate")

# Get metabolites in Stickland fermentation pathway
metabolome_aa <- metabolome[, c('infection',stickland_aa)]
metabolome_carboxy <- metabolome[, c('infection',stickland_carboxy)]
resistant_metabolome_aa <- resistant_metabolome[, stickland_aa]
resistant_metabolome_carboxy <- resistant_metabolome[, stickland_carboxy]
rm(metabolome, resistant_metabolome, stickland_aa, stickland_carboxy)

#------------#

# Reformat metabolite names
colnames(metabolome_aa) <- c('abx','infection',
                             "Proline","4-Hydroxyproline","Prolyl-Hydroxyproline",
                             'Thioproline','N-Acetylthreonine','N-Acetylserine',
                             "Glycine","Alanine","N-Acetylalanine","Valine","N-Acetylvaline",
                             "Leucine","N-Acetylleucine","Isoleucine","N-Acetylisoleucine",
                             "Phenylalanine","N-Acetylphenylalanine","N-Acetyltyrosine",
                             "Tyrosine","N-Acetylglycine")
colnames(resistant_metabolome_aa) <- c('abx','infection',
                                       "Proline","4-Hydroxyproline","Prolyl-Hydroxyproline",
                                       'Thioproline','N-Acetylthreonine','N-Acetylserine',
                                       "Glycine","Alanine","N-Acetylalanine","Valine","N-Acetylvaline",
                                       "Leucine","N-Acetylleucine","Isoleucine","N-Acetylisoleucine",
                                       "Phenylalanine","N-Acetylphenylalanine","N-Acetyltyrosine",
                                       "Tyrosine","N-Acetylglycine")
colnames(metabolome_carboxy) <- c('abx','infection',
                                  "Isocaproate","Lactate","5-Aminovalerate",
                                  "2-3-Dihydroxy-isovalerate",
                                  "2-4-Hydroxyphenyl-propionate","X3..3.hydroxyphenyl.propionate",
                                  "3-4-Hydroxyphenyl-propionate",
                                  "Isovalerate","alpha-Hydroxyisovalerate")
colnames(resistant_metabolome_carboxy) <- c('abx','infection',
                                            "Isocaproate","Lactate","5-Aminovalerate",
                                            "2-3-Dihydroxy-isovalerate",
                                            "2-4-Hydroxyphenyl-propionate","X3..3.hydroxyphenyl.propionate",
                                            "3-4-Hydroxyphenyl-propionate",
                                            "Isovalerate","alpha-Hydroxyisovalerate")

#------------#

# Subset groups to mock and infected
metabolome_aa_mock <- subset(metabolome_aa, infection == 'mock')
metabolome_aa_mock$infection <- NULL
metabolome_aa_mock_strep <- subset(metabolome_aa_mock, abx == 'streptomycin')
metabolome_aa_mock_strep$abx <- NULL
metabolome_aa_mock_cef <- subset(metabolome_aa_mock, abx == 'cefoperazone')
metabolome_aa_mock_cef$abx <- NULL
metabolome_aa_mock_clinda <- subset(metabolome_aa_mock, abx == 'clindamycin')
metabolome_aa_mock_clinda$abx <- NULL
rm(metabolome_aa_mock)
metabolome_aa_infected <- subset(metabolome_aa, infection == '630')
metabolome_aa_infected$infection <- NULL
metabolome_aa_infected_strep <- subset(metabolome_aa_infected, abx == 'streptomycin')
metabolome_aa_infected_strep$abx <- NULL
metabolome_aa_infected_cef <- subset(metabolome_aa_infected, abx == 'cefoperazone')
metabolome_aa_infected_cef$abx <- NULL
metabolome_aa_infected_clinda <- subset(metabolome_aa_infected, abx == 'clindamycin')
metabolome_aa_infected_clinda$abx <- NULL
rm(metabolome_aa_infected)

metabolome_carboxy_mock <- subset(metabolome_carboxy, infection == 'mock')
metabolome_carboxy_mock$infection <- NULL
metabolome_carboxy_mock_strep <- subset(metabolome_carboxy_mock, abx == 'streptomycin')
metabolome_carboxy_mock_strep$abx <- NULL
metabolome_carboxy_mock_cef <- subset(metabolome_carboxy_mock, abx == 'cefoperazone')
metabolome_carboxy_mock_cef$abx <- NULL
metabolome_carboxy_mock_clinda <- subset(metabolome_carboxy_mock, abx == 'clindamycin')
metabolome_carboxy_mock_clinda$abx <- NULL
rm(metabolome_carboxy_mock)
metabolome_carboxy_infected <- subset(metabolome_carboxy, infection == '630')
metabolome_carboxy_infected$infection <- NULL
metabolome_carboxy_infected <- subset(metabolome_carboxy, infection == '630')
metabolome_carboxy_infected$infection <- NULL
metabolome_carboxy_infected_strep <- subset(metabolome_carboxy_infected, abx == 'streptomycin')
metabolome_carboxy_infected_strep$abx <- NULL
metabolome_carboxy_infected_cef <- subset(metabolome_carboxy_infected, abx == 'cefoperazone')
metabolome_carboxy_infected_cef$abx <- NULL
metabolome_carboxy_infected_clinda <- subset(metabolome_carboxy_infected, abx == 'clindamycin')
metabolome_carboxy_infected_clinda$abx <- NULL
rm(metabolome_carboxy_infected)

rm(metadata, metabolome_aa, metabolome_keto, metabolome_carboxy)

#-------------------------------------------------------------------------------------------------------------------------#




# Move inside function
# Find significant differences
infection_aa_pval <- c()
for (i in 1:ncol(metabolome_aa_mock)){infection_aa_pval[i] <- wilcox.test(metabolome_aa_mock[,i], metabolome_aa_infected[,i], exact=FALSE)$p.value}
infection_aa_pval <- round(p.adjust(infection_aa_pval, method='BH'), 3)
infection_carboxy_pval <- c()
for (i in 1:ncol(metabolome_carboxy_mock)){infection_carboxy_pval[i] <- wilcox.test(metabolome_carboxy_mock[,i], metabolome_carboxy_infected[,i], exact=FALSE)$p.value}
infection_carboxy_pval <- round(p.adjust(infection_carboxy_pval, method='BH'), 3)

res_mock_aa_pval <- c()
for (i in 1:ncol(metabolome_aa_mock)){res_mock_aa_pval[i] <- wilcox.test(metabolome_aa_mock[,i], resistant_metabolome_aa[,i], exact=FALSE)$p.value}
res_mock_aa_pval <- round(p.adjust(res_mock_aa_pval, method='BH'), 3)
res_630_aa_pval <- c()
for (i in 1:ncol(metabolome_aa_infected)){res_630_aa_pval[i] <- wilcox.test(metabolome_aa_infected[,i], resistant_metabolome_aa[,i], exact=FALSE)$p.value}
res_630_aa_pval <- round(p.adjust(res_630_aa_pval, method='BH'), 3)
res_mock_carboxy_pval <- c()
for (i in 1:ncol(metabolome_carboxy_mock)){res_mock_carboxy_pval[i] <- wilcox.test(metabolome_carboxy_mock[,i], resistant_metabolome_carboxy[,i], exact=FALSE)$p.value}
res_mock_carboxy_pval <- round(p.adjust(res_mock_carboxy_pval, method='BH'), 3)
res_630_carboxy_pval <- c()
for (i in 1:ncol(metabolome_carboxy_infected)){res_630_carboxy_pval[i] <- wilcox.test(metabolome_carboxy_infected[,i], resistant_metabolome_carboxy[,i], exact=FALSE)$p.value}
res_630_carboxy_pval <- round(p.adjust(res_630_carboxy_pval, method='BH'), 3)

x <- 1
proline_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                  res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
hydroxyproline_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                         res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
glycine_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                  res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
alanine_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                  res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
acetylalanine_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                        res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
valine_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                 res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
acetylvaline_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                       res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
leucine_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                  res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
acetylleucine_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                        res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
isoleucine_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                     res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
acetylisoleucine_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                           res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
phenylalanine_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                        res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
acetylphenylalanine_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                              res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
acetyltyrosine_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                         res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
isocaproate_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                      res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
lactate_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                  res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
x <- x + 1
aminovalerate_pval <- c(infection_aa_pval[x],infection_carboxy_pval[x],res_mock_aa_pval[x],
                        res_630_aa_pval[x],res_mock_carboxy_pval[x],res_630_carboxy_pval[x])
rm(infection_aa_pval,infection_carboxy_pval,res_mock_aa_pval,
   res_630_aa_pval,res_mock_carboxy_pval,res_630_carboxy_pval, x)

#-------------------------------------------------------------------------------------------------------------------------#

# Function for repeated plotting
metabolitePlot <- function(resistant, 
                           strep_mock, strep_630,
                           cef_mock, cef_630,
                           clinda_mock, clinda_630,
                           index){
  # Get metabolite name
  metabolite <- colnames(resistant)[index]
  # Find y-maximum
  yLimit <- round(max(max(resistant[,index]), 
                      max(strep_mock[,index]), max(strep_630[,index]),
                      max(cef_mock[,index]), max(cef_630[,index]), 
                      max(clinda_mock[,index]), max(clinda_630[,index]))) + 3
  # Calculate p-values
  res_strep_mock_pval <- round(wilcox.test(resistant[,index], strep_mock[,index], exact=FALSE)$p.value, 3)
  res_strep_630_pval <- round(wilcox.test(resistant[,index], strep_630[,index], exact=FALSE)$p.value, 3)
  strep_mock_630_pval <- round(wilcox.test(strep_mock[,index], strep_630[,index], exact=FALSE)$p.value, 3)
  res_cef_mock_pval <- round(wilcox.test(resistant[,index], cef_mock[,index], exact=FALSE)$p.value, 3)
  res_cef_630_pval <- round(wilcox.test(resistant[,index], cef_630[,index], exact=FALSE)$p.value, 3)
  cef_mock_630_pval <- round(wilcox.test(cef_mock[,index], cef_630[,index], exact=FALSE)$p.value, 3)
  res_clinda_mock_pval <- round(wilcox.test(resistant[,index], clinda_mock[,index], exact=FALSE)$p.value, 3)
  res_clinda_630_pval <- round(wilcox.test(resistant[,index], clinda_630[,index], exact=FALSE)$p.value, 3)
  clinda_mock_630_pval <- round(wilcox.test(clinda_mock[,index], clinda_630[,index], exact=FALSE)$p.value, 3)
  
  par(mar=c(3,4,1,1), xpd=FALSE, las=1, mgp=c(3,0.7,0))
  stripchart(resistant[,index], at=0.5, vertical=T, pch=21, lwd=2,
             xaxt='n', bg=noabx_col, ylim=c(0,yLimit), xlim=c(0.25,5.25),
             cex=1.5, ylab='', method='jitter', jitter=0.1)
  
  stripchart(strep_mock[,index], at=1.5, vertical=T, pch=21, lwd=2,
             xaxt='n', bg=strep_col, ylim=c(0,yLimit), xlim=c(0.25,5.25),
             cex=1.5, ylab='', method='jitter', jitter=0.1, add=TRUE)
  stripchart(strep_630[,index], at=2, vertical=T, pch=21, lwd=2,
             xaxt='n', bg=strep_col, ylim=c(0,yLimit), xlim=c(0.25,5.25),
             cex=1.5, ylab='', method='jitter', jitter=0.1, add=TRUE)
  
  stripchart(cef_mock[,index], at=3, vertical=T, pch=21, lwd=2,
             xaxt='n', bg=cef_col, ylim=c(0,yLimit), xlim=c(0.25,5.25),
             cex=1.5, ylab='', method='jitter', jitter=0.1, add=TRUE)
  stripchart(cef_630[,index], at=3.5, vertical=T, pch=21, lwd=2,
             xaxt='n', bg=cef_col, ylim=c(0,yLimit), xlim=c(0.25,5.25),
             cex=1.5, ylab='', method='jitter', jitter=0.1, add=TRUE)
  
  stripchart(clinda_mock[,index], at=4.5, vertical=T, pch=21, lwd=2,
             xaxt='n', bg=clinda_col, ylim=c(0,yLimit), xlim=c(0.25,5.25),
             cex=1.5, ylab='', method='jitter', jitter=0.1, add=TRUE)
  stripchart(clinda_630[,index], at=5.5, vertical=T, pch=21, lwd=2,
             xaxt='n', bg=clinda_col, ylim=c(0,yLimit), xlim=c(0.25,5.25),
             cex=1.5, ylab='', method='jitter', jitter=0.1, add=TRUE)
  abline(v=c(1,2.5,4), lty=5)
  box()
  mtext(text=expression(paste('Scaled Intensity (',log[10],')')), side=2, cex=1.1, las=0, padj=-2.5)
  legend('topright', legend=metabolite, pt.cex=0, bty='n', cex=0.9)

  mtext(c('CDI:','Group:'), side=1, at=-0.1, padj=c(0.3,2.5), cex=0.7, xpd=TRUE)
  mtext(c('-','-','+','-','+','-','+'), side=1, 
        at=c(0.5,1.5,2,3,3.5,4.5,5), padj=0.3, cex=1.1)
  mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin'), side=1, 
        at=c(0.5,1.75,3.25,4.75), padj=2, cex=0.9)

  # Medians
  segments(x0=c(0.3,1.3,1.8,2.8,3.3,4.3,4.8), x1=c(0.7,1.7,2.2,3.2,3.7,4.7,5.2),
           y0=c(median(resistant[,index]),
                median(subset(aminovalerate_strep, infection=='mock')[,2]), median(subset(aminovalerate_strep, infection=='infected')[,2]),
                median(subset(aminovalerate_cef, infection=='mock')[,2]), median(subset(aminovalerate_cef, infection=='infected')[,2]),
                median(subset(aminovalerate_clinda, infection=='mock')[,2]), median(subset(aminovalerate_clinda, infection=='infected')[,2])),
           y1=c(median(resistant[,index]),
                median(strep_mock[,index]), median(strep_630[,index]),
                median(cef_mock[,index]), median(cef_630[,index]),
                median(clinda_mock[,index]), median(clinda_630[,index])),
           lwd=3)
  segments(x0=c(1.5,3,4.5), y0=yLimit-2, x1=c(2,3.5,5), y1=yLimit-2, lwd=2)
  
  # Add significance
  #text(x=c(3.5,6.5,9.5), y=5.2, '*', font=2, cex=2)
  #mtext(c(rep('*',5),'n.s.'), side=3, adj=c(0.25,0.35,
  #                                          0.54,0.64,
  #                                          0.83,0.95), padj=c(rep(0.4,5),-0.1), font=2, cex=c(rep(1.6,5),1.1), col='chartreuse3') # Untreated vs Mock significance

}



# Plot the figure
pdf(file=fig3_plot, width=12, height=6)
layout(matrix(c(1,2),
              nrow=1, ncol=2, byrow=TRUE))



metabolitePlot(resistant, 
               strep_mock, strep_630,
               cef_mock, cef_630,
               clinda_mock, clinda_630,
               1)




# Aminovalerate
pdf(file=plot_g, width=6, height=4)
par(mar=c(3.5,5,1.5,1), xpd=FALSE, las=1, mgp=c(3,0.7,0))
stripchart(substrate~infection, data=aminovalerate_untreated, vertical=T, pch=19, 
           xaxt='n', yaxt='n', col='gray40', ylim=c(0,6), xlim=c(0.5,10.5),
           cex=1.5, ylab='', method='jitter', jitter=0.15, cex.lab=1.2)
stripchart(substrate~infection, data=aminovalerate_strep, vertical=T, pch=19, at=c(3,4),
           xaxt='n', yaxt='n', col=strep_col, ylim=c(0,6), xlim=c(0.5,10.5),
           cex=1.5, ylab='', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=aminovalerate_cef, vertical=T, pch=19, at=c(6,7),
           xaxt='n', yaxt='n', col=cef_col, ylim=c(0,6), xlim=c(0.5,10.5),
           cex=1.5, ylab='', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
stripchart(substrate~infection, data=aminovalerate_clinda, vertical=T, pch=19, at=c(9,10),
           xaxt='n', yaxt='n', col=clinda_col, ylim=c(0,6), xlim=c(0.5,10.5),
           cex=1.5, ylab='', method='jitter', jitter=0.15, cex.lab=1.2, add=TRUE)
axis(side=2, at=c(0:6), labels=c('0.0','1.0','2.0','3.0', '4.0','5.0','6.0'), cex.axis=1.2)
box()
mtext(text=expression(paste('Scaled Intensity (',log[10],')')), side=2, cex=1.2, las=0, padj=-2.5)
abline(v=c(2,5,8,11), lty=2, col='gray35')
mtext(c('CDI:','Group:'), side=1, at=-0.7, padj=c(0.3,2.5), cex=0.7)
mtext(c('-','-','+','-','+','-','+'), side=1, 
      at=c(1,3,4,6,7,9,10), padj=0.3, cex=1.1)
mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin'), side=1, 
      at=c(1,3.5,6.5,9.5), padj=2, cex=0.9)
legend('topright', legend='5-aminovalerate', pt.cex=0, bty='n', cex=0.9)
segments(x0=c(0.6,2.6,3.6,5.6,6.6,8.6,9.6), x1=c(1.4,3.4,4.4,6.4,7.4,9.4,10.4),
         y0=c(median(aminovalerate_untreated[,2]),
              median(subset(aminovalerate_strep, infection=='mock')[,2]), median(subset(aminovalerate_strep, infection=='infected')[,2]),
              median(subset(aminovalerate_cef, infection=='mock')[,2]), median(subset(aminovalerate_cef, infection=='infected')[,2]),
              median(subset(aminovalerate_clinda, infection=='mock')[,2]), median(subset(aminovalerate_clinda, infection=='infected')[,2])),
         y1=c(median(aminovalerate_untreated[,2]),
              median(subset(aminovalerate_strep, infection=='mock')[,2]), median(subset(aminovalerate_strep, infection=='infected')[,2]),
              median(subset(aminovalerate_cef, infection=='mock')[,2]), median(subset(aminovalerate_cef, infection=='infected')[,2]),
              median(subset(aminovalerate_clinda, infection=='mock')[,2]), median(subset(aminovalerate_clinda, infection=='infected')[,2])),
         lwd=3)
segments(x0=c(3,6,9), y0=5, x1=c(4,7,10), y1=5, lwd=2)
text(x=c(3.5,6.5,9.5), y=5.2, '*', font=2, cex=2)
mtext(c(rep('*',5),'n.s.'), side=3, adj=c(0.25,0.35,
                                          0.54,0.64,
                                          0.83,0.95), padj=c(rep(0.4,5),-0.1), font=2, cex=c(rep(1.6,5),1.1), col='chartreuse3') # Untreated vs Mock significance
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
