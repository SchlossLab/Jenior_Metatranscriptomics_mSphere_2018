

pca_axes <- function(pcaFile, clusterFile, axesFile){
  
  PCA <- read.csv(pcaFile,header=TRUE,row.names=1)
  Clusters <- read.csv(clusterFile,header=FALSE,row.names=1)
  PCA.df <- data.frame(x=PCA[,1],y=PCA[,2],c=Clusters$V2)
  colnames(PCA.df) <- c('PCA1','PCA2','cluster')
  PCA.df$contig <- rownames(PCA.df)
  
  write.table(PCA.df, file=axesFile, sep='\t', quote=FALSE, row.names=FALSE)
}


strep_clusterFile <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metagenome/concoct/Streptomycin.concoct_output/clustering_gt1000.csv'
strep_pcaFile <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metagenome/concoct/Streptomycin.concoct_output/PCA_transformed_data_gt1000.csv'
pca_axes(strep_pcaFile, strep_clusterFile, '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metagenome/concoct/strep_clusters.tsv')

clinda_clusterFile <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metagenome/concoct/Clindamycin.concoct_output/clustering_gt1000.csv'
clinda_pcaFile <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metagenome/concoct/Clindamycin.concoct_output/PCA_transformed_data_gt1000.csv'
pca_axes(clinda_pcaFile, clinda_clusterFile, '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metagenome/concoct/clinda_clusters.tsv')

cef_clusterFile <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metagenome/concoct/Cefoperazone.concoct_output/clustering_gt1000.csv'
cef_pcaFile <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metagenome/concoct/Cefoperazone.concoct_output/PCA_transformed_data_gt1000.csv'
pca_axes(cef_pcaFile, cef_clusterFile, '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metagenome/concoct/cef_clusters.tsv')

noabx_clusterFile <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metagenome/concoct/Conventional.concoct_output/clustering_gt1000.csv'
noabx_pcaFile <- '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metagenome/concoct/onventional.concoct_output/PCA_transformed_data_gt1000.csv'
pca_axes(noabx_pcaFile, noabx_clusterFile, '~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/data/metagenome/concoct/noabx_clusters.tsv')



