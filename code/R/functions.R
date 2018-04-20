
# Set sorking directory
setwd('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/')

# Load dependencies
deps <- c('vegan', 'shape', 'plotrix', 'reshape2', 'GMD', 'randomForest', 'AUCRF','RColorBrewer', 'gplots','viridis', 'scales','Hmisc','VennDiagram')
for (dep in deps){
  if (dep %in% installed.packages()[,"Package"] == FALSE){
    install.packages(as.character(dep), quiet=TRUE);
  } 
  library(dep, verbose=FALSE, character.only=TRUE)
}
rm(dep)

# Set seed for RNG
set.seed(8619)

# Conserved colors across studies and figures
strep_col <- '#D37A1F'
cef_col <- '#3A9CBC'
clinda_col <- '#A40019'
noabx_col <- 'gray30'
gf_col <- 'forestgreen'
heat_palette <- viridis(n=200)

# Calculate variances of rows
RowVar <- function(x) {
  rowSums((x - rowMeans(x))^2)/(dim(x)[2] - 1)
}

# Remove columns with low variance
rm_lowVar <- function(data_table) {
  vars <- apply(data_table, 2, var)
  keep <- which(vars > as.vector(quantile(vars, na.rm = TRUE))[2])
  keep_table <- data_table[,keep]
  return(keep_table)
}

# Filter out columns that have values in at least 3 samples (ignores first column if needed)
filter_table <- function(data) {
  
  drop <- c()
  if (class(data[,1]) != 'character') {
    if (sum(data[,1] != 0) < 3) {
      drop <- c(drop, colnames(data)[1])
    }
  }
  
  for (index in 2:ncol(data)) {
    if (sum(data[,index] != 0) < 3) {
      drop <- c(drop, colnames(data)[index])
    }
  }
  
  filtered_data <- data[,!(colnames(data) %in% drop)]
  return(filtered_data)
}

# Calculate median and 95% confidence for non-normal, large datasets
# Based on: Conover, W.J. (1980) Practical Nonparametric Statistics John Wiley and Sons, New York.
conf_interval <- function(data) {
  
  data_median <- median(data)
  data <- sort(unique(data))
  n <- length(data)
  q <- 0.5
  nq <- n * q
  conf_range <- 1.96 * sqrt(n * q * (1 - q))
  j <- ceiling(nq - conf_range)
  k <- ceiling(nq + conf_range)
  lower_95 <- data[j]
  upper_95 <- data[k]
  
  return(c(lower_95, data_median, upper_95))
}

# Neatly merge 2 matices with shared row names
clean_merge <- function(data_1, data_2){
  
  clean_merged <- merge(data_1, data_2, by='row.names', all.y=TRUE)
  rownames(clean_merged) <- clean_merged$Row.names
  clean_merged$Row.names <- NULL
  
  return(clean_merged)
}

# Calculates distance of a point from x=y in 2-d space
dist_xy <- function(x) {
  v1 <- c(20,20) - x
  v2 <- c(0,0) - c(20,20)
  m <- cbind(v1, v2)
  distance <- abs(det(m)) / sqrt(sum(v1 * v1))
  return(distance)
} 

# Stepwise rarefaction analysis (order of magnitude)
stepRarefy <- function(abundVect){
  options(warn=-1)
  rareVect <- c()
  subVect <- c(1)
  sub_max <- 10^ceiling(log10(signif(as.numeric(sum(abundVect)), digits=1)))
  current <- 1
  while (current != sub_max)
  {
    current <- current * 10
    subVect <- c(subVect, current)
  }
  for (x in 1:length(subVect)) {
    rareVect[x] <- as.numeric(sum(as.vector(rrarefy(abundVect, sample=subVect[x])) != 0))
  }
  options(warn=0)
  return(rareVect)
}

# Plot logarithmic tick marks on axes
minor.ticks.axis <- function(ax, n, t.ratio=0.5, mn, mx,...){
  
  lims <- par("usr")
  if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]
  
  major.ticks <- pretty(lims,n=5)
  if(missing(mn)) mn <- min(major.ticks)
  if(missing(mx)) mx <- max(major.ticks)
  
  major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]
  
  axis(ax, at=major.ticks, labels=rep('', length(major.ticks)), las=1)
  
  n <- n+2
  minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
  minors <- minors[-c(1,n)]
  
  minor.ticks = c(outer(minors,major.ticks,`+`))
  minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]
  
  axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
}

abxRF <- function(training_data){
  
  mTries <- round(sqrt(ncol(training_data) - 1))
  
  # Finish reformatting data
  abx <- training_data$abx
  abx <- as.factor(abx)
  abx <- droplevels(abx)
  training_data$abx <- NULL
  
  # Run random forest and get MDA values
  set.seed(906801)
  modelRF <- randomForest(abx ~ ., data=training_data, 
                          importance=TRUE, replace=FALSE, err.rate=TRUE, mtry=mTries)
  #featRF <- importance(modelRF, type=1)
  #featRF <- as.data.frame(featRF)
  #featRF$features <- rownames(featRF)
  return(modelRF)
  
  # Filter to significant features (Strobl 2002) and sort
  #sigFeatRF <- as.data.frame(subset(featRF, 
  #                                  featRF$MeanDecreaseAccuracy > (abs(min(featRF$MeanDecreaseAccuracy)))))
  #sigFeatRF <- sigFeatRF[order(-sigFeatRF$MeanDecreaseAccuracy),]
  # Subset trainng data to significant features
  #finalFeat <- training_data[,which(rownames(sigFeatRF) %in% colnames(training_data))]
  #rm(training_data, sigFeatRF)
  #return(sigFeatRF)
}

# 2 levels only!
aucrfInfection <- function(training_data){
  
  # Format levels of infection for AUCRF
  levels <- as.vector(unique(training_data$infection))
  training_data$infection <- as.character(training_data$infection)
  training_data$infection[which(training_data$infection==levels[1])] <- 0
  training_data$infection[which(training_data$infection==levels[2])] <- 1
  training_data$infection <- as.factor(as.numeric(training_data$infection))
  rm(levels)
  
  # Run AUCRF with reproduceable parameters
  set.seed(906801)
  data_RF <- AUCRF(infection ~ ., data=training_data, pdel=0.05, k0=5, ranking='MDA')
  return(data_RF)
}

# 2 levels only!
aucrfSusceptibility <- function(training_data){
  
  # Format levels of susceptibility for AUCRF
  colnames(training_data) <- make.names(colnames(training_data))
  levels <- as.vector(unique(training_data$susceptibility))
  training_data$susceptibility <- as.character(training_data$susceptibility)
  training_data$susceptibility[which(training_data$susceptibility==levels[1])] <- 0
  #print(paste('0 = ', levels[1], sep=''))
  training_data$susceptibility[which(training_data$susceptibility==levels[2])] <- 1
  #print(paste('1 = ', levels[2], sep=''))
  training_data$susceptibility <- as.factor(as.numeric(training_data$susceptibility))
  rm(levels)
  
  # Run AUCRF with reproduceable parameters
  set.seed(906801)
  data_RF <- AUCRF(susceptibility ~ ., data=training_data, pdel=0.05, k0=5, ranking='MDA')
  return(data_RF)
}

# Reads and formats network data for plotting
format_network <- function(community_importances, color_pallette){
  community_importances <- community_importances[!community_importances$importance == 0.0,]
  community <- community_importances[,1:2]
  
  # Format metabolic network data
  raw_graph <- graph.data.frame(community, directed=TRUE)
  simple_graph <- simplify(raw_graph)
  
  # Break importances into categorical variable that corresponds to 1-10 rgb scale
  color_increment <- max(community_importances$importances) / 10
  edge_colors <- c()
  edge_sizes <- c()
  col_index <- 1
  for (index in community_importances$importances){
    if (index < color_increment) {
      edge_colors[col_index] <- '#001AE5'
      sizes[col_index] <- 0.05
      col_index <- col_index + 1
    }
    else if (index < (2 * color_increment)) {
      edge_colors[col_index] <- '#1317CE'
      sizes[col_index] <- 0.15
      col_index <- col_index + 1
    }
    else if (index < (3 * color_increment)) {
      edge_colors[col_index] <- '#2614B7'
      sizes[col_index] <- 0.25
      col_index <- col_index + 1
    }
    else if (index < (4 * color_increment)) {
      edge_colors[col_index] <- '#3912A0'
      sizes[col_index] <- 0.35
      col_index <- col_index + 1
    }
    else if (index < (5 * color_increment)) {
      edge_colors[col_index] <- '#4C0F89'
      sizes[col_index] <- 0.45
      col_index <- col_index + 1
    }
    else if (index < (6 * color_increment)) {
      edge_colors[col_index] <- '#5F0D72'
      sizes[col_index] <- 0.55
      col_index <- col_index + 1
    }
    else if (index < (7 * color_increment)) {
      edge_colors[col_index] <- '#720A5B'
      sizes[col_index] <- 0.65
      col_index <- col_index + 1
    }
    else if (index < (8 * color_increment)) {
      edge_colors[col_index] <- '#850744'
      sizes[col_index] <- 0.75
      col_index <- col_index + 1
    }
    else if (index < (9 * color_increment)) {
      edge_colors[col_index] <- '#98052D'
      sizes[col_index] <- 0.85
      col_index <- col_index + 1
    }
    else {
      edge_colors[col_index] <- '#AB0216'
      sizes[col_index] <- 0.95
      col_index <- col_index + 1
    }
  }
  
  # Color the graph
  V(simple_graph)$color <- color_pallette
  E(simple_graph)$color <- edge_colors
  # increase edge size - corresponding with color
  
  # format node names???
  
  return(simple_graph)
}


# Generates plot for significant differences in metabolite concentrations
multiStripchart <- function(plot_file, metabolome1, metabolome2, pvalues, oob, group1, group2, 
                            treatment_col1, treatment_col2, titleStr, titleCol, formattedNames, xLabel){
  
  pdf(file=plot_file, width=4, height=ncol(metabolome1)*1.5)
  layout(matrix(c(1:(ncol(metabolome1)+2)), nrow=(ncol(metabolome1)+2), ncol=1, byrow = TRUE))
  
  par(mar=c(0.2, 0, 0, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE)
  plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
  
  text(x=-4.5, y=-4.5, labels=titleStr, cex=1.2, font=2, col=titleCol) 
  
  legend('bottomright', legend=c(group1, group2), bty='n',
         pt.bg=c(treatment_col1, treatment_col2), pch=21, cex=1.2, pt.cex=2, ncol=2)
  
  par(mar=c(0.2, 2, 0.2, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE, yaxs='i')
  for(i in c(1:(ncol(metabolome1)))){
    xmax <- ceiling(max(c(max(metabolome1[,i]), max(metabolome2[,i]))))
    while(xmax %% 5 != 0 ){xmax <- xmax + 1}
    if (xmax > 1000) {while(xmax %% 100 != 0 ){xmax <- xmax + 1}} else if (xmax > 70){while(xmax %% 10 != 0 ){xmax <- xmax + 1}}
    plot(0, type='n', xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,xmax), ylim=c(0.3,1.8))
    stripchart(at=1.2, jitter(metabolome1[,i], amount=1e-5), 
               pch=21, bg=treatment_col1, method='jitter', jitter=0.12, cex=2, add=TRUE)
    stripchart(at=0.66, jitter(metabolome2[,i], amount=1e-5), 
               pch=21, bg=treatment_col2, method='jitter', jitter=0.12, cex=2, add=TRUE)
    box()
    legend('topright', legend=do.call(expression, formattedNames[i]), pch=1, cex=1.3, pt.cex=0, bty='n')
    if (xmax <= 10) {
      text(x=seq(0,xmax,1), y=0.42, labels=seq(0,xmax,1), cex=1)
      axis(1, at=seq(0,xmax,1), NA, cex.axis=0.8, tck=0.015)
    } else if (xmax > 1000){
      text(x=seq(0,xmax,200), y=0.42, labels=seq(0,xmax,200), cex=1)
      axis(1, at=seq(0,xmax,200), NA, cex.axis=0.8, tck=0.015)
    } else if (xmax > 500){
      text(x=seq(0,xmax,100), y=0.42, labels=seq(0,xmax,100), cex=1)
      axis(1, at=seq(0,xmax,100), NA, cex.axis=0.8, tck=0.015)
    } else if (xmax > 100){
      text(x=seq(0,xmax,50), y=0.42, labels=seq(0,xmax,50), cex=1)
      axis(1, at=seq(0,xmax,50), NA, cex.axis=0.8, tck=0.015)
    } else if (xmax > 50){
      text(x=seq(0,xmax,10), y=0.42, labels=seq(0,xmax,10), cex=1)
      axis(1, at=seq(0,xmax,10), NA, cex.axis=0.8, tck=0.015)
    } else {
      text(x=seq(0,xmax,5), y=0.42, labels=seq(0,xmax,5), cex=1)
      axis(1, at=seq(0,xmax,5), NA, cex.axis=0.8, tck=0.015)
    }
    segments(median(metabolome1[,i]), 1.03, median(metabolome1[,i]), 1.37, lwd=2.5)
    segments(median(metabolome2[,i]), 0.49, median(metabolome2[,i]), 0.83, lwd=2.5)
    if (pvalues[i] <= 0.05){
      mtext('*', side=4, font=2, cex=1.6, padj=0.6)
    } else {
      mtext('n.s.', side=4, cex=0.9)
    }
  }
  
  par(mar=c(0, 0, 0, 0))
  plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
  text(x=0, y=4, labels=xLabel, cex=1.4)
  text(x=8, y=4.5, labels=paste('OOB Error = ', oob, '%',sep=''))
  dev.off()
}

# Function for panels of metabolite differences
metabolitePlot <- function(resistant, 
                           strep_mock, strep_630,
                           cef_mock, cef_630,
                           clinda_mock, clinda_630,
                           index, panelLab){
  # Get metabolite name
  metabolite <- colnames(resistant)[index]
  # Make sure everything is numeric
  resistant <- data.frame(apply(resistant, 2, function(x) as.numeric(as.character(x))))
  strep_mock <- data.frame(apply(strep_mock, 2, function(x) as.numeric(as.character(x))))
  strep_630 <- data.frame(apply(strep_630, 2, function(x) as.numeric(as.character(x))))
  cef_mock <- data.frame(apply(cef_mock, 2, function(x) as.numeric(as.character(x))))
  cef_630 <- data.frame(apply(cef_630, 2, function(x) as.numeric(as.character(x))))
  clinda_mock <- data.frame(apply(clinda_mock, 2, function(x) as.numeric(as.character(x))))
  clinda_630 <- data.frame(apply(clinda_630, 2, function(x) as.numeric(as.character(x))))
  # Find y-maximum
  yLimit <- round(max(max(resistant[,index]), 
                      max(strep_mock[,index]), max(strep_630[,index]),
                      max(cef_mock[,index]), max(cef_630[,index]), 
                      max(clinda_mock[,index]), max(clinda_630[,index]))) + 3
  # Calculate p-values
  res_strep_mock_pval <- round(wilcox.test(resistant[,index], strep_mock[,index], exact=FALSE)$p.value, 3)
  res_strep_630_pval <- round(wilcox.test(resistant[,index], strep_630[,index], exact=FALSE)$p.value, 3)
  res_cef_mock_pval <- round(wilcox.test(resistant[,index], cef_mock[,index], exact=FALSE)$p.value, 3)
  res_cef_630_pval <- round(wilcox.test(resistant[,index], cef_630[,index], exact=FALSE)$p.value, 3)
  res_clinda_mock_pval <- round(wilcox.test(resistant[,index], clinda_mock[,index], exact=FALSE)$p.value, 3)
  res_clinda_630_pval <- round(wilcox.test(resistant[,index], clinda_630[,index], exact=FALSE)$p.value, 3)
  res_pval <- c(res_strep_mock_pval, res_strep_630_pval, res_cef_mock_pval, 
                res_cef_630_pval, res_clinda_mock_pval, res_clinda_630_pval)
  res_pval <- p.adjust(res_pval, method='BH')
  res_strep_mock_pval <- res_pval[1]
  res_strep_630_pval <- res_pval[2]
  res_cef_mock_pval <- res_pval[3]
  res_cef_630_pval <- res_pval[4]
  res_clinda_mock_pval <- res_pval[5]
  res_clinda_630_pval <- res_pval[6]
  strep_mock_630_pval <- round(wilcox.test(strep_mock[,index], strep_630[,index], exact=FALSE)$p.value, 3)
  cef_mock_630_pval <- round(wilcox.test(cef_mock[,index], cef_630[,index], exact=FALSE)$p.value, 3)
  clinda_mock_630_pval <- round(wilcox.test(clinda_mock[,index], clinda_630[,index], exact=FALSE)$p.value, 3)
  
  par(mar=c(3,4,1.5,1), xpd=FALSE, las=1, mgp=c(2.5,0.7,0))
  stripchart(resistant[,index], at=0.5, vertical=T, pch=21, lwd=2,
             xaxt='n', bg=noabx_col, ylim=c(0,yLimit), xlim=c(0.25,5.25),
             cex=2, ylab='', method='jitter', jitter=0.1, cex.axis=1.4)
  
  stripchart(strep_mock[,index], at=1.5, vertical=T, pch=21, lwd=2,
             xaxt='n', bg=strep_col, ylim=c(0,yLimit), xlim=c(0.25,5.25),
             cex=2, ylab='', method='jitter', jitter=0.1, cex.axis=1.4, add=TRUE)
  stripchart(strep_630[,index], at=2, vertical=T, pch=21, lwd=2,
             xaxt='n', bg=strep_col, ylim=c(0,yLimit), xlim=c(0.25,5.25),
             cex=2, ylab='', method='jitter', jitter=0.1, cex.axis=1.4, add=TRUE)
  
  stripchart(cef_mock[,index], at=3, vertical=T, pch=21, lwd=2,
             xaxt='n', bg=cef_col, ylim=c(0,yLimit), xlim=c(0.25,5.25),
             cex=2, ylab='', method='jitter', jitter=0.1, cex.axis=1.4, add=TRUE)
  stripchart(cef_630[,index], at=3.5, vertical=T, pch=21, lwd=2,
             xaxt='n', bg=cef_col, ylim=c(0,yLimit), xlim=c(0.25,5.25),
             cex=2, ylab='', method='jitter', jitter=0.1, cex.axis=1.4, add=TRUE)
  
  stripchart(clinda_mock[,index], at=4.5, vertical=T, pch=21, lwd=2,
             xaxt='n', bg=clinda_col, ylim=c(0,yLimit), xlim=c(0.25,5.25),
             cex=2, ylab='', method='jitter', jitter=0.1, cex.axis=1.4, add=TRUE)
  stripchart(clinda_630[,index], at=5, vertical=T, pch=21, lwd=2,
             xaxt='n', bg=clinda_col, ylim=c(0,yLimit), xlim=c(0.25,5.25),
             cex=2, ylab='', method='jitter', jitter=0.1, cex.axis=1.4, add=TRUE)
  abline(v=c(1,2.5,4), lty=5)
  mtext(text=expression(paste('Scaled Intensity (',log[10],')')), side=2, cex=1.1, las=0, padj=-1.5)
  mtext(panelLab, side=2, line=2, las=2, adj=1, padj=-8, cex=1.6, font=2)
  mtext(c('CDI:','Group:'), side=1, at=-0.3, padj=c(0.35,2.55), cex=0.7, xpd=TRUE)
  mtext(c('-','-','+','-','+','-','+'), side=1, 
        at=c(0.5,1.5,2,3,3.5,4.5,5), padj=0.3, cex=1.1)
  mtext(c('No Antibiotics','Streptomycin','Cefoperazone','Clindamycin'), side=1, 
        at=c(0.5,1.75,3.25,4.75), padj=2, cex=0.8)
  
  # Medians
  segments(x0=c(0.3,1.3,1.8,2.8,3.3,4.3,4.8), x1=c(0.7,1.7,2.2,3.2,3.7,4.7,5.2),
           y0=c(median(resistant[,index]),
                median(strep_mock[,index]), median(strep_630[,index]),
                median(cef_mock[,index]), median(cef_630[,index]),
                median(clinda_mock[,index]), median(clinda_630[,index])),
           y1=c(median(resistant[,index]),
                median(strep_mock[,index]), median(strep_630[,index]),
                median(cef_mock[,index]), median(cef_630[,index]),
                median(clinda_mock[,index]), median(clinda_630[,index])),
           lwd=3)
  segments(x0=c(1.5,3,4.5), y0=yLimit-2, x1=c(2,3.5,5), y1=yLimit-2, lwd=2)
  
  # Add significance
  # VS Resistant
  res_strep_mock_pval[is.nan(res_strep_mock_pval)] <- 0
  if (res_strep_mock_pval <= 0.05) {
    mtext('*', font=2, cex=2, side=3, padj=0.4, adj=0.25, col='chartreuse4')
  } else {
    mtext('n.s.', cex=1, side=3, padj=-0.2, adj=0.25, col='chartreuse4')
  }
  res_strep_630_pval[is.nan(res_strep_630_pval)] <- 0
  if (res_strep_630_pval <= 0.05) {
    mtext('*', font=2, cex=2, side=3, padj=0.4, adj=0.35, col='chartreuse4')
  } else {
    mtext('n.s.', cex=1, side=3, padj=-0.2, adj=0.35, col='chartreuse4')
  }
  res_cef_mock_pval[is.nan(res_cef_mock_pval)] <- 0
  if (res_cef_mock_pval <= 0.05) {
    mtext('*', font=2, cex=2, side=3, padj=0.4, adj=0.54, col='chartreuse4')
  } else {
    mtext('n.s.', cex=1, side=3, padj=-0.2, adj=0.54, col='chartreuse4')
  }
  res_cef_630_pval[is.nan(res_cef_630_pval)] <- 0
  if (res_cef_630_pval <= 0.05) {
    mtext('*', font=2, cex=2, side=3, padj=0.4, adj=0.64, col='chartreuse4')
  } else {
    mtext('n.s.', cex=1, side=3, padj=-0.2, adj=0.64, col='chartreuse4')
  }
  res_clinda_mock_pval[is.nan(res_clinda_mock_pval)] <- 0
  if (res_clinda_mock_pval <= 0.05) {
    mtext('*', font=2, cex=2, side=3, padj=0.4, adj=0.83, col='chartreuse4')
  } else {
    mtext('n.s.', cex=1, side=3, padj=-0.2, adj=0.83, col='chartreuse4')
  }
  res_clinda_630_pval[is.nan(res_clinda_630_pval)] <- 0
  if (res_clinda_630_pval <= 0.05) {
    mtext('*', font=2, cex=2, side=3, padj=0.4, adj=0.93, col='chartreuse4')
  } else {
    mtext('n.s.', cex=1, side=3, padj=-0.2, adj=0.93, col='chartreuse4')
  }
  # Within pretreatments
  strep_mock_630_pval[is.nan(strep_mock_630_pval)] <- 0
  if (strep_mock_630_pval <= 0.05) {
    text(x=1.75, y=yLimit-1.7, '*', font=2, cex=2.3)
  } else {
    text(x=1.75, y=yLimit-1.7, 'n.s.', cex=1.5)
  }
  cef_mock_630_pval[is.nan(cef_mock_630_pval)] <- 0
  if (cef_mock_630_pval <= 0.05) {
    text(x=3.25, y=yLimit-1.7, '*', font=2, cex=2.3)
  } else {
    text(x=3.25, y=yLimit-1.7, 'n.s.', cex=1.5)
  }
  clinda_mock_630_pval[is.nan(clinda_mock_630_pval)] <- 0
  if (clinda_mock_630_pval <= 0.05) {
    text(x=4.75, y=yLimit-1.7, '*', font=2, cex=2.3)
  } else {
    text(x=4.75, y=yLimit-1.7, 'n.s.', cex=1.5)
  }
  legend('topright', legend=metabolite, pt.cex=0, cex=1.3, box.lwd=1, box.col="black", bg="white")
  box(lwd=2)
}

