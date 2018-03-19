
# Set sorking directory
setwd('~/Desktop/Repositories/Jenior_Metatranscriptomics_PLOSPathogens_2017/')

# Load dependencies
deps <- c('vegan', 'shape', 'plotrix', 'reshape2', 'GMD', 'randomForest', 'RColorBrewer', 'gplots','viridis', 'scales','Hmisc','VennDiagram')
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


# Feature selection with Random Forest (requires 'randomForest' package)
featureselect_RF <- function(training_data, classes){
  
  attach(training_data)
  levels <- as.vector(unique(training_data[,classes]))
  subfactor_1 <- round(length(rownames(training_data[which(training_data[,classes]==levels[1]),])) * 0.623)
  subfactor_2 <- round(length(rownames(training_data[which(training_data[,classes]==levels[2]),])) * 0.623)
  factor <- max(c(round(subfactor_1 / subfactor_2), round(subfactor_2 / subfactor_1))) * 3
  
  # Breiman (2001). Random Forests. Machine Learning.
  n_trees <- round(length(colnames(training_data)) - 1) * factor
  m_tries <- round(sqrt(length(colnames(training_data)) - 1))
  data_randomForest <- randomForest(classes~., 
                                    data=training_data, importance=TRUE, replace=FALSE, 
                                    err.rate=TRUE, ntree=n_trees, mtry=m_tries)
  detach(training_data)
  
  # Examine OOB error
  print(data_randomForest)
  
  # Parse features for significance and sort
  features_RF <- importance(data_randomForest, type=1)
  final_features_RF <- subset(features_RF, features_RF > abs(min(features_RF)))
  final_features_RF <- final_features_RF[!(rownames(final_features_RF) == classes),]
  final_features_RF <- as.data.frame(final_features_RF)
  final_features_RF$feature <- rownames(final_features_RF)
  colnames(final_features_RF) <- c('MDA','feature')
  
  return(final_features_RF)
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
metabolite_stripchart <- function(plot_file, metabolome1, metabolome2, pvalues, mda, 
                                  oob, group1, group2, fig_title, treatment_col, fig_label){
  
  pdf(file=plot_file, width=4, height=ncol(metabolome1)*1.5)
  layout(matrix(c(1:(ncol(metabolome1)+2)), nrow=(ncol(metabolome1)+2), ncol=1, byrow = TRUE))
  
  par(mar=c(0.2, 0, 0, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE)
  plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
  text(x=-10.2, y=-3, labels=fig_label, cex=2.4, xpd=TRUE, font=2)
  legend('bottomright', legend=c(group1, group2), bty='n',
         pt.bg='black', pch=c(1,16), cex=1.2, pt.cex=2, lwd=c(1,1.5), ncol=2)
  text(x=-4.5, y=-4.5, labels=fig_title, cex=1.2, font=2, col=treatment_col)
  
  par(mar=c(0.2, 2, 0.2, 1), mgp=c(2.3, 0.75, 0), xpd=FALSE, yaxs='i')
  for(i in c(1:(ncol(metabolome1)))){
    xmax <- ceiling(max(c(max(metabolome1[,i]), max(metabolome2[,i]))))
    while(xmax %% 5 != 0 ){xmax <- xmax + 1}
    if (xmax > 1000) {while(xmax %% 100 != 0 ){xmax <- xmax + 1}} else if (xmax > 70){while(xmax %% 10 != 0 ){xmax <- xmax + 1}}
    plot(0, type='n', xlab='', ylab='', xaxt='n', yaxt='n', xlim=c(0,xmax), ylim=c(0.3,1.8))
    stripchart(at=1.2, jitter(metabolome1[,i], amount=1e-5), 
               pch=1, bg=treatment_col, method='jitter', jitter=0.12, cex=2, lwd=2, add=TRUE)
    stripchart(at=0.66, jitter(metabolome2[,i], amount=1e-5), 
               pch=16, bg=treatment_col, method='jitter', jitter=0.12, cex=2, add=TRUE)
    metabolite <- paste(colnames(metabolome1)[i], ' [',as.character(round(mda[i],3)),']', sep='')
    legend('topright', legend=metabolite, pch=1, cex=1.3, pt.cex=0, bty='n')
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
    if (pvalues[i] < 0.001){
      mtext('*', side=4, font=2, cex=1.2, padj=0.5)
    } else if (pvalues[i] <= 0.01){
      mtext('*', side=4, font=2, cex=1.2, padj=0.5)
    } else if (pvalues[i] <= 0.05){
      mtext('*', side=4, font=2, cex=1.2, padj=0.5)
    } else {
      mtext('n.s.', side=4)
    }
  }
 
  par(mar=c(0, 0, 0, 0))
  plot(0, type='n', axes=FALSE, xlab='', ylab='', xlim=c(-10,10), ylim=c(-5,5))
  text(x=0, y=4, labels=expression(paste('Scaled Intensity (',log[10],')')), cex=1.4)
  text(x=8, y=4.5, labels=paste('OOB Error = ', as.character(oob), '%',sep=''))
  
  dev.off()
}

