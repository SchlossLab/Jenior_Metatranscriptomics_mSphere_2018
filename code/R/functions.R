
# Set sorking directory
setwd('~/Desktop/Repositories/Jenior_Metatranscriptomics_2016/')

# Load dependencies
deps <- c('vegan', 'shape', 'plotrix', 'reshape2', 'GMD', 'randomForest')
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
noabx_col <- 'gray40'

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
featureselect_RF <- function(training_data, feature){
  
  attach(training_data)
  levels <- as.vector(unique(training_data[,feature]))
  subfactor_1 <- round(length(rownames(training_data[which(training_data[,feature]==levels[1]),])) * 0.623)
  subfactor_2 <- round(length(rownames(training_data[which(training_data[,feature]==levels[2]),])) * 0.623)
  factor <- max(c(round(subfactor_1 / subfactor_2), round(subfactor_2 / subfactor_1))) * 3
  
  # Breiman (2001). Random Forests. Machine Learning.
  n_trees <- round(length(colnames(training_data)) - 1) * factor
  m_tries <- round(sqrt(length(colnames(training_data)) - 1))
  data_randomForest <- randomForest(training_data[,feature]~., 
                                    data=training_data, importance=TRUE, replace=FALSE, 
                                    err.rate=TRUE, ntree=n_trees, mtry=m_tries)
  detach(training_data)
  
  # Examine OOB error
  print(data_randomForest)
  
  # Parse features for significance and sort
  features_RF <- importance(data_randomForest, type=1)
  final_features_RF <- subset(features_RF, features_RF > abs(min(features_RF)))
  final_features_RF <- final_features_RF[!(rownames(final_features_RF) == feature),]
  final_features_RF <- as.data.frame(final_features_RF)
  
  return(final_features_RF)
}

# get KEGG organism codes from gene annotation row names
get_kegg_org <- function(mappings){
  genes <- strsplit(rownames(mappings), ':')
  orgs <- lapply(genes, `[[`, 1)
  mappings$org_code <- orgs
  return(mappings)
}

