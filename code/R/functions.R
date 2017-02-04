
# Neatly merge 2 matices with shared row names
clean_merge <- function(data_1, data_2){
  
  clean_merged <- merge(data_1, data_2, by = 'row.names')
  rownames(clean_merged) <- clean_merged$Row.names
  clean_merged$Row.names <- NULL
  
  return(clean_merged)
}


# Plot logarithmic tick marks on axes
minor.ticks.axis <- function(ax,n,t.ratio=0.5,mn,mx,...){
  
  lims <- par("usr")
  if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]
  
  major.ticks <- pretty(lims,n=5)
  if(missing(mn)) mn <- min(major.ticks)
  if(missing(mx)) mx <- max(major.ticks)
  
  major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]
  
  axis(ax,at=major.ticks, las=1)
  
  n <- n+2
  minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
  minors <- minors[-c(1,n)]
  
  minor.ticks = c(outer(minors,major.ticks,`+`))
  minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]
  
  
  axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
}


# Feature selection with Random Forest
featureselect_RF <- function(training_data, feature){
  
  # Load package
  if ('randomForest' %in% installed.packages()[,'Package'] == FALSE){
    install.packages('randomForest', quiet=TRUE)}
  library('randomForest', verbose=FALSE, character.only=TRUE)
  
  # Set parameters
  set.seed(6189)
  attach(training_data)
  levels <- as.vector(unique(training_data[,feature]))
  subfactor_1 <- round(length(rownames(training_data[which(training_data[,feature]==levels[1]),])) * 0.623)
  subfactor_2 <- round(length(rownames(training_data[which(training_data[,feature]==levels[2]),])) * 0.623)
  factor <- max(c(round(subfactor_1 / subfactor_2), round(subfactor_2 / subfactor_1))) * 3
  
  # Breiman (2001). Random Forests. Machine Learning.
  n_tree <- round(length(colnames(training_data)) - 1) * factor
  m_try <- round(sqrt(length(colnames(training_data)) - 1))
  data_randomForest <- randomForest(training_data[,feature]~., 
                                    data=training_data, importance=TRUE, replace=FALSE, 
                                    do.trace=100, err.rate=TRUE, ntree=n_tree, mtry=m_try)
  detach(training_data)
  
  # Parse features for significance and sort
  features_RF <- importance(data_randomForest, type=1)
  final_features_RF <- subset(features_RF, features_RF > abs(min(features_RF)))
  final_features_RF <- final_features_RF[!(rownames(final_features_RF) == feature),]
  final_features_RF <- as.data.frame(final_features_RF)
  
  return(final_features_RF)
}