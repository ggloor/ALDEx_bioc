# Data structure returned by aldex.clr function
# The output returned is a list (x) that contains Monte-Carlo instances of
# the centre log-ratio transformed values for each sample
# sample IDs: names(x)
# number of features: length(x[[1]][,1])
# number of Monte-Carlo Dirichlet instances: length(x[[1]][1,])
# feature names: rownames(x[[1]])

# INVOCATION
# conditions is using selex dataset
# conditions <- c(rep("N", 7), rep("S",7)
# x.tt <- aldex.ttest(x,conditions)

# returns a dataframe of expected P and fdr statistics for each feature
aldex.ttest <- function(clr, conditions, paired.test=FALSE, hist.plot=FALSE) {
  
  # get dimensions, names, etc from the input data
  smpl.ids <- getSampleIDs(clr)
  feature.names <- getFeatureNames(clr)
  feature.number <- numFeatures(clr)
  mc.instances <- numMCInstances(clr)
  
  conditions <- as.factor( conditions )
  levels     <- levels( conditions )
  
  if ( length( conditions ) !=  numConditions(clr) ){
    stop(paste("mismatch btw 'length(conditions)' and 'length(names(clr))'. len(condtitions):",
               length(conditions),"len(names(clr)):",numConditions(clr)))}
  if ( length( levels ) != 2 ) stop("only two condition levels are currently supported")
  
  # generate the comparison sets from the condition levels
  levels <- vector( "list", length( levels ) )
  names( levels ) <- levels( conditions )
  sets <- names(levels)
  setAsBinary <- as.numeric(conditions == sets[1])
  setA <- which(conditions == sets[1])
  setB <- which(conditions == sets[2])
  
  # set up the t-test result containers
  wi.p.matrix <- as.data.frame(matrix(1, nrow = feature.number, ncol = mc.instances))
  wi.BH.matrix <- wi.p.matrix # duplicate result container
  we.p.matrix <- wi.p.matrix # duplicate result container
  we.BH.matrix <- wi.p.matrix # duplicate result container
  
  # mc.i is the i-th Monte-Carlo instance
  print("running tests for each MC instance:")
  mc.all <- getMonteCarloInstances(clr)
  for(mc.i in 1:mc.instances){
    
    numTicks <- progress(mc.i, mc.instances, numTicks)
    
    # generate a matrix of i-th Monte-Carlo instance, columns are samples, rows are features
    t.input <- sapply(mc.all, function(y){y[, mc.i]})
    
    wi.p.matrix[, mc.i] <- wilcox.fast(t.input, setAsBinary, paired.test)
    wi.BH.matrix[, mc.i] <- p.adjust(wi.p.matrix[, mc.i], method = "BH")
    
    we.p.matrix[, mc.i] <- t.fast(t.input, setAsBinary, paired.test)
    we.BH.matrix[, mc.i] <- p.adjust(we.p.matrix[, mc.i], method = "BH")
  }
  
  if(hist.plot == TRUE){
    
    par(mfrow=c(2,2))
    hist(we.p.matrix[,1], breaks=99, main="Welch's P values Instance 1")
    hist(wi.p.matrix[,1], breaks=99, main="Wilcoxon P values Instance 1")
    hist(we.BH.matrix[,1], breaks=99, main="Welch's BH values Instance 1")
    hist(wi.BH.matrix[,1], breaks=99, main="Wilcoxon BH values Instance 1")
    par(mfrow=c(1,1))
  }
  
  # get the Expected values of p, q and lfdr
  we.ep <- rowMeans(we.p.matrix) # rowMeans is faster than apply()!!
  we.eBH <- rowMeans(we.BH.matrix)
  wi.ep <- rowMeans(wi.p.matrix)
  wi.eBH <- rowMeans(wi.BH.matrix)
  
  z <- data.frame(we.ep, we.eBH, wi.ep, wi.eBH)
  rownames(z) <- getFeatureNames(clr)
  return(z)
}
