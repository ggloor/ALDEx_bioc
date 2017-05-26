###########################
# a first run at a glm function for ALDEx2 output

# April 14th, 2014

# Initially, I wrote the function to extract model coefficients and 95%CI, but decided that
# we would probably want to do the pairwise comparisons (with proper p-value correction)
# for the features that have significant p-values by glm or Kruskal-Wallis. If people are interested
# in the coefficients, I can reinsert them, but it ups computation time.

# There is probably room for optimization to speed things up, but for the moment this works

##################
# for glms
##################

# Data structure returned by clr_test_function.r function
# The output returned is a list (x) that contains Monte-Carlo instances of
# the centre log-ratio transformed values for each sample
# sample IDs: names(x)
# number of features: length(x[[1]][,1])
# number of Monte-Carlo Dirichlet instances: length(x[[1]][1,])
# feature names: rownames(x[[1]])

# INVOCATION
# conditions is using selex dataset
# conditions <- c(rep("N", 7), rep("S",7)
# x.glm <- aldex.glm(x,conditions)

# returns a dataframe of expected P and fdr statistics for each feature
aldex.glm <- function(clr, conditions, useMC=FALSE){
  
  # make sure that the multicore package is in scope and return if available
  is.multicore = FALSE
  
  if ("BiocParallel" %in% rownames(installed.packages()) & useMC){
    print("multicore environment is OK -- using the BiocParallel package")
    #require(BiocParallel)
    is.multicore = TRUE
  }
  else {
    print("operating in serial mode")
  }
  
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
  
  # generate the comparison sets from the condition levels
  levels <- vector( "list", length( levels ) )
  names( levels ) <- levels( conditions )
  sets <- names(levels)
  
  # set up the glm results containers
  glm.matrix.p <- as.data.frame(matrix(1, nrow = feature.number, ncol = mc.instances))
  glm.matrix.pBH <- glm.matrix.p # duplicate result container
  kw.p.matrix <- glm.matrix.p # duplicate result container
  kw.pBH.matrix <- glm.matrix.p # duplicate result container
  
  # mc.i is the i-th Monte-Carlo instance
  print("running tests for each MC instance:")
  mc.all <- getMonteCarloInstances(clr)
  for(mc.i in 1:mc.instances){
    
    numTicks <- progress(mc.i, mc.instances, numTicks)
    
    # generate a matrix of i-th Monte-Carlo instance, columns are samples, rows are features
    t.input <- sapply(mc.all, function(y){y[, mc.i]})
    
    # do glms on each feature and make a list of glm outputs
    x <-
      apply(t.input, 1, function(yy) {
        glm(as.numeric(yy) ~ factor(conditions))
      })
    
    # calculate p-values for generalized linear model
    if(is.multicore == TRUE){
      pps <- bplapply(x, drop1, test = "Chis")
    }else{
      pps <- lapply(x, drop1, test = "Chis")
    }
    glm.matrix.p[, mc.i] <- sapply(pps, function(x){x[[5]][2]})
    glm.matrix.pBH[, mc.i] <- as.numeric(p.adjust(glm.matrix.p[, mc.i], method = "BH"))
    
    # calculate p-values for Kruskal Wallis test
    kw.p.matrix[, mc.i] <-
      apply(t.input, 1, function(yy){
        kruskal.test(yy, g = factor(conditions))[[3]]
      })
    kw.pBH.matrix[, mc.i] <- as.numeric(p.adjust(kw.p.matrix[, mc.i], method = "BH"))
    
  }
  
  # get the Expected values of p, q and lfdr
  glm.ep <- rowMeans(glm.matrix.p) # rowMeans is faster than apply()!!
  glm.eBH <- rowMeans(glm.matrix.pBH)
  kw.ep <- rowMeans(kw.p.matrix)
  kw.eBH <- rowMeans(kw.pBH.matrix)
  
  z <- data.frame(kw.ep, kw.eBH, glm.ep, glm.eBH)
  rownames(z) <- getFeatureNames(clr)
  return(z)
}
