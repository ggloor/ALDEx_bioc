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

#returns a dataframe of expected P and fdr statistics for each feature

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
    feature.number <- numFeatures(clr)
    mc.instances <- numMCInstances(clr)
    feature.names <- getFeatureNames(clr)

  conditions <- as.factor( conditions )
  levels     <- levels( conditions )

  if ( length( conditions ) !=  numConditions(clr) )  stop("mismatch btw 'length(conditions)' and 'length(names(clr))'")

  levels <- vector( "list", length( levels ) )
  names( levels ) <- levels( conditions )
  sets <- names(levels)
  lsets <- length(sets)

  # set up the glm results containers

  # p-values and BH p-values
  glm.matrix.p <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
  glm.matrix.pBH <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)

  # for Kruskal Wallis p
  kw.p.matrix = matrix(data = NA, nrow = feature.number, ncol = mc.instances)
  kw.pBH.matrix = matrix(data = NA, nrow = feature.number, ncol = mc.instances)

  #########
  # this is where optimization would probably be helpful to possibly avoid this for loop
  #########
  #mc.i is the monte carlo instance
  for(mc.i in 1:mc.instances){
    print(mc.i)

    #generate a matrix of each Monte-Carlo instance, columns are samples, rows are features
    t.input <- sapply(getMonteCarloInstances(clr), function(y){y[,mc.i]})

    # do glms on each feature
    # make a list of glm outputs
    x <- apply(t.input, 1, function(yy) {
      glm(as.numeric(yy) ~ factor(conditions))
    })

    # p-valuess
    if (is.multicore == TRUE){
        pps <- bplapply(x, drop1, test = "Chis"  )
    } else {
        pps <- lapply(x, drop1, test = "Chis")
    }
    glm.matrix.p[, mc.i] <- sapply(pps, function(x){x[[5]][2]})
    glm.matrix.pBH[, mc.i] <- as.numeric(p.adjust(glm.matrix.p[, mc.i], method = "BH"))

    # Kruskal Wallis
    kw.p.matrix[, mc.i] <- t(apply(t.input, 1, function(yy){
      kruskal.test(yy, g = factor(conditions))[[3]]
    }))
    kw.pBH.matrix[, mc.i] <- as.numeric(p.adjust(kw.p.matrix[, mc.i], method = "BH"))

  }

  #get the Expected values of p, q and lfdr
  glm.ep <- apply(glm.matrix.p, 1, mean)
  glm.eBH <- apply(glm.matrix.pBH, 1, mean)

  kw.ep <- apply(kw.p.matrix, 1, mean)
  kw.eBH <- apply(kw.pBH.matrix, 1, mean)

  z <- data.frame(kw.ep, kw.eBH, glm.ep, glm.eBH)
  rownames(z) <- getFeatureNames(clr)
  return(z)

}
