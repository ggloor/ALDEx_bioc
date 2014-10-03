####################################
# to add to the aldex options
# a function that runs parametric Pearson's Product moment correlations
# and nonparametric Spearman's rank correlations
# returns the mean cor and rho estimates as well as p-values and BH asjusted p-values
#Arianne Albert, April 23, 2014
#####

#######################################
# a correlation function to correlate a continuous variable with the relative abundances

aldex.corr <- function(clr, covar){ 
  # covar is the continuous variable with which to run correlations
  
  # get dimentions, names, etc from the input data
  smpl.ids <- names(clr)
  feature.number <- length(clr[[1]][,1])
  mc.instances <- length(clr[[1]][1,])
  feature.names <- rownames(clr[[1]])
  
  if ( length( covar ) !=  length(names(clr)) )  stop("mismatch btw 'length(covar)' and 'length(names(clr))'")
  if ( is.numeric( covar ) != TRUE) stop("covar is not numeric")
  
  # set up the correlation results containers
  
  # Pearson correlation cor values, p-values and BH p-values
  pearson.matrix.cor <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
  pearson.matrix.p <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
  pearson.matrix.pBH <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
  
  # Spearman rank correlation rho values, p-values and BH p-values
  spearman.matrix.rho <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
  spearman.matrix.p <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
  spearman.matrix.pBH <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
  
  #mc.i is the monte carlo instance
  for(mc.i in 1:mc.instances){
    print(mc.i)
    
    #generate a matrix of each Monte-Carlo instance, columns are samples, rows are features
    t.input <- sapply(clr, function(y){y[,mc.i]})
    
    # do Pearson correlations on each feature
    # make a list of correlation outputs
    x <- apply(t.input, 1, function(yy) { 
      cor.test(as.numeric(yy), covar)
    })
    
    # p-values
    pearson.matrix.p[, mc.i] <- sapply(x, function(x) x[[3]])
    pearson.matrix.pBH[, mc.i] <- as.numeric(p.adjust(pearson.matrix.p[, mc.i], method = "BH"))
    # cor values
    pearson.matrix.cor[, mc.i] <- sapply(x, function(x) x[[4]])
    
    # Spearman correlations
    x <- apply(t.input, 1, function(yy) { 
      cor.test(as.numeric(yy), covar, method = "spearman")
    })
    
    # p-values
    spearman.matrix.p[, mc.i] <- sapply(x, function(x) x[[3]])
    spearman.matrix.pBH[, mc.i] <- as.numeric(p.adjust(spearman.matrix.p[, mc.i], method = "BH"))
    # rho values
    spearman.matrix.rho[, mc.i] <- sapply(x, function(x) x[[4]])
    
  }
  
  
  #get the Expected values of p, cor, rho, and lfdr
  pearson.ecor <- apply(pearson.matrix.cor, 1, mean)
  pearson.ep <- apply(pearson.matrix.p, 1, mean)
  pearson.eBH <- apply(pearson.matrix.pBH, 1, mean)
  
  spearman.erho <- apply(spearman.matrix.rho, 1, mean)
  spearman.ep <- apply(spearman.matrix.p, 1, mean)
  spearman.eBH <- apply(spearman.matrix.pBH, 1, mean)
  
  z <- data.frame(pearson.ecor, pearson.ep, pearson.eBH, spearman.erho, spearman.ep, spearman.eBH)
  rownames(z) <- rownames(clr[[1]])
  return(z)
  
}