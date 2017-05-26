library(ALDEx2)

# Old aldex.ttest function used to test new function
aldex.ttest.old <- function(clr, conditions, paired.test=FALSE, hist.plot=FALSE) {
  
  # get dimensions, names, etc from the input data
  #   smpl.ids <- names(clr)
  #   feature.number <- length(clr[[1]][,1])
  #   mc.instances <- length(clr[[1]][1,])
  #   feature.names <- rownames(clr[[1]])
  smpl.ids <- getSampleIDs(clr)
  feature.number <- numFeatures(clr)
  mc.instances <- numMCInstances(clr)
  feature.names <- getFeatureNames(clr)
  
  conditions <- as.factor( conditions )
  levels     <- levels( conditions )
  
  levels <- vector( "list", length( levels ) )
  names( levels ) <- levels( conditions )
  sets <- names(levels)
  
  #generate the comparison sets from the condition levels
  setA <- which(conditions == sets[1])
  setB <- which(conditions == sets[2])
  
  # set up the t-test result containers
  we.p.matrix =  matrix(data=NA, nrow = feature.number, ncol = mc.instances)
  we.BH.matrix =  matrix(data=NA, nrow = feature.number, ncol = mc.instances) #benjamini-hochberg
  wi.BH.matrix =  matrix(data=NA, nrow = feature.number, ncol = mc.instances)
  wi.p.matrix =  matrix(data=NA, nrow = feature.number, ncol = mc.instances)
  
  #mc.i is the monte carlo instance
  for(mc.i in 1:mc.instances){
    
    #generate a matrix of each Monte-Carlo instance, columns are samples, rows are features
    t.input <- sapply(getMonteCarloInstances(clr), function(y){y[,mc.i]})
    
    # do the Wilcoxon tests on each feature
    wi.p.matrix[,mc.i] <- t(apply(t.input, 1, function(t.input){as.numeric(wilcox.test(x=t.input[setA],y=t.input[setB])[3])}))
    wi.BH.matrix[,mc.i] <- as.numeric(p.adjust(wi.p.matrix[,mc.i], method="BH"))
    
    # do the welch's test on each feature
    we.p.matrix[,mc.i] <- t(apply(t.input, 1, function(t.input){as.numeric(t.test(x=t.input[setA],y=t.input[setB], paired=paired.test)[3])}))
    we.BH.matrix[,mc.i] <- as.numeric(p.adjust(we.p.matrix[,mc.i], method="BH"))
    
  }
  if (hist.plot == TRUE) {
    par(mfrow=c(2,2))
    hist(we.p.matrix[,1], breaks=99, main="Welch's P values Instance 1")
    hist(wi.p.matrix[,1], breaks=99, main="Wilcoxon P values Instance 1")
    hist(we.BH.matrix[,1], breaks=99, main="Welch's BH values Instance 1")
    hist(wi.BH.matrix[,1], breaks=99, main="Wilcoxon BH values Instance 1")
    par(mfrow=c(1,1))
  }
  #get the Expected values of p, q and lfdr
  we.ep <- apply(we.p.matrix, 1, mean)
  we.eBH <- apply(we.BH.matrix,1,mean)
  wi.ep <- apply(wi.p.matrix, 1, mean)
  wi.eBH <- apply(wi.BH.matrix,1,mean)
  
  z <- data.frame(we.ep, we.eBH, wi.ep, wi.eBH)
  rownames(z) <- getFeatureNames(clr)
  return(z)
}

data(selex)
group <- c(rep("A", 7), rep("B", 7))
clr <- ALDEx2::aldex.clr(selex[1:100,], group, mc.samples = 128)

test_that("new faster alex.ttest matches old function", {
  
  expect_equal(
    aldex.ttest.old(clr, group),
    aldex.ttest(clr, group)
  )
})

data(selex)
dat <- selex
for(i in 1:6){
  fakecol <- data.frame(sample(selex[,1]))
  colnames(fakecol) <- paste0("fake", i)
  dat <- as.data.frame(cbind(dat, fakecol))
}
group <- c(rep("A", 10), rep("B", 10))
clr <- ALDEx2::aldex.clr(dat[1:100,], group, mc.samples = 128)

test_that("new faster alex.ttest matches old function", {
  
  expect_equal(
    aldex.ttest.old(clr, group),
    aldex.ttest(clr, group)
  )
})
