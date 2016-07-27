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

#returns a dataframe of expected P and fdr statistics for each feature
aldex.ttest <- function(clr, conditions, paired.test=FALSE, hist.plot=FALSE) {

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

    if ( length( conditions ) !=  numConditions(clr) )  stop(paste("mismatch btw 'length(conditions)' and 'length(names(clr))'. len(condtitions):",length(conditions),"len(names(clr)):",numConditions(clr)))

    if ( length( levels ) != 2 ) stop("only two condition levels are currently supported")

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

