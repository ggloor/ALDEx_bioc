### -------------------------------------------------------------------------
### aldex.clr
###

setGeneric("getMonteCarloInstances", function(.object) standardGeneric("getMonteCarloInstances"))

setGeneric("getSampleIDs", function(.object) standardGeneric("getSampleIDs"))

setGeneric("numFeatures", function(.object) standardGeneric("numFeatures"))

setGeneric("getFeatures", function(.object) standardGeneric("getFeatures"))

setGeneric("numMCInstances", function(.object) standardGeneric("numMCInstances"))

setGeneric("getFeatureNames", function(.object) standardGeneric("getFeatureNames"))

setGeneric("getReads", function(.object) standardGeneric("getReads"))

setGeneric("numConditions", function(.object) standardGeneric("numConditions"))

setGeneric("getMonteCarloReplicate", function(.object, i=-1) standardGeneric("getMonteCarloReplicate"))

setGeneric("aldex.clr", function(reads, conds, mc.samples=128, denom="all", verbose=FALSE, useMC=FALSE) standardGeneric("aldex.clr"), signature=c("reads"))
