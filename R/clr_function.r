#  invocation:
#  use selex dataset from ALDEx2 library
#  x <- aldex.clr( reads, conditions, mc.samples=128, denom="all", verbose=FALSE, useMC=FALSE )
#  this function generates the centre log-ratio transform of Monte-Carlo instances
#  drawn from the Dirichlet distribution.

aldex.clr.function <- function( reads, conds, mc.samples=128, denom="all", verbose=FALSE, useMC=FALSE, summarizedExperiment=NULL ) {

# INPUT
# The 'reads' data.frame MUST have row
# and column names that are unique, and
# looks like the following:
#
#              T1a T1b  T2  T3  N1  N2
#   Gene_00001   0   0   2   0   0   1
#   Gene_00002  20   8  12   5  19  26
#   Gene_00003   3   0   2   0   0   0
#       ... many more rows ...
#
# ---------------------------------------------------------------------

# OUTPUT
# The output returned is a list (x) that contains Monte-Carlo instances of
# the centre log-ratio transformed values for each sample
# Access to values
# sample IDs: names(x)
# number of features (genes, OTUs): length(x[[1]][,1])
# number of Monte-Carlo Dirichlet instances: length(x[[1]][1,])
# feature names: rownames(x[[1]])

    # Fully validate and coerce the data into required formats
	# coerce SummarizedExperiment reads into data.frame
	if (summarizedExperiment) {
		reads <- data.frame(as.list(assays(reads,withDimnames=TRUE)))
		if (verbose) {
			message("converted SummarizedExperiment read count object into data frame")
		}
	}

  if(missing(conds)){

    message("no conditions provided: forcing denom = 'all'")
    message("no conditions provided: forcing conds = 'NA'")
    denom <- "all"
    conds <- rep("NA", ncol(reads))

  }#else{

  #   # add special handling for model.matrix input
  #   if(class(conds) == "matrix"){
  #     message("conditions provided as matrix: selecting first column for aldex.clr")
  #     conds <- conds[,1]
  #   }
  #
  #   # ncol df and length(c) must be equal
  #   if(ncol(reads) != length(conds)){
  #     stop("mismatch between number of samples and condition vector")
  #   }
  #
  #   # reorder the samples and conditions by level
  #   conds <- conds[order(conds)]
  #   reads <- data.frame(reads[,order(conds)])
  # }

    # make sure that the multicore package is in scope and return if available
    has.BiocParallel <- FALSE
    if ("BiocParallel" %in% rownames(installed.packages()) & useMC){
        message("multicore environment is is OK -- using the BiocParallel package")
        #require(BiocParallel)
        has.BiocParallel <- TRUE
    }
    else {
        message("operating in serial mode")
    }

    # make sure that mc.samples is an integer, despite it being a numeric type value
    mc.samples <- as.numeric(as.integer(mc.samples))

    #  remove all rows with reads less than the minimum set by minsum
    minsum <- 0

    # remove any row in which the sum of the row is 0
    z <- as.numeric(apply(reads, 1, sum))
    reads <- as.data.frame( reads[(which(z > minsum)),]  )

    if (verbose) message("removed rows with sums equal to zero")


    #  SANITY CHECKS ON THE DATA INPUT
    if ( any( round(reads) != reads ) ) stop("not all reads are integers")
    if ( any( reads < 0 ) )             stop("one or more reads are negative")

    for ( col in names(reads) ) {
        if ( any( ! is.finite( reads[[col]] ) ) )  stop("one or more reads are not finite")
    }

    if ( length(rownames(reads)) == 0 ) stop("rownames(reads) cannot be empty")
    if ( length(colnames(reads)) == 0 ) stop("colnames(reads) cannot be empty")

    if ( length(rownames(reads)) != length(unique(rownames(reads))) ) stop ("row names are not unique")
    if ( length(colnames(reads)) != length(unique(colnames(reads))) ) stop ("col names are not unique")
    if ( mc.samples < 128 ) warning("values are unreliable when estimated with so few MC smps")

    # add a prior expection to all remaining reads that are 0
    # this should be by a Count Zero Multiplicative approach, but in practice
    # this is not necessary because of the large number of features
    prior <- 0.5

    # This extracts the set of features to be used in the geometric mean computation
    # returns a list of features
    feature.subset <- aldex.set.mode(reads, conds, denom)

    if ( length(feature.subset[[1]]) == 0 ) stop("No low variance, high abundance features in common between conditions\nPlease choose another denomiator.")

    reads <- reads + prior

if (verbose == TRUE) message("data format is OK")

    # ---------------------------------------------------------------------
    # Generate a Monte Carlo instance of the frequencies of each sample via the Dirichlet distribution,
    # returns frequencies for each feature in each sample that are consistent with the
    # feature count observed as a proportion of the total counts per sample given
    # technical variation (i.e. proportions consistent with error observed when resequencing the same library)

    nr <- nrow( reads )
    rn <- rownames( reads )

    #this returns a list of proportions that are consistent with the number of reads per feature and the
    #total number of reads per sample

    # environment test, runs in multicore if possible
    if (has.BiocParallel){
        p <- bplapply( reads ,
            function(col) {
                q <- t( rdirichlet( mc.samples, col ) ) ;
                rownames(q) <- rn ;
                q })
        names(p) <- names(reads)
    }
    else{
        p <- lapply( reads ,
            function(col) {
                q <- t( rdirichlet( mc.samples, col ) ) ;
                rownames(q) <- rn ; q } )
    }

    # sanity check on the data, should never fail
    for ( i in 1:length(p) ) {
            if ( any( ! is.finite( p[[i]] ) ) ) stop("non-finite frequencies estimated")
    }

if (verbose == TRUE) message("dirichlet samples complete")

    # ---------------------------------------------------------------------
    # Take the log2 of the frequency and subtract the geometric mean log2 frequency per sample
    # i.e., do a centered logratio transformation as per Aitchison

    # apply the function over elements in a list, that contains an array

    # DEFAULT
    if (is.list(feature.subset)) {
        # ZERO only
        feat.result <- vector("list", length(unique(conds))) # Feature Gmeans
        condition.list <- vector("list", length(unique(conds)))    # list to store conditions

        for (i in 1:length(unique(conds)))
        {
            condition.list[[i]] <- which(conds == unique(conds)[i]) # Condition list
            feat.result[[i]] <- lapply( p[condition.list[[i]]], function(m) {
                apply(log2(m), 2, function(x){mean(x[feature.subset[[i]]])})
            })
        }
        set.rev <- unlist(feat.result, recursive=FALSE) # Unlist once to aggregate samples
        p.copy <- p
        for (i in 1:length(set.rev))
        {
            p.copy[[i]] <- as.data.frame(p.copy[[i]])
            p[[i]] <- apply(log2(p.copy[[i]]),1, function(x){ x - (set.rev[[i]])})
            p[[i]] <- t(p[[i]])
        }
        l2p <- p    # Save the set in order to generate the aldex.clr variable
    } else if (is.vector(feature.subset)){
        # Default ALDEx2, iqlr, user defined, lvha
        if (has.BiocParallel){
            l2p <- bplapply( p, function(m) {
                apply( log2(m), 2, function(col) { col - mean(col[feature.subset]) } )
            })
            names(l2p) <- names(p)
        }
        else{
            l2p <- lapply( p, function(m) {
                apply( log2(m), 2, function(col) { col - mean(col[feature.subset]) } )
            })
        }
    }  else {
        message("the denominator is not recognized, use a different denominator")
    }

    # sanity check on data
    for ( i in 1:length(l2p) ) {
        if ( any( ! is.finite( l2p[[i]] ) ) ) stop("non-finite log-frequencies were unexpectedly computed")
    }
if (verbose == TRUE) message("transformation complete")

    return(new("aldex.clr",reads=reads,mc.samples=mc.samples,conds=conds,denom=feature.subset,verbose=verbose,useMC=useMC,analysisData=l2p))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setMethod("getMonteCarloInstances", signature(.object="aldex.clr"), function(.object) .object@analysisData)

setMethod("getSampleIDs", signature(.object="aldex.clr"), function(.object) names(.object@analysisData))

setMethod("getFeatures", signature(.object="aldex.clr"), function(.object) .object@analysisData[[1]][,1])

setMethod("numFeatures", signature(.object="aldex.clr"), function(.object) length(.object@analysisData[[1]][,1]))

setMethod("numMCInstances", signature(.object="aldex.clr"), function(.object) length(.object@analysisData[[1]][1,]))

setMethod("getFeatureNames", signature(.object="aldex.clr"), function(.object) rownames(.object@analysisData[[1]]))

setMethod("getReads", signature(.object="aldex.clr"), function(.object) .object@reads)

setMethod("numConditions", signature(.object="aldex.clr"), function(.object) length(names(.object@analysisData)))

setMethod("getMonteCarloReplicate", signature(.object="aldex.clr",i="numeric"), function(.object,i) .object@analysisData[[i]])

setMethod("getMonteCarloSample", signature(.object="aldex.clr",i="numeric"), function(.object,i) sapply(.object@analysisData, function(x) {x[,i]}))

setMethod("getDenom", signature(.object="aldex.clr"), function(.object) .object@denom)

setMethod("aldex.clr", signature(reads="data.frame"), function(reads, conds, mc.samples=128, denom="all", verbose=FALSE, useMC=FALSE) aldex.clr.function(reads, conds, mc.samples, denom, verbose, useMC, summarizedExperiment=FALSE))

setMethod("aldex.clr", signature(reads="matrix"), function(reads, conds, mc.samples=128, denom="all", verbose=FALSE, useMC=FALSE) aldex.clr.function(as.data.frame(reads), conds, mc.samples, denom, verbose, useMC, summarizedExperiment=FALSE))

setMethod("aldex.clr", signature(reads="RangedSummarizedExperiment"), function(reads, conds, mc.samples=128, denom="all", verbose=FALSE, useMC=FALSE) aldex.clr.function(reads, conds, mc.samples, denom, verbose, useMC, summarizedExperiment=TRUE))
