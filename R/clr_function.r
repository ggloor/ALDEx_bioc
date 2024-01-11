#' Compute an \code{aldex.clr} Object
#' Generate Monte Carlo samples of the Dirichlet distribution for each sample.
#' Convert each instance using a centered log-ratio transform.
#' This is the input for all further analyses.
#' @aliases aldex.clr aldex.clr,data.frame-method aldex.clr,matrix-method aldex.clr,RangedSummarizedExperiment-method
#' @param reads A \code{data.frame} or \code{RangedSummarizedExperiment} object containing
#' non-negative integers only and with unique names for all rows and columns,
#' where each row is a different gene and each column represents a sequencing
#' read-count sample. Rows with 0 reads in each sample are deleted prior to
#' analysis.
#' @param conds A \code{vector} containing a descriptor for the samples, allowing them to
#' be grouped and compared.
#' @param mc.samples The number of Monte Carlo instances to use to estimate the underlying
#' distributions; since we are estimating central tendencies, 128 is usually
#' sufficient, but larger numbers may be needed with small sample sizes.
#' @param denom An \code{any} variable (all, iqlr, zero, lvha, median, user) indicating
#' features to use as the denominator for the Geometric Mean calculation
#' The default "all" uses the geometric mean abundance of all features.
#' Using "median" returns the median abundance of all features.
#' Using "iqlr" uses the features that are between the first and third
#' quartile of the variance of the clr values across all samples.
#' Using "zero" uses the non-zero features in each grop
#' as the denominator. This approach is an extreme case where there are
#' many nonzero features in one condition but many zeros in another. Using
#' "lvha" uses features that have low variance (bottom quartile) and high
#' relative abundance (top quartile in every sample). It is also
#' possible to supply a vector of row indices to use as the denominator.
#' Here, the experimentalist is determining a-priori which rows are thought
#' to be invariant. In the case of RNA-seq, this could include ribosomal
#' protein genes and and other house-keeping genes. This should be used
#' with caution because the offsets may be different in the original data
#' and in the data used by the function because features that are 0 in all
#' samples are removed by \code{aldex.clr}.
#' @param verbose Print diagnostic information while running. Useful only for debugging
#' if fails on large datasets.
#' @param useMC Use multicore by default (FALSE). Multi core processing will be attempted
#' with the BiocParallel package. Serial processing will be used if this is
#' not possible. In practice serial and multicore are nearly the same speed
#' because of overhead in setting up the parallel processes.
#' @param gamma Use scale simulation if not NULL. If a matrix is supplied, scale simulation
#' will be used assuming that matrix denotes the scale samples. If a numeric is
#' supplied, scale simulation will be applied by relaxing the geometric mean
#' assumption with the numeric representing the standard deviation of the
#' scale distribution.
#' @param entropy=FALSE use differential entropy
#' @param summarizedExperiment must be set to TRUE if input data are in this format.
#'
#' @return The object produced by the \code{clr} function contains the log-ratio transformed
#' values for each Monte-Carlo Dirichlet instance, which can be accessed through
#' \code{getMonteCarloInstances(x)}, where \code{x} is the \code{clr} function output.
#' Each list element is named by the sample ID. \code{getFeatures(x)} returns the
#' features, \code{getSampleIDs(x)} returns sample IDs, and \code{getFeatureNames(x)}
#' returns the feature names.
#'
#'    # The 'reads' data.frame or
#'    # RangedSummarizedExperiment object should
#'    # have row and column names that are unique,
#'    # and looks like the following:
#'    #
#'    #              T1a T1b  T2  T3  N1  N2  Nx
#'    #   Gene_00001   0   0   2   0   0   1   0
#'    #   Gene_00002  20   8  12   5  19  26  14
#'    #   Gene_00003   3   0   2   0   0   0   1
#'    #       ... many more rows ...
#'
#'    data(selex)
#'    #subset for efficiency
#'    selex <- selex[1201:1600,]
#'    conds <- c(rep("NS", 7), rep("S", 7))
#'    x <- aldex.clr(selex, conds, mc.samples=4, gamma=NULL, verbose=FALSE)
#' @export
aldex.clr.function <- function( reads, conds, mc.samples=128, denom="all", 
    verbose=FALSE, useMC=FALSE, summarizedExperiment=NULL, gamma = NULL,
    entropy=FALSE) {
#  invocation:
#  use selex dataset from ALDEx2 library
#  x <- aldex.clr( reads, conds, mc.samples=128, denom="all", verbose=FALSE, useMC=FALSE )
#  this function generates the centre log-ratio transform of Monte-Carlo instances
#  drawn from the Dirichlet distribution.

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
  # make sure the conditions vector or matrix is reasonable
  if(missing(conds)){

    message("no conditions provided: forcing denom = 'all'")
    message("no conditions provided: forcing conds = 'NA'")
    denom <- "all"
    conds <- rep("NA", ncol(reads))

  }

# if a model matrix is supplied, then aldex.effect is not valid
# force the use of either all for the denominator
# or
# the use of a user-supplied denominator
  if(is(conds, "vector")){
      message("conditions vector supplied")
  }    
  else if(is(conds, "matrix")  & all(round(conds) == conds)){
    message("integer matrix provided")
    if(is.vector(denom, mode="numeric")){
      message("user-defined denominator used")
    } else if (denom == "all"){
      message("using all features for denominator")
    } else {
      stop("please supply the desired vector of indices for the denominator")
    }
#     if(conds.col == 0){
#       message("conditions provided as matrix: selecting first column for aldex.clr")
#       conds <- as.character(conds[,1])
#     }else{
#       message("conditions provided as matrix: user selected column for aldex.clr")
#       if(is.numeric(conds.col)){
#         print(conds[,conds.col])
#         conds <- as.vector(conds[,conds.col])
#         prints(conds)
#         print(length(conds))
#       }
#     }
  }

  if(ncol(reads) != length(conds) & !is(conds, "matrix")){
    print(length(conds))
    print(ncol(reads))
    stop("mismatch between number of samples and condition vector")
  }

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
    if(is.null(gamma)){
      feature.subset <- aldex.set.mode(reads, conds, denom)
      if ( length(feature.subset[[1]]) == 0 ) stop("No low variance, high abundance features in common between conditions\nPlease choose another denomiator.")
    } else{
      feature.subset <- vector()
    }

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
    # Add scale samples (if desired)
    # Checking the size of the scale samples
    
    if(!is.null(gamma)){
      message("aldex.scaleSim: adjusting samples to reflect scale uncertainty.")
      l2p <- list()
      if(length(gamma) == 1){ ##Add uncertainty around the scale samples
        
        ## grabbing samples from the default scale model
        if(verbose) message("sampling from the default scale model.")
        scale_samples <- default.scale.model(gamma, conds, p, mc.samples)
        ## negative here because calculating CLR
        for(i in 1:length(p)){
          l2p[[i]] <- sweep(log2(p[[i]]), 2,  scale_samples[i,], "-")
          
        }
        scale_samples <- -1*scale_samples #-1 because of the "-" above. This just matters for what is returned.
      } else if(length(gamma) >1 & is.null(dim(gamma))){ ##Vector case/scale sim + senstitivity
        warning("A vector was supplied for scale.samples. To run a sensitivity analysis, use 'aldex.senAnalysis()'.")
        stop("Please supply either a single value or a matrix.")
      } else{ ##User input of scale samples
        Q <- nrow(gamma)
        N <- ncol(gamma)
        if(Q != ncol(reads) | N != mc.samples){
          stop("Scale samples are of incorrect size!")
        }
        if(verbose) message("using user specified scale samples.")
        for(i in 1:length(p)){
        ## positive here because W = W(||) + W(perp)
          l2p[[i]] <- sweep(log2(p[[i]]), 2,  log2(gamma[i,]), "+")
        }
        scale_samples <- log2(gamma)
        
      }
      names(l2p) <- names(p)
    }
    
    # ---------------------------------------------------------------------
    # Take the log2 of the frequency and subtract the geometric mean log2 frequency per sample
    # i.e., do a centered logratio transformation as per Aitchison
    
    # apply the function over elements in a list, that contains an array
    if(is.null(gamma)){
      scale_samples <- NULL
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
        # denom[1] is put in explicitly for the user-defined denominator case
        if (has.BiocParallel){
          if (denom[1] != "median"){
            l2p <- bplapply( p, function(m) {
              apply( log2(m), 2, function(col) { col - mean(col[feature.subset]) } )
             })
          } else if (denom[1] == "median"){
            l2p <- bplapply( p, function(m) {
              apply( log2(m), 2, function(col) { col - median(col[feature.subset]) } )
            })
          }
          names(l2p) <- names(p)
        }
        else{
          if (denom[1] != "median"){
            l2p <- lapply( p, function(m) {
              if(entropy == FALSE) {
                apply( log2(m), 2, function(col) { col - mean(col[feature.subset]) } )
              } else if (entropy == TRUE) {
                apply( m, 2, function(col) { log2(col) - sum(-1*col*log2(col)) } )
              } 
            })
          } else if (denom[1] == "median"){
            l2p <- lapply( p, function(m) {
              apply( log2(m), 2, function(col) { col - median(col[feature.subset]) } )
            })
          }
        }
      }  else {
        message("the denominator is not recognized, use a different denominator")
      }
      
      # sanity check on data
      for ( i in 1:length(l2p) ) {
        if ( any( ! is.finite( l2p[[i]] ) ) ) stop("non-finite log-frequencies were unexpectedly computed")
      }
      if (verbose == TRUE) message("transformation complete")
    }
    
    return(new("aldex.clr",reads=reads,mc.samples=mc.samples,conds=conds,denom=feature.subset,verbose=verbose,useMC=useMC,dirichletData=p,analysisData=l2p, scaleSamps = scale_samples))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setMethod("getMonteCarloInstances", signature(.object="aldex.clr"), function(.object) .object@analysisData)

setMethod("getDirichletInstances", signature(.object="aldex.clr"), function(.object) .object@dirichletData)

setMethod("getSampleIDs", signature(.object="aldex.clr"), function(.object) names(.object@analysisData))

setMethod("getFeatureNames", signature(.object="aldex.clr"), function(.object) rownames(.object@analysisData[[1]]))

setMethod("getFeatures", signature(.object="aldex.clr"), function(.object) .object@analysisData[[1]][,1])

setMethod("numFeatures", signature(.object="aldex.clr"), function(.object) length(.object@analysisData[[1]][,1]))

setMethod("numMCInstances", signature(.object="aldex.clr"), function(.object) length(.object@analysisData[[1]][1,]))

setMethod("getReads", signature(.object="aldex.clr"), function(.object) .object@reads)

setMethod("numConditions", signature(.object="aldex.clr"), function(.object) length(names(.object@analysisData)))

setMethod("getDirichletReplicate", signature(.object="aldex.clr",i="numeric"), function(.object,i) .object@dirichletData[[i]])

setMethod("getDirichletSample", signature(.object="aldex.clr",i="numeric"), function(.object,i) sapply(.object@dirichletData, function(x) {x[,i]}))

setMethod("getMonteCarloReplicate", signature(.object="aldex.clr",i="numeric"), function(.object,i) .object@analysisData[[i]])

setMethod("getMonteCarloSample", signature(.object="aldex.clr",i="numeric"), function(.object,i) sapply(.object@analysisData, function(x) {x[,i]}))

setMethod("getDenom", signature(.object="aldex.clr"), function(.object) .object@denom)

setMethod("getConditions", signature(.object="aldex.clr"), function(.object) .object@conds)

setMethod("getScaleSamples", signature(.object="aldex.clr"), function(.object) .object@scaleSamps)

setMethod("aldex.clr", signature(reads="data.frame"), function(reads, conds, mc.samples=128, denom="all", verbose=FALSE, useMC=FALSE, gamma, entropy) aldex.clr.function(reads, conds, mc.samples, denom, verbose, useMC, summarizedExperiment=FALSE, gamma, entropy))

setMethod("aldex.clr", signature(reads="matrix"), function(reads, conds, mc.samples=128, denom="all", verbose=FALSE, useMC=FALSE, gamma, entropy) aldex.clr.function(as.data.frame(reads), conds, mc.samples, denom, verbose, useMC, summarizedExperiment=FALSE, gamma, entropy))

setMethod("aldex.clr", signature(reads="RangedSummarizedExperiment"), function(reads, conds, mc.samples=128, denom="all", verbose=FALSE, useMC=FALSE, gamma, entropy) aldex.clr.function(reads, conds, mc.samples, denom, verbose, useMC, summarizedExperiment=TRUE, gamma, entropy))
