# iqlr_features.r
# Author: Jia Rong Wu, Greg Gloor
#
# USAGE: source("iqlr_features.R")
#
# DESCRIPTION:
# Functions for declaring a set of features to be used in clr_function.r
# The output returned is a list of indicies representing the features to be used
# in either the center log-ratio or the inter quantile log-ratio transformation.
#
# INPUT:
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
# The 'conds' input can be a simple vector
# i.e.
# conds <- c(rep("A", 10), rep("B", 10))
#
# The 'input' input can be "zero", "iqlr", "all", "lvha", "" or numeric for advanced
# users.
#
# 'input' defaults to "all" if either no parameters or incorrect parameters
# given
# a numeric 'input' variable is only for users who understand the set of
# invariant features in their dataset

##### Determines the mode to be used in ALDEx2
##### Default is all features
aldex.set.mode <- function(reads, conds, denom="all")
{

    if (is.character(denom))
    {
        if (denom == "zero") {
            message("computing zero removal")
            features <- zero.features(reads,conds)
        } else if (denom == "iqlr") {
            message("computing iqlr centering")
            features <- iqlr.features(reads,conds)
        } else if (denom == "all" | denom == "") {
            message("computing center with all features")
            features <- all.features(reads,conds)
        } else if (denom == "lvha" ) {
            message("computing center with housekeeping features")
            features <- house.features(reads,conds)
        } else {
            message(paste("denom: '", denom, "' unrecognized. Using all features.", sep=""))
            features <- all.features(reads,conds)
        }
        return(features)
    } else {

        # Prevent input of null features for centering
        if (length(denom) == 0)
        {
            stop("Number of features can not be 0. Please check the input vector.")
        }
        else
        {
            # Return user specified indicies for centering
            features <- denom
            return(features)
        }
    }
}

##### Returns a list of "n" vectors where n is the number of conditions
##### Each vector represents the indicies of the IQLR features for the
##### conditions
iqlr.features <- function(reads, conds)
{
    neg.indicies <- vector("list", length(unique(conds)))

    ##### Adjust all reads with prior of 0.5
    reads <- reads + 0.5

    nr <- nrow( reads )
    rn <- rownames( reads )

    ##### Generate the CLR of the DATA and get variance of the CLR
    reads.clr <- t(apply(reads, 2, function(x){log2(x) - mean(log2(x))}))
 # per condition typical variance
	for(i in 1:length(unique(conds))){
	  these.rows <- which(conds == unique(conds)[i])

	  reads.var <- apply(reads.clr[these.rows,], 2, function(x){var(x)})
	  reads.qtl <- quantile(unlist(reads.var))

	  ##### Get the indicies of the "invariant set" features
	  this.set <- which(
		  (reads.var < (reads.qtl[4])) & (reads.var > (reads.qtl[2]))
	  )
	  neg.indicies[[i]] <- this.set

    }
 # total dataset typical variance
    reads.var <- apply(reads.clr, 2, function(x){var(x)})
	reads.qtl <- quantile(unlist(reads.var))
	  this.set <- which(
		  (reads.var < (reads.qtl[4])) & (reads.var > (reads.qtl[2]))
	  )
	neg.indicies[[length(unique(conds)) + 1]] <- this.set


    invariant.set <- Reduce(intersect, neg.indicies)

    if(!length(invariant.set)) stop("No intersecting features have typical variance")
    return(as.vector(invariant.set))
}

##### Returns a list of vectors with the low variance, high
##### abundance features. Indistinguishable from user supplied

house.features <- function(reads, conds)
{
  neg.indicies <- vector("list", length(unique(conds)))

  invariant.set.list <- vector("list",
  length(unique(conds)))

  # clr transform
  reads.clr <- t(apply(reads + 0.5, 2,
  function(x){log2(x) - mean(log2(x))}))

  # per-condition offsets found
  for(i in 1:length(unique(conds))){
	these.rows <- which(conds ==
	  unique(conds)[i])

  # find the least variable
	reads.var <- apply(reads.clr[these.rows,],2, function(x){var(x)})
	var.set <- which(reads.var
	  < quantile(unlist(reads.var))[2])

	# find the most relative abundant
	# top quartile in each sample
	quantile.sample <- apply(reads.clr, 1, quantile)

	abund.set <- apply(reads.clr, 2, function(x)
      sum(x > quantile.sample[4,]) == length(quantile.sample[4,]))

	invariant.set.list[[i]] <-
	  intersect(var.set, which(abund.set == TRUE))
  }
  # find the least variable
	reads.var <- apply(reads.clr[],2, function(x){var(x)})
	var.set <- which(reads.var
	  < quantile(unlist(reads.var))[2])

	# find the most relative abundant
	# top quartile in each sample
	quantile.sample <- apply(reads.clr, 1, quantile)

	abund.set <- apply(reads.clr, 2, function(x)
      sum(x > quantile.sample[4,]) == length(quantile.sample[4,]))

	invariant.set.list[[length(unique(conds))+1]] <-
	  intersect(var.set, which(abund.set == TRUE))
  # get the intersect of all conditions
  # successive operations on the list elements

  invariant.set <-
	Reduce(intersect, invariant.set.list)

  if(!length(invariant.set)) stop("No intersecting features are low variance and high relative abundance")
  return(invariant.set)
}


##### Returns a list of "n" vectors where n is the number of conditions
##### Each vector represents the indicies of NONZERO features for that condition
zero.features <- function(reads, conds)
{
    sample.indices <- as.numeric(seq(1, length(conds),1))
    feature.indices <- as.numeric(seq(1, nrow(reads),1))

    indicies <- vector("list", length(unique(conds)))
    neg.indicies <- vector("list", length(unique(conds)))
    condition.list <- vector("list", length(unique(conds)))

    # Iterate through conditions and check indicies of features that are zero in
    # one condition but not the other
    for (i in 1: length(unique(conds)))
    {
        condition.list[[i]] <- which(conds == unique(conds)[i]) # Condition list
        sub.set <- as.numeric(condition.list[[i]])
        different.conds <- setdiff(sample.indices, sub.set)
        for (j in 1:nrow(reads))
        {
            if ( (sum(reads[j,sub.set]) == 0)
                & (sum(reads[j,different.conds]) > 0) )
            {
                indicies[[i]] <- c(indicies[[i]],j)            # Save the indexes
            }
        }

        # Set the negation of the indicies to retain the nonzero features
        neg.indicies[[i]] <- setdiff(feature.indices, indicies[[i]])
    }

    return(neg.indicies)
}


##### Returns a vector containing the indicies of every feature
all.features <- function(reads, conds)
{
    feature.indices <- as.numeric(seq(1, nrow(reads),1))
    return(feature.indices)
}

