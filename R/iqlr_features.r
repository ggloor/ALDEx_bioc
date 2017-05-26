# iqlr_features.r
# Author: Jia Rong Wu
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
# The 'input' input can be "zero", "iqlr", "all", "" or numeric for advanced
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
            print("computing zero removal")
            features <- zero.features(reads,conds)
        } else if (denom == "iqlr") {
            print("computing iqlr centering")
            features <- iqlr.features(reads,conds)
        } else if (denom == "all" | denom == "") {
            print("computing center with all features")
            features <- all.features(reads,conds)
        } else {
            print(paste("denom: '", denom, "' unrecognized. Using all features.", sep=""))
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
            features <- custom.features(denom, conds)
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

    #####  remove all rows with reads less than the minimum set by minsum
    minsum <- 0

    ##### remove any row in which the sum of the row is 0
    z <- as.numeric(apply(reads, 1, sum))
    reads <- as.data.frame( reads[(which(z > minsum)),]  )

    ##### Adjust all reads with prior of 0.5
    reads <- reads + 0.5

    nr <- nrow( reads )
    rn <- rownames( reads )

    ##### Generate the CLR of the DATA and get variance of the CLR
    reads.clr <- t(apply(reads, 2, function(x){log2(x) - mean(log2(x))}))
    reads.var <- apply(reads.clr, 2, function(x){var(x)})
    reads.qtl <- quantile(unlist(reads.var))

    ##### Get the indicies of the "invariant set" features
    invariant.set <- which(
        (reads.var < (reads.qtl[4])) & (reads.var > (reads.qtl[2]))
    )

    for (i in 1:length(unique(conds)))
    {
        neg.indicies[[i]] <- invariant.set
    }

    return(neg.indicies)
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

##### Allows the user to submit a custom set of indicies for centering
##### Returns a list of "n" vectors where "n" is the number of conditions
custom.features <- function(denom, conds)
{

    custom.indicies <- vector("list", length(unique(conds)))

    # Set up the format that downstream clr_function.r is expecting
    for (i in 1: length(unique(conds)))
    {
        custom.indicies[[i]] <- denom
    }

    return(custom.indicies)
}


##### Returns a vector containing the indicies of every feature
all.features <- function(reads, conds)
{
    feature.indices <- as.numeric(seq(1, nrow(reads),1))
    return(feature.indices)
}

