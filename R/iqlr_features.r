
#' Identify set of denominator features for log-ratio calculation
#'
#' @name aldex.set.mode
#' @title identify set of denominator features for log-ratio calculation
#' @description calculate the features that are to be used as the denominator for the Geometric Mean calculation in clr_function.R
#'
#' @usage aldex.set.mode(reads, conds, denom="all")
#' @param reads A data frame containing the samples and features per sample.
#' @param conds A vector describing which samples belong to what condition.
#' @param denom Character argument specifying which indicies to return.
#' 'all' returns all features in both conditons.
#' 'zero' returns the nonzero count features per condition.
#' 'iqlr' returns the features whose variance falls within the
#' inter-quantile range of the CLR-transformed data.
#' In cases of malformed or null queries, input defaults to 'all'.
#' Additionally, the input can be a numeric vector, which contains
#' a set of row indicies to center the data against. Only for advanced
#' users who can pre-determine the invariant set of features within
#' their data. Check that the same number of features are in the
#' input and output datasets.
#' @return Outputs a vector containing indices per condition, or a single vector in some cases.
#' @details An explicit example for two conditions is shown in the `Examples' below.
#' @references Please use the citation given by \code{citation(package="ALDEx")}.
#' @author Jia Rong Wu
#' @seealso \code{\link{aldex.clr}}, \code{\link{aldex.ttest}}, \code{\link{aldex.effect}}, \code{\link{selex}}
#' @examples
#' # x is the output of the \code{x <- clr(data, mc.samples)} function
#' # conditions is a description of the data
#' # for the selex dataset, conditions <- c(rep("N", 7), rep("S", 7))
#' # input can be "all", "iqlr", "zero" or numeric for advanced users
#' data(selex)
#' #subset for efficiency
#' selex <- selex[1201:1600,]
#' conds <- c(rep("NS", 7), rep("S", 7))
#' x <- aldex.clr(selex, conds, mc.samples=2, denom="all")
#' @export
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
        } else if (denom == "all" | denom == "" | denom == "median") {
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
    if(min(reads) == 0) reads <- reads + 0.5 

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
    if(length(invariant.set) <= 5) stop("Five or less members in intersecting feature set")

    return(as.vector(invariant.set))
}

##### Returns a list of vectors with the low variance, high
##### abundance features. Indistinguishable from user supplied

house.features <- function(reads, conds)
{
  neg.indicies <- vector("list", length(unique(conds)))

  invariant.set.list <- vector("list",
  length(unique(conds)))

  if(min(reads) == 0) {
    reads.n0 <- cmultRepl(t(reads), label=0, method="CZM",
     suppress.print=TRUE)
  } else {
    reads.n0 = reads
  }

  # clr transform
  reads.clr <- t(apply(reads.n0, 1, function(x){log2(x) - mean(log2(x))}))

  # per-condition offsets found
  for(i in 1:length(unique(conds))){
	these.rows <- which(conds == unique(conds)[i])

  # find the least variable
	reads.var <- apply(reads.clr[these.rows,],2, function(x){var(x)})
	var.set <- which(reads.var < quantile(unlist(reads.var))[2])

	# find the most relative abundant
	# top quartile in each sample
	rab.all <- apply(reads.clr[these.rows,], 2, function(x)
	    sum(x/length(these.rows)))
	abund.set <- which(rab.all > quantile(unlist(rab.all))[4])

	invariant.set.list[[i]] <-
	  intersect(var.set, abund.set)
  }

  # get the intersect of all conditions
  # successive operations on the list elements

  invariant.set <-
	Reduce(intersect, invariant.set.list)

  if(!length(invariant.set)) stop("No intersecting features are low variance and high relative abundance")
  if(length(invariant.set) <= 5) stop("Five or less members in intersecting feature set")

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

