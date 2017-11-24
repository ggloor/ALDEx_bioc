aitchison.mean <- function( n, log=FALSE ) {

    # Input is a vector of non-negative integer counts.
    # Output is a probability vector of expected frequencies.
    # If log-frequencies are requested, the uninformative subspace is removed.

    n <- round( as.vector( n, mode="numeric" ) )
    if ( any( n < 0 ) ) stop("counts cannot be negative")

    a <- n + 0.5
    sa <- sum(a)

    log.p <- digamma(a) - digamma(sa)
    log.p <- log.p - mean(log.p)

    if ( log ) return(log.p)

    p <- exp( log.p - max(log.p) )
    p <- p / sum(p)
    return(p)
}


#<<BEGIN>>
#copied from mc2d R package
#licenced GPL>=2
#should be compatable with AGPL3
rdirichlet <- function (n, alpha)
#ISALIAS ddirichlet
#--------------------------------------------
{
  if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) can not be negative in rtriang")

  if(is.vector(alpha)) alpha <- t(alpha)
  l <- dim(alpha)[2]
  x <- matrix(rgamma(l * n, t(alpha)), ncol = l, byrow=TRUE)  # Gere le recycling
  return(x / rowSums(x))
}

# coerce data into dataframe
# reorder samples by condition

coerce.data <- function(reads, conds){
  # ncol df and length(c) must be equal
  if(ncol(reads) != length(conds)) stop("mismatch between number of samples and condition vector")

  # make an output list
  return( list(c.out <- conds[order(conds)], df.out <- data.frame(reads[,order(conds)])) )
}
