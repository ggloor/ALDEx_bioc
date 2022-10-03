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

# modified from the R-Help mailing list
# https://stat.ethz.ch/pipermail/r-help/2000-December/009561.html
# modified to test for rational inputs first by gg Oct 09, 2018
# this should be clean of any GPL license

rdirichlet<-function(n,a)
  ## pick n random deviates from the Dirichlet function with shape parameters a
{
    if(length(n) > 1 || length(n) < 1 || n < 1) stop("n must be a single positive integer value")
    if(length(a) < 2) stop("a must be a vector of numeric value")
    n <- floor(n)
    l<-length(a);
    x<-matrix(rgamma(l*n,a),ncol=l,byrow=TRUE);
    return(x/rowSums(x) );
}
