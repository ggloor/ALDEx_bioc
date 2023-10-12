#' Generate a differential scale matrix by group
#'
#' Takes as input the conditions vector, dispersion paramenter, starting
#' scale values and the number of random instances. The ratio between the 
#' scale values is key; setting mu = c(1,1.2) will have the same effect on
#' the analysis as a value of mu =(c(0.5,0.6). The function returns a matrix 
#' of scale values of the same dimension as the number of samples
#' and the number of mc.samples used by the aldex() or aldex.clr()
#' function.
#' 
#' @param gamma - the base gamma value for the sdlog parameter of rlnorm
#' @param mu - pair of values, or a vector of values one for each sample 
#' @param conditions - the conditions vector for the dataset
#' @param log - scale ratio in log2 (TRUE) or as simple ratio (FALSE)
#' @param mc.samples - the number of Monte-Carlo instances used by aldex()
#' 
#' @return returns a matrix of gamma values that are used as an estimate
#' of the scale for the aldex.clr() function.
#' This allows different scale and gamma values to be applied to each group
#' and can move the centre of mass of the data if required. The example dataset
#' has very extreme differences in scale. Most often these are likely in the
#' range of 10-15% of the base scale.
#'
#' @references Please use the citation given by \code{citation(package="ALDEx")}.
#' 
#' @author Greg Gloor, Michelle Pistner Nixon
#' 
#' @seealso \code{\link{aldex.clr}}, \code{\link{aldex}}
#' 
#' @examples
#'
#' # conditions is a vector describing the data
#' data(selex)
#' # subset for efficiency
#' conds <- c(rep("NS", 7), rep("S", 7))
#' mu.in <- c(1,50) # 50-fold difference in scale between groups
#' mu.vec <- aldex.makeScaleMatrix(1, mu.in, conds, log=TRUE, mc.samples=128)
#' 
#' @export
aldex.makeScaleMatrix <- function(gamma, mu, conditions, log=TRUE, mc.samples=128){
  ## new scale model
  # mu is the scale value and should be in units of fold difference in G
  # gamma is the dispersion paramenter
  if(!is.numeric(mu)) stop('the values in mu must be numeric')
  if(!is.vector(mu)) stop('the values in mu must be a vector')
  
  if(length(mu) == 2){
    mu1 = mu[1]
    mu2 = mu[2]    
	mu.vec <- gsub(levels(factor(conditions))[1], mu1, conditions)
    mu.vec <- as.numeric(gsub(levels(factor(conditions))[2], mu2, mu.vec))
  } else if (length(mu) == length(conditions)){
    mu.vec <- mu
  } else { 
    stop('provide either a groupwise or samplewise vector of scales')
  }
  gamma=gamma
  if(log == F){
    mu.vec <- log2(mu.vec)
  }

  # note: it is the log2 difference between mu1 and mu2 that is key here
  # eg; mu1=1, mu2=1.15 is equivalent to mu1=4, mu2=4.6
  # log2(1)=0, log2(1.15)~0.2; log2(4)=2, log2(4.6)~2.2
  
  return(t( sapply(mu.vec, FUN = function(mu) rlnorm(mc.samples, mu, gamma) ))/log(2))
  #return( 2^(-1 * log2(mat)) )
}
