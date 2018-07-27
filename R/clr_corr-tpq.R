#' Calculate correlation with a continuous variable
#' 
#' \code{aldex.corr} calculates the expected values for the correlation
#'  between each feature and a continuous variable, using data returned
#'  returned by \code{aldex.clr}.
#' 
#' @param clr An \code{ALDEx2} object. The output of \code{aldex.clr}.
#' @inheritParams aldex
#' @param ... Arguments passed to \code{cor.test}.
#' 
#' @return Returns a data.frame of the average
#'  coefficients and their p-values for each feature,
#'  with FDR appended as a \code{BH} column.
#' 
#' @author Thom Quinn
#' 
#' @seealso
#'  \code{\link{aldex}},
#'  \code{\link{aldex.clr}},
#'  \code{\link{aldex.ttest}},
#'  \code{\link{aldex.kw}},
#'  \code{\link{aldex.glm}},
#'  \code{\link{aldex.effect}},
#'  \code{\link{aldex.corr}},
#'  \code{\link{selex}}
#'  
#' @references Please use the citation given by
#'  \code{citation(package="ALDEx2")}.
#'  
#' @examples
#' data(selex)
#' #subset for efficiency
#' selex <- selex[1201:1600,]
#' x <- aldex.clr(selex, 1:14)
#' corr.test <- aldex.corr(x)
aldex.corr <- function(clr, verbose=FALSE, ...){
  
  # Use clr conditions slot instead of input
  conditions <- clr@conds
  
  lr2cor <- function(lr, conditions, ...){
    
    if(!is.numeric(conditions)){
      
      stop("Please define the aldex.clr object for a numeric 'conditions'.")
    }
    
    if(nrow(lr) != length(conditions)){
      
      stop("Incorrect length for 'conditions' argument.")
    }
    
    cors <- apply(lr, 2, function(x){
      cor.test(x, conditions, ...)
    })
    
    r <- sapply(cors, getElement, "statistic")
    p <- sapply(cors, getElement, "p.value")
    BH <- p.adjust(p, method = "BH")
    
    data.frame(r, p, BH,
               row.names = colnames(lr))
  }
  
  # Keep a running sum of lr2cor instances
  if(verbose) message("running tests for each MC instance:")
  mc <- ALDEx2::getMonteCarloInstances(clr)
  k <- ALDEx2::numMCInstances(clr)
  r <- 0
  for(i in 1:k){
    
    if(verbose) numTicks <- progress(i, k, numTicks)
    mci_lr <- t(sapply(mc, function(x) x[, i]))
    r <- r + lr2cor(mci_lr, conditions, ...)
  }
  
  r / k # return expected
}
