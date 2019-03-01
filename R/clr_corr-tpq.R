#' Calculate correlation with a continuous variable
#'
#' \code{aldex.corr} calculates the expected values for the correlation
#'  between each feature and a continuous variable, using data returned
#'  returned by \code{aldex.clr} and a vector of the continuous variable.
#'  By default uses pearson but method="kendall" or "spearman" can be
#'  passed to the cor.test function.
#'
#' @param clr An \code{ALDEx2} object. The output of \code{aldex.clr}.
#' @param cont.var A continuous numeric vector
#' @inheritParams aldex
#' @param ... Arguments passed to \code{cor.test}.
#'
#' @return Returns a data.frame of the average
#'  coefficients and their p-values for each feature,
#'  with FDR appended as a \code{BH} column.
#'
#' @author Thom Quinn, Greg Gloor
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
#' conds <- c(rep("N", 7), rep("S",7))
#' cont.var <-  1:14
#' x <- aldex.clr(selex, conds)
#' corr.test <- aldex.corr(x, cont.var)
aldex.corr <- function(clr, cont.var, verbose=FALSE, ...){

  # Use var instead of clr conditions slot
  conditions <- cont.var

  lr2cor <- function(lr, conditions, ...){

    if(!is.numeric(conditions)){

      stop("Please define the 'cont.var' argument as a numeric.")
    }

    if(nrow(lr) != length(conditions)){

      stop("Incorrect length for 'cont.var' argument.")
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
