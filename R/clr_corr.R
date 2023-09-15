#' Calculate correlation with a continuous variable
#'
#' \code{aldex.corr} calculates the expected values for the correlation
#'  between each feature and a continuous variable, using data returned
#'  returned by \code{aldex.clr} and a vector of the continuous variable.
#'  Returns results of Pearson, Spearman and Kendall tests.
#'
#' @param clr An \code{ALDEx2} object. The output of \code{aldex.clr}.
#' @param cont.var A continuous numeric vector
#'
#' @return Returns a data.frame of the average
#'  Pearson, Spearman and Kendall
#'  coefficients and their p-values for each feature,
#'  with FDR appended as a \code{BH} column.
#'
#' @author Arianne Albert, Greg Gloor, Thom Quinn
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
#' selex <- selex[1:50,]
#' conds <- c(rep("N", 7), rep("S",7))
#' cont.var <-  c(rep(1,7), rep(2,7))
#' x <- aldex.clr(selex, conds, mc.samples=16)
#' corr.test <- aldex.corr(x, cont.var)

#######################################
# a correlation function to correlate a continuous variable with the relative abundances

aldex.corr <- function(clr, cont.var){
  # cont.var is the continuous variable with which to run correlations

  # get dimensions, names, etc from the input data
  smpl.ids <- getSampleIDs(clr)
  feature.number <- numFeatures(clr)
  mc.instances <- numMCInstances(clr)
  feature.names <- getFeatureNames(clr)

  if ( length( cont.var ) !=  length(smpl.ids) )  stop("mismatch btw 'length(cont.var)' and 'length(smpl.ids)'")
  if ( is.numeric( cont.var ) != TRUE) stop("cont.var is not numeric")

  # set up the correlation results containers

  # Pearson correlation cor values, p-values and BH p-values
  pearson.matrix.cor <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
  pearson.matrix.p <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
  pearson.matrix.pBH <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)

  # Spearman rank correlation rho values, p-values and BH p-values
  spearman.matrix.rho <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
  spearman.matrix.p <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
  spearman.matrix.pBH <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)

  # kendall correlation rho values, p-values and BH p-values
  kendall.matrix.tau <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
  kendall.matrix.p <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)
  kendall.matrix.pBH <- matrix(data = NA, nrow = feature.number, ncol = mc.instances)

  #mc.i is the monte carlo instance
  for(mc.i in 1:mc.instances){
    # print(mc.i)

    #generate a matrix of each Monte-Carlo instance, columns are samples, rows are features
    t.input <- sapply(getMonteCarloInstances(clr), function(y){y[,mc.i]})

    # do Pearson correlations on each feature
    # make a list of correlation outputs
    x <- apply(t.input, 1, function(yy) {
      cor.test(as.numeric(yy), cont.var)
    })

    # p-values
    pearson.matrix.p[, mc.i] <- sapply(x, function(x) x[[3]])
    pearson.matrix.pBH[, mc.i] <- as.numeric(p.adjust(pearson.matrix.p[, mc.i], method = "BH"))
    # cor values
    pearson.matrix.cor[, mc.i] <- sapply(x, function(x) x[[4]])

    # Spearman correlations
    x <- apply(t.input, 1, function(yy) {
      cor.test(as.numeric(yy), cont.var, method = "spearman")
    })

    # p-values
    spearman.matrix.p[, mc.i] <- sapply(x, function(x) x[[3]])
    spearman.matrix.pBH[, mc.i] <- as.numeric(p.adjust(spearman.matrix.p[, mc.i], method = "BH"))
    # rho values
    spearman.matrix.rho[, mc.i] <- sapply(x, function(x) x[[4]])


    # kendall correlations
    x <- apply(t.input, 1, function(yy) {
      cor.test(as.numeric(yy), cont.var, method = "kendall")
    })

    # p-values
    kendall.matrix.p[, mc.i] <- sapply(x, function(x) x[[3]])
    kendall.matrix.pBH[, mc.i] <- as.numeric(p.adjust(kendall.matrix.p[, mc.i], method = "BH"))
    # rho values
    kendall.matrix.tau[, mc.i] <- sapply(x, function(x) x[[4]])

  }


  #get the Expected values of p, cor, rho, and lfdr
  pearson.ecor <- apply(pearson.matrix.cor, 1, mean)
  pearson.ep <- apply(pearson.matrix.p, 1, mean)
  pearson.eBH <- apply(pearson.matrix.pBH, 1, mean)

  spearman.erho <- apply(spearman.matrix.rho, 1, mean)
  spearman.ep <- apply(spearman.matrix.p, 1, mean)
  spearman.eBH <- apply(spearman.matrix.pBH, 1, mean)

  kendall.etau <- apply(kendall.matrix.tau, 1, mean)
  kendall.ep <- apply(kendall.matrix.p, 1, mean)
  kendall.eBH <- apply(kendall.matrix.pBH, 1, mean)

  z <- data.frame(pearson.ecor, pearson.ep, pearson.eBH, spearman.erho, spearman.ep,
    spearman.eBH, kendall.etau, kendall.ep, kendall.eBH)
  rownames(z) <- rownames(getMonteCarloInstances(clr)[[1]])
  return(z)

}
