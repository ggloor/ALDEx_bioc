#' Compute an \code{aldex} Object
#'
#' @description
#' Welcome to the \code{ALDEx2} package!
#'
#' The \code{aldex} function is a wrapper that performs log-ratio transformation
#'  and statistical testing in a single line of code. Specifically, this function:
#'  (a) generates Monte Carlo samples of the Dirichlet distribution for each sample,
#'  (b) converts each instance using a log-ratio transform, then (c) returns test
#'  results for two sample (Welch's t, Wilcoxon) or multi-sample (glm, Kruskal-Wallace)
#'  tests. This function also estimates effect size for two sample analyses.
#'
#' @details
#' See "Examples" below for a description of the sample input.
#'
#' @param reads A non-negative, integer-only \code{data.frame} or \code{matrix}
#'  with unique names for all rows and columns. Rows should contain genes and
#'  columns should contain sequencing read counts (i.e., sample vectors).
#'  Rows with 0 reads in each sample are deleted prior to analysis.
#' @param conditions A character vector. A description of the data structure used
#'  for testing. Typically, a vector of group labels. For \code{aldex.glm}, use
#'  a \code{model.matrix}.
#' @param mc.samples An integer. The number of Monte Carlo samples to use when
#'  estimating the underlying distributions. Since we are estimating central tendencies,
#'  128 is usually sufficient.
#' @param denom A character string. Indicates which features to retain as the
#'  denominator for the Geometric Mean calculation. Using "iqlr" accounts for data
#'  with systematic variation and centers the features on the set features that have
#'  variance that is between the lower and upper quartile of variance. Using "zero"
#'  is a more extreme case where there are many non-zero features in one condition but
#'  many zeros in another. In this case the geometric mean of each group is calculated
#'  using the set of per-group non-zero features.
#' @param test A character string. Indicates which tests to perform. "t" runs
#'  Welch's t and Wilcoxon tests. "kw" runs Kruskal-Wallace and glm tests.
#'  "glm" runs a generalized linear model using a \code{model.matrix}.
#'  "corr" runs a correlation test using \code{cor.test}.
#' @param iterate A boolean. Toggles whether to iteratively perform a test. For example,
#'  this will use the results from an initial "t" routine to seed the reference
#'  (i.e., denominator of Geometric Mean calculation) for a second "t" routine.
#' @param effect A boolean. Toggles whether to calculate abundances and effect sizes.
#'  Applies to \code{test = "t"} and \code{test = "iterative"}.
#' @param include.sample.summary A boolean. Toggles whether to include median clr
#'  values for each sample. Applies to \code{effect = TRUE}.
#' @param paired.test A boolean. Toggles whether to do paired-sample tests.
#'  Applies to \code{effect = TRUE} and \code{test = "t"}.
#' @param verbose A boolean. Toggles whether to print diagnostic information while
#'  running. Useful for debugging errors on large datasets. Applies to
#'  \code{effect = TRUE}.
#' @param gamma A numeric. The standard deviation on the within sample variation.
#' @param bayesEst A boolean. Do we use the Bayes estimate for testing?
#' @param ... Arguments to embedded method (e.g., \code{glm} or \code{cor.test}).
#'
#' @return Returns a number of values that depends on the set of options.
#'  See the return values of aldex.ttest, aldex.kw, aldex.glm, and aldex.effect
#'  for explanations and examples.
#'
#' @author Greg Gloor, Andrew Fernandes, and Matt Links contributed to
#'  the original package. Thom Quinn added the "glm" test method, the
#'  "corr" test method, and the "iterate" procedure.
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
#' # The 'reads' data.frame should have row
#' # and column names that are unique, and
#' # looks like the following:
#' #
#' #              T1a T1b  T2  T3  N1  N2  Nx
#' #   Gene_00001   0   0   2   0   0   1   0
#' #   Gene_00002  20   8  12   5  19  26  14
#' #   Gene_00003   3   0   2   0   0   0   1
#' #   Gene_00004  75  84 241 149 271 257 188
#' #   Gene_00005  10  16   4   0   4  10  10
#' #   Gene_00006 129 126 451 223 243 149 209
#' #       ... many more rows ...
#'
#' data(selex)
#' selex <- selex[1201:1600,] # subset for efficiency
#' conds <- c(rep("NS", 7), rep("S", 7))
#' x <- aldex(selex, conds, mc.samples=2, denom="all",
#'            test="t", effect=TRUE, paired.test=FALSE)
aldex <- function(reads, conditions, mc.samples=128, test="t", effect=TRUE,
                  include.sample.summary=FALSE, verbose=FALSE,
                  denom="all", iterate=FALSE, gamma = NULL, bayesEst = TRUE, ...){

  if(missing(conditions)) stop("The 'conditions' argument is needed for this analysis.")

  # wrapper function for the entire set of
  message("aldex.clr: generating Monte-Carlo instances and clr values")
  x <- aldex.clr(reads=reads, conds=conditions, mc.samples=mc.samples,
                 denom=denom, verbose=verbose, useMC=FALSE, gamma = gamma)

  if(test == "t") {

    message("aldex.ttest: doing t-test")
    x.tt <- aldex.ttest(x, paired.test=FALSE, hist.plot=FALSE, verbose=verbose, bayesEst = bayesEst)

  }else if(test == "kw"){

    message("aldex.glm: doing Kruskal-Wallace and glm test (ANOVA-like)")
    x.tt <- aldex.kw(x)

  }else if(test == "glm"){

    message("aldex.glm: doing glm test based on a model matrix")
    x.tt <- aldex.glm(x, ...)

  }else if(test == "cor" | test == "corr"){

    message("aldex.corr: doing correlation with a continuous variable")
    x.tt <- aldex.corr(x, ...)

  }else{

    stop("argument 'test' not recognized")
  }

  if(iterate){

    message("iterate: seeding a second test")
    x.BHonly <- x.tt[,grepl("BH", colnames(x.tt)), drop = FALSE]
    nonDE.i <- as.logical(apply(x.BHonly > .05, 1, prod))
    if(sum(nonDE.i) == 0) stop("no non-DE references found")
    x.tt <- aldex(reads, conditions, mc.samples=mc.samples, test=test, effect=effect,
                  include.sample.summary=include.sample.summary, verbose=verbose,
                  denom=nonDE.i, iterate=FALSE, ...)
  }

  if(test == "t" && effect && !iterate){

    message("aldex.effect: calculating effect sizes")
    x.effect <- aldex.effect(x, include.sample.summary=include.sample.summary, verbose=verbose, paired.test=paired.test)
    z <- data.frame(x.effect, x.tt, check.names=F)

  }else{

    z <- data.frame(x.tt)
  }

  return(z)
}
