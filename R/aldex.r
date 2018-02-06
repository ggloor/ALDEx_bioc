#' Compute an \code{aldex} Object
#' 
#' @description
#' Welcome to the \code{ALDEx2} package!
#' 
#' The \code{aldex} function is a wrapper that performs log-ratio transformation
#'  and statistical testing in a single line of code. Specifically, this function:
#'  (a) generates Monte Carlo samples of the Dirichlet distribution for each sample,
#'  (b) converts each instance using a log-ratio transform, then (c) returns test
#'  results for two sample (Welch's t, Wilcoxon) or multi-sample (glm, Kruskal Wallace)
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
#'  for testing. Typically, a vector of group labels.
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
#' @param test A character string. Indicates which tests to perform. "t" calls
#'  Welch's t and Wilcoxon tests. "glm" calls Kruskal Wallace and glm tests.
#'  "iterative" uses the results from an initial "t" routine to seed the denominator
#'  (i.e., for the Geometric Mean calculation) of a second "t" routine.
#' @param effect A boolean. Toggles whether to calculate abundances and effect sizes.
#'  Applies to \code{test = "t"} and \code{test = "iterative"}.
#' @param include.sample.summary A boolean. Toggles whether to include median clr
#'  values for each sample. Applies to \code{effect = TRUE}.
#' @param verbose A boolean. Toggles whether to print diagnostic information while
#'  running. Useful for debugging errors on large datasets. Applies to
#'  \code{effect = TRUE}.
#' 
#' @return Returns a number of values that depends on the set of options.
#'  See the return values of aldex.ttest, aldex.glm, and aldex.effect for
#'  explanations and example.
#'  
#' @author Greg Gloor, Andrew Fernandes, and Matt Links contributed to
#'  the original package. Thom Quinn added the "iterative" test method.
#'  
#' @references Please use the citation given by \code{citation(package="ALDEx")}.
#' 
#' @seealso \code{\link{aldex.ttest}}, \code{\link{aldex.glm}},
#'  \code{\link{aldex.effect}}, \code{\link{aldex.corr}},
#'  \code{\link{selex}}
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
#'            test="t", effect=FALSE)
aldex <- function(reads, conditions, mc.samples=128, test="t",
                  effect=TRUE, include.sample.summary=FALSE, verbose=FALSE, denom="all"){
  
  if(missing(conditions)) stop("The 'conditions' argument is needed for this analysis.")
  
  # wrapper function for the entire set of
  print("aldex.clr: generating Monte-Carlo instances and clr values")
  x <- aldex.clr(reads=reads, conds=conditions, mc.samples=mc.samples,
                 denom=denom, verbose=verbose, useMC=FALSE)
  
  if(test == "iterative"){
    print("aldex.ttest: doing t-test")
    x.tt <- aldex.ttest(x, conditions, paired.test=FALSE)
    print("aldex.ttest: seeding a second t-test")
    nonDE.i <- which(rownames(reads) %in% rownames(x.tt[x.tt$wi.eBH > .05 | x.tt$we.eBH > .05, ]))
    if(length(nonDE.i) == 0) stop("no non-DE references found")
    x.tt <- aldex(reads, conditions, mc.samples=mc.samples, test="t",
                  effect=effect, include.sample.summary=include.sample.summary,
                  verbose=verbose, denom=nonDE.i)
  }else if(test == "t") {
    print("aldex.ttest: doing t-test")
    x.tt <- aldex.ttest(x, conditions, paired.test=FALSE)
  }else if(test == "glm"){
    print("aldex.glm: doing Kruskal Wallace and glm test")
    x.tt <- aldex.glm(x, conditions)
  }else{
    stop("argument 'test' not recognized")
  }
  
  if(effect == TRUE && test == "t"){
    print("aldex.effect: calculating effect sizes")
    x.effect <- aldex.effect(x, conditions,
                             include.sample.summary=include.sample.summary, verbose=verbose)
    z <- data.frame(x.effect, x.tt)
  }else{
    z <- data.frame(x.tt)
  }
  
  return(z)
}
