
#' Calculate effect sizes and differences between all contrasts for the aldex.glm model matrix
#'
#' @name aldex.glm.effect
#' @title calculate effect sizes and differences between all contrasts for the aldex.glm model matrix
#' @description
#' data for this function is saved in a list with entries named by contrast
#' determines the median clr abundance of the feature in all samples and in groups
#' determines the median difference between the two groups
#' determines the median variation within each two group
#' determines the effect size, which is the median of the ratio of the between group difference and the larger of the variance within groups
#'
#' @usage aldex.glm.effect(clr, verbose = TRUE, include.sample.summary = FALSE,
#'   useMC=FALSE, CI=FALSE)
#' @param clr The data output of \code{aldex.clr}
#' @param verbose Print diagnostic information while running. Useful only for debugging if fails on large datasets
#' @param include.sample.summary include median clr values for each sample, defaults to FALSE
#' @param useMC use multicore by default (FALSE)
#' @param CI give effect 95\% confidence intervals, defaults to FALSE
#' @return A dataframe with the following information:
#' \describe{
#' \item{rab.all}{
#' a vector containing the median clr value for each feature
#' }
#' \item{rab.win.conditionA}{
#' a vector containing the median clr value for each feature in condition A
#' }
#' \item{rab.win.conditionB}{
#' a vector containing the median clr value for each feature in condition B
#' }
#' \item{diff.btw}{
#' a vector containing the per-feature median difference between condition A and B
#' }
#' \item{diff.win}{
#' a vector containing the per-feature maximum median difference between Dirichlet instances within conditions
#' }
#' \item{effect}{
#' a vector containing the per-feature effect size
#' }
#' \item{overlap}{
#' a vector containing the per-feature proportion of effect size that is 0 or less
#' }
#' }
#' @details An explicit example for two conditions is shown in the `Examples' below.
#' @references Please use the citation given by \code{citation(package="ALDEx")}.
#' @author Greg Gloor, Andrew Fernandes, Matt Links
#' @seealso \code{\link{aldex.clr}}, \code{\link{aldex.effect}}, \code{\link{aldex.ttest}}, \code{\link{aldex.glm}}, \code{\link{selex}}
#' @examples
#' # x is the output of the \code{x <- clr(data, mc.samples)} function
#' # conditions is a description of the data
#' # for the selex dataset, conditions <- c(rep("N", 7), rep("S", 7))
#' data(selex)
#' #subset for efficiency
#' selex <- selex[1201:1600,]
#' covariates <- data.frame("A" = sample(0:1, 14, replace = TRUE),
#' "B" = c(rep(0, 7), rep(1, 7)),
#' "Z" = sample(c(1,2,3), 14, replace=TRUE))
#' mm <- model.matrix(~ A + Z + B, covariates)
#' x <- aldex.clr(selex, mm, mc.samples=8, denom="all")
#' glm.effect <- aldex.glm.effect(x)
#'
#' @export 
aldex.glm.effect <- function(clr, verbose=TRUE, include.sample.summary=FALSE, useMC=FALSE, CI=FALSE){
# simple wrapper for clr_effect with more than two conditions
# returns the clr_effect output as a list when there are more than two conditions
# for use with the aldex.glm function

  if (is.vector(clr@conds)) {
    stop("only a single condition vector detected\n  use aldex.effect instead")
  } else if (is.matrix(clr@conds)) {
	effect.out <- list()
    names <- colnames(clr@conds)
    names <- names[-1]
    for(name in names){
      conds=clr@conds[,name]      
      conditions <- as.factor( conds )
      levels     <- levels( factor( conds) )

      if ( length( conds ) !=  numConditions(clr) ) stop("mismatch btw 'length(conditions)' and 'ncol(reads)'")
      if ( length( levels ) != 2 ) {
        warning("only two condition levels are currently supported\neffect not calculated for ",name)
        next
      }

      effect.out[[name]] <- aldex.effect(clr, glm.conds=conds, verbose=verbose, include.sample.summary=include.sample.summary, useMC=useMC, CI=CI)
    }
    return(effect.out)
  } else {
    stop("please check that an appropriate condition matrix was supplied to aldex.clr")
  }
}
