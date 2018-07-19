#' Calculate the Kruskal-Wallis test and glm ANOVA statistics
#' 
#' \code{aldex.kw} calculates the expected values of the Kruskal-Wallis
#'  test and a glm ANOVA on the data returned by \code{aldex.clr}.
#' 
#' @param clr An \code{ALDEx2} object. The output of \code{aldex.clr}.
#' @inheritParams aldex
#' @param useMC Toggles whether to use multi-core.
#' 
#' @return Returns a \code{data.frame} with the following information:
#' \item{kw.ep}{ a vector containing the expected p-value of the Kruskal-Wallis test
#'  for each feature }
#' \item{kw.eBH}{ a vector containing the corresponding expected value of the
#'  Benjamini-Hochberg corrected p-value for each feature }
#' \item{glm.ep}{ a vector containing the expected p-value of the glm ANOVA
#'  for each feature }
#' \item{glm.eBH}{ a vector containing the corresponding expected value of the
#'  Benjamini-Hochberg corrected p-value for each feature }
#' 
#' @author Arianne Albert
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
#' conds <- c(rep("A", 4), rep("B", 3), rep("C", 7))
#' x <- aldex.clr(selex, conds, mc.samples=1, denom="all")
#' kw.test <- aldex.kw(x, conds)
aldex.kw <- function(clr, conditions, useMC=FALSE){
  
  conditions <- clr@conds
  
  # make sure that the multicore package is in scope and return if available
  is.multicore = FALSE
  
  if ("BiocParallel" %in% rownames(installed.packages()) & useMC){
    print("multicore environment is OK -- using the BiocParallel package")
    #require(BiocParallel)
    is.multicore = TRUE
  }
  else {
    print("operating in serial mode")
  }
  
  # get dimensions, names, etc from the input data
  smpl.ids <- getSampleIDs(clr)
  feature.names <- getFeatureNames(clr)
  feature.number <- numFeatures(clr)
  mc.instances <- numMCInstances(clr)
  
  conditions <- as.factor( conditions )
  levels     <- levels( conditions )
  
  if ( length( conditions ) !=  numConditions(clr) ){
    stop(paste("mismatch btw 'length(conditions)' and 'length(names(clr))'. len(condtitions):",
               length(conditions),"len(names(clr)):",numConditions(clr)))}
  
  # generate the comparison sets from the condition levels
  levels <- vector( "list", length( levels ) )
  names( levels ) <- levels( conditions )
  sets <- names(levels)
  
  # set up the glm results containers
  glm.matrix.p <- as.data.frame(matrix(1, nrow = feature.number, ncol = mc.instances))
  glm.matrix.pBH <- glm.matrix.p # duplicate result container
  kw.p.matrix <- glm.matrix.p # duplicate result container
  kw.pBH.matrix <- glm.matrix.p # duplicate result container
  
  # mc.i is the i-th Monte-Carlo instance
  print("running tests for each MC instance:")
  mc.all <- getMonteCarloInstances(clr)
  for(mc.i in 1:mc.instances){
    
    numTicks <- progress(mc.i, mc.instances, numTicks)
    
    # generate a matrix of i-th Monte-Carlo instance, columns are samples, rows are features
    t.input <- sapply(mc.all, function(y){y[, mc.i]})
    
    # do glms on each feature and make a list of glm outputs
    x <-
      apply(t.input, 1, function(yy) {
        glm(as.numeric(yy) ~ factor(conditions))
      })
    
    # calculate p-values for generalized linear model
    if(is.multicore == TRUE){
      pps <- bplapply(x, drop1, test = "Chis")
    }else{
      pps <- lapply(x, drop1, test = "Chis")
    }
    glm.matrix.p[, mc.i] <- sapply(pps, function(x){x[[5]][2]})
    glm.matrix.pBH[, mc.i] <- as.numeric(p.adjust(glm.matrix.p[, mc.i], method = "BH"))
    
    # calculate p-values for Kruskal Wallis test
    kw.p.matrix[, mc.i] <-
      apply(t.input, 1, function(yy){
        kruskal.test(yy, g = factor(conditions))[[3]]
      })
    kw.pBH.matrix[, mc.i] <- as.numeric(p.adjust(kw.p.matrix[, mc.i], method = "BH"))
    
  }
  
  # get the Expected values of p, q and lfdr
  glm.ep <- rowMeans(glm.matrix.p) # rowMeans is faster than apply()!!
  glm.eBH <- rowMeans(glm.matrix.pBH)
  kw.ep <- rowMeans(kw.p.matrix)
  kw.eBH <- rowMeans(kw.pBH.matrix)
  
  z <- data.frame(kw.ep, kw.eBH, glm.ep, glm.eBH)
  rownames(z) <- getFeatureNames(clr)
  return(z)
}
