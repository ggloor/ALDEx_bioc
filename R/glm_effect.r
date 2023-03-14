# simple wrapper for clr_effect with more than two conditions
# returns the clr_effect output as a list when there are more than two conditions
# for use with the aldex.glm function

aldex.glm.effect <- function(clr, verbose=TRUE, include.sample.summary=FALSE, useMC=FALSE, CI=FALSE){

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
