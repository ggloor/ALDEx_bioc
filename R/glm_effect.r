# simple wrapper for clr_effect with more than two conditions
# returns the clr_effect output as a list when there are more than two conditions
# for use with the aldex.glm function

aldex.glm.effect <- function(clr, verbose=TRUE, include.sample.summary=FALSE, useMC=FALSE, CI=FALSE){

  if (is.vector(clr@conds)) {
    stop("only a single condition vector detected\n    use aldex.effect instead")
  } else if (is.matrix(clr@conds)) {
    effect.out <- list()
    names <- colnames(clr@conds[,2:ncol(clr@conds)])
    for(i in 1:length(names)){
      conds=clr@conds[,i+1]
      conditions <- as.factor( conds )
      levels     <- levels( factor( conds) )

      if ( length( conds ) !=  numConditions(clr) ) stop("mismatch btw 'length(conditions)' and 'ncol(reads)'")
      if ( length( levels ) != 2 ) {
        warning("only two condition levels are currently supported\neffect not calculated for ",names[i])
        next
      }

      effect.out[[names[i]]] <- aldex.effect(clr, glm.conds=conds)
    }
    return(effect.out)
  } else {
    stop("please check that an appropriate condition matrix was supplied to aldex.clr")
  }
}
