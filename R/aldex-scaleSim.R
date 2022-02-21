###Sensitivity analysis
###Reuse MonteCarlo samples

##WANT: Repeatedly call the aldex function buy able to reuse samples

##Helper function to calculate the geometric mean
gm <- function(x, na.rm = TRUE){
  exp(sum(log2(x[x > 0]), na.rm=na.rm) / length(x))
}