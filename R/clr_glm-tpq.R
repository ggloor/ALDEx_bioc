#' Calculate glm test statistics using a \code{model.matrix}
#'
#' \code{aldex.glm} calculates the expected values for each coefficient of a
#'  glm model on the data returned by \code{aldex.clr}. This function
#'  requires the user to define a model with \code{model.matrix}.
#'
#' @param clr An \code{ALDEx2} object. The output of \code{aldex.clr}.
#' @param fdr.method A string ("BH" or "holm") denoting which method to use to adjust p-values. Default is "holm"

#' @inheritParams aldex
#' @param ... Arguments passed to \code{glm}.
#' @return Returns a data.frame of the average
#'  coefficients and their p-values for each feature,
#'  with FDR appended as a \code{holm} column.
#'
#' @author Thom Quinn, Michelle Pistner
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
#' covariates <- data.frame("A" = sample(0:1, 14, replace = TRUE),
#'                          "B" = c(rep(0, 7), rep(1, 7)))
#' mm <- model.matrix(~ A + B, covariates)
#' x <- aldex.clr(selex, mm, mc.samples=4, denom="all")
#' glm.test <- aldex.glm(x)
#' glm.eff <- aldex.glm.effect(x)
#' aldex.glm.plot(glm.test, eff=glm.eff, contrast='B', type='MW', post.hoc='holm')
#'
aldex.glm <- function(clr, verbose=FALSE, fdr.method = "holm", ...){
  
  if(!(fdr.method %in% c("holm", "BH", "fdr"))){
    stop("Method to adjust p-values not supported.")
  }

  # Use clr conditions slot instead of input
  conditions <- clr@conds

  # Keep a running sum of lr2glm instances
  # verbose was throwing an error 'the condition has length > 1'
  if(verbose[1] == TRUE) message("running tests for each MC instance:")
  mc <- ALDEx2::getMonteCarloInstances(clr)
  k <- ALDEx2::numMCInstances(clr)
  r <- 0
  r.p.lower <- 0
  r.p.upper <- 0
  r.bh.lower <- 0
  r.bh.upper <- 0
  for(i in 1:k){

    if(verbose[1] == TRUE ){ numTicks <- progress(i, k, numTicks) }
    mci_lr <- t(sapply(mc, function(x) x[, i]))
    mod <- lr2glm(mci_lr, conditions, fdr.method = fdr.method, ...)
    r <- r + mod$df
    r.p.lower <- r.p.lower + mod$df.p.lower
    r.p.upper <- r.p.upper + mod$df.p.greater
    r.bh.lower <- r.bh.lower + mod$df.bh.lower
    r.bh.upper <- r.bh.upper + mod$df.bh.greater
  }
  # simplify names for plotting
  colnames(r) <- gsub("\\(Intercept)", 'Intercept:', colnames(r))
  colnames(r) <- gsub('model.', '', colnames(r))
  colnames(r) <- gsub(' Estimate', ':Est', colnames(r))
  colnames(r) <- gsub(' t value', ':t.val', colnames(r))
  colnames(r) <- gsub(' Std. Error', ':SE', colnames(r))
  colnames(r) <- gsub(" Pr\\(.+\\)", ':pval', colnames(r))

  
  ## filtering through the both the pvalus and adjusted pvals
  r.p <- r.p.lower / k
  for(j in 1:ncol(r.p)){
    tmp <- cbind(r.p.upper[,j]/k, r.p.lower[,j]/k)
    r.p[,j] <- apply(tmp, 1, min)
  }
  
  
  r.p.adj <- r.bh.lower / k ## creating a matrix to store results
  ## runs through the variables and selects the appropriate two-sided test
  for(j in 1:ncol(r.p.adj)){
      tmp <- cbind(r.bh.upper[,j]/k, r.bh.lower[,j]/k)
      r.p.adj[,j] <- apply(tmp, 1, min)
    }
  ##cleaning up column names for p-value matrices
  colnames(r.p) <- gsub(" Pr\\(.+\\)", ':pval', colnames(r.p))
  colnames(r.p) <- gsub("model.", '', colnames(r.p))
  colnames(r.p.adj) <- sub(" Pr\\(.+\\)", ':pval', colnames(r.p.adj))
  colnames(r.p.adj) <- sub("model.", '', colnames(r.p.adj))
  
  
  cbind(r/k, r.p, r.p.adj) # return expected
}

# declaring this once provides a 40% speedup
# Extract coefficients and p-values
extract <- function(model){
  x <- coef(summary(model))
  coefs <- lapply(1:nrow(x), function(i){
	y <- x[i,,drop=FALSE]
	colnames(y) <- paste(rownames(y), colnames(y))
	y})
  do.call("cbind", coefs)
}

lr2glm <- function(lr, conditions, fdr.method = "holm", ...){
  
  if( !is(conditions, "matrix") &
      !("assign" %in% names(attributes(conditions)))){
    
    stop("Please define the aldex.clr object for a model.matrix 'conditions'.")
  }
  
  if(nrow(lr) != nrow(conditions)){
    
    stop("Input data and 'model.matrix' should have same number of rows.")
  }
  
  # Build the glm models
  model. <- conditions
  glms <- apply(lr, 2, function(x){
    glm(x ~ model., ...)
  })
  dof <- glm(lr[,1]~model.)$df.residual
  # Combine to make data.frame
  extracts <- lapply(glms, extract)
  df <- do.call("rbind", extracts)
  rownames(df) <- colnames(lr)
  df <- as.data.frame(df)
  
  # Adjusting to one-sided p-values
  colsToAdjust <- which(grepl("Pr\\(>",colnames(df)))
  for(j in colsToAdjust){
    df[,j] <- pt(df[,(j-1)], df = dof,lower.tail = FALSE)
  }
  df.p.greater <- 2*df[,colsToAdjust]
  df.p.lower <- 2*(1-df[,colsToAdjust])
  
  ##Making sure no p-values are over 1
  df.p.greater <- apply(df.p.greater, c(1,2), FUN = function(x) min(x,1))
  df.p.lower <- apply(df.p.lower, c(1,2), FUN = function(x) min(x,1))
  
  
  # Create new data.frame for p-values and FDR
  pvals <- colnames(df)[grepl("Pr\\(>", colnames(df))]
  df.bh.greater <- df[,pvals]
  colnames(df.bh.greater) <- paste0(colnames(df.bh.greater), ".padj")
  
  df.bh.lower <- 1-df[,pvals]
  colnames(df.bh.lower) <- paste0(colnames(df.bh.lower), ".padj")
  for(j in 1:ncol(df.bh.greater)){
    df.bh.greater[,j] <- p.adjust(2*df.bh.greater[,j], method=fdr.method)
    df.bh.lower[,j] <- p.adjust(2*df.bh.lower[,j], method=fdr.method)
  }
  
  # Merge results with FDR
  list(df = df[,-colsToAdjust], df.p.greater = df.p.greater, df.p.lower = df.p.lower, df.bh.greater = df.bh.greater, df.bh.lower = df.bh.lower)
}
