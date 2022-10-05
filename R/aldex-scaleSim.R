#' Sensitivity analysis using scale simulation
#' 
#' Performs scale simulation over a range of values for gamma
#' Dirichlet samples are reused for computational convenience.
#' 
#' @param aldex_clr An `aldex.clr` object
#' @param gamma A vector of positive numeric components. Used as the standard deviation of the scale simulation model.
#' @param bayesEst A boolean. Do you want to use the Bayesian hypothesis testing method?
#' @inheritParams aldex
#' @return A list of results. Each element corresponds to a single result for a given value of gamma
#' @export
aldex.senAnalysis <- function(aldex_clr, gamma, test="t", effect=TRUE,
                              include.sample.summary=FALSE, verbose=FALSE,
                              iterate=FALSE, bayesEst = FALSE, ...){
  gamma <- sort(gamma)
  sen_results <- list()
  p <- getDirichletInstances(aldex_clr)
  mc.samples <- ncol(p[[1]])
  conds <- getConditions(aldex_clr)
  for(j in 1:length(gamma)){
    p <- getDirichletInstances(aldex_clr)
    l2p <- list()
    conds_mat <- matrix(conds, nrow = length(p))
    conds_mat <- apply(conds_mat, 2, FUN = function(vec) as.numeric(as.factor(vec)))
    conds_mat <- apply(conds_mat, 2, FUN = function(vec) vec - mean(vec))##Centering
    col_var <- gamma^2/apply(conds_mat, 2, var)
    scale_samples <- matrix(NA, length(p), mc.samples)
    
    for(i in 1:length(p)){
      geo_means <- log(apply(p[[i]],2,gm))
      noise <- sapply(col_var, FUN = function(sd){stats::rnorm(mc.samples, 0, sqrt(sd))})
      noise_mean <- rowMeans(noise)
      
      scale_samples[i,] <- geo_means + noise_mean
      l2p[[i]] <- sweep(log2(p[[i]]), 2,  scale_samples[i,], "-")
    }
    names(l2p) <- names(aldex_clr@dirichletData)
    x <-  new("aldex.clr",reads=clr@reads,mc.samples=clr@mc.samples,conds=clr@conds,
              denom=getDenom(aldex_clr),verbose=verbose,useMC=FALSE,dirichletData=getDirichletInstances(aldex_clr),analysisData=l2p, scaleSamps = gamma[j])
    if(test == "t") {
      
      message("aldex.ttest: doing t-test")
      x.tt <- aldex.ttest(x, paired.test=FALSE, hist.plot=FALSE, verbose=verbose, bayesEst = bayesEst)
      
    }else if(test == "kw"){
      
      message("aldex.glm: doing Kruskal-Wallace and glm test (ANOVA-like)")
      x.tt <- aldex.kw(x)
      
    }else if(test == "glm"){
      
      message("aldex.glm: doing glm test based on a model matrix")
      x.tt <- aldex.glm(x, bayesEst = bayesEst, ...)
      
    }else if(test == "cor" | test == "corr"){
      
      message("aldex.corr: doing correlation with a continuous variable")
      x.tt <- aldex.corr(x, ...)
      
    }else{
      
      stop("argument 'test' not recognized")
    }
    if(test == "t" && effect && !iterate){
      
      message("aldex.effect: calculating effect sizes")
      x.effect <- aldex.effect(x, include.sample.summary=include.sample.summary, verbose=verbose)
      z <- data.frame(x.effect, x.tt, check.names=F)
      
    }else{
      
      z <- data.frame(x.tt)
    }
    sen_results[[j]] <- z
  }
  names(sen_results) = paste0("gamma_", gamma)
  return(sen_results)
}

#' Create alpha diagram for scale simulation sensitivity result
#' 
#' @param sen_results A list return by aldex.senAnalysis()
#' @param test A character string. What test was used to calculate the results
#' @param thresh A numeric between 0 and 1. What threshold should be used for significance?
#' @param glmVar If `test = "glm"`, what variable do you want plotted?
#' @param bayesEst A boolean. Do you want to use the Bayesian hypothesis testing method?
#' @return A plot object
#' @export
plot_alpha <- function(sen_results, test = "t", thresh = 0.05, glmVar = NULL, bayesEst = FALSE){
  if(thresh < 0 | thresh > 1){
    stop("Please return a valid value for threshold")
  }
  
  gamma <- as.numeric(sub("gamma_", "", names(sen_results)))
  B <- matrix(NA, nrow = length(sen_results), ncol = dim(sen_results[[1]])[1])
  pvals <- matrix(NA, nrow = length(sen_results), ncol = dim(sen_results[[1]])[1])
  
  if(test == "t"){
    for(i in 1:length(sen_results)){
      B[i,] <- sen_results[[i]]$effect
      if(bayesEst){
        pvals[i,] <- sen_results[[i]]$p.val
      } else{
        pvals[i,] <- sen_results[[i]]$we.eBH
      }
    }
  } else if(test == "glm"){
    if(is.null(glmVar)){
      stop("Please supply what variable you want to plot!")
    }
    nameEffect = names(sen_results[[1]])[stringr::str_detect(names(sen_results[[1]]), paste0(glmVar, ".Estimate"))]
    namePval = names(sen_results[[1]])[stringr::str_detect(names(sen_results[[1]]), paste0(glmVar, ".Pr...t...BH"))]
    for(i in 1:length(sen_results)){
      B[i,] <- sen_results[[i]][,nameEffect]
      pvals[i,] <- sen_results[[i]][,namePval]
    }
  } else{
    stop("Test not supported by plot_alpha!")
  }
  
  
  P = as.data.frame(pvals)
  P$gamma = gamma
  P = P[,c(ncol(P), 1:(ncol(P)-1))]
  P = reshape(P, direction = "long", idvar = "gamma", varying = list(2:ncol(P)), times = names(P)[2:ncol(P)], v.names = "pval", timevar = "Sequence")

  P.toLabel = P[P$pval < thresh, ]
  P.toLabel = P.toLabel[order(P.toLabel$gamma, decreasing = TRUE),]
  P.toLabel = unique(P.toLabel$Sequence)
  P.toLabel = as.numeric(sub("V","",P.toLabel))
  

  B_graph <- as.data.frame(B)
  B_graph$gamma <- gamma
  B_graph <- B_graph[,c(ncol(B_graph), 1:(ncol(B_graph)-1))]
  B_graph <- reshape(B_graph, direction = "long", idvar = "gamma", varying = list(2:ncol(B_graph)), times = names(B_graph)[2:ncol(B_graph)], v.names = "Effect", timevar = "Sequence")
  B_graph <- merge(B_graph, P, by = c("gamma", "Sequence"))
  B_graph$Sequence <- sub("V", "", B_graph$Sequence)
  
  ##Switching the graph around
  B_graph$Effect = -B_graph$Effect
  B_graph$Sequence = as.numeric(B_graph$Sequence)
  B_graph = B_graph[order(B_graph$Sequence, B_graph$gamma),]
  p <- xyplot(Effect~gamma, data = B_graph, groups = Sequence, col = "grey", type = "l", xlab = "Gamma", ylab = "Effect")
  
  B_thresh <- B_graph[B_graph$pval <= thresh,]
  p2 <- direct.label(xyplot(Effect~gamma, data = B_thresh, groups = Sequence, type = "l", xlab = "Gamma", col= "black", ylab = "Effect"), "last.points")
  p3 <- xyplot(rep(0, length(unique(B_graph$gamma)))~unique(B_graph$gamma),  col = "red", type = "l", lty = "dashed", xlab = "Gamma", ylab = "Effect")
  p + as.layer(p2) + as.layer(p3)
  

}

gm <- function(x, na.rm = TRUE){
  exp(mean(log(x[x > 0]), na.rm=na.rm))
}
