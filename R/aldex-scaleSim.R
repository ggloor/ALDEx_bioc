#' Sensitivity analysis using scale simulation
#' 
#' Performs scale simulation over a range of values for lambda.
#' Dirichlet samples are reused for computational convenience.
#' 
#' @param aldex_clr An `aldex.clr` object
#' @param lambda A vector of positive numeric components. Used as the standard deviation of the scale simulation model.
#' @inheritParams aldex
#' @return A list of results. Each element corresponds to a single result for a given value of lambda
#' @export
aldex.senAnalysis <- function(aldex_clr, lambda, test="t", effect=TRUE,
                              include.sample.summary=FALSE, verbose=FALSE,
                              iterate=FALSE, ...){
  lambda <- sort(lambda)
  sen_results <- list()
  for(j in 1:length(lambda)){
    p <- getDirichletInstances(aldex_clr)
    l2p <- list()
    for(i in 1:length(p)){
      gm_sample <- log(apply(p[[i]],2,gm))
      scale_for_sample <- sapply(gm_sample, FUN = function(mu){stats::rlnorm(1, mu, lambda[j])})
      l2p[[i]] <- sweep(log2(p[[i]]), 2,  log2(scale_for_sample), "-")
    }
    names(l2p) <- names(aldex_clr@dirichletData)
    x <-  new("aldex.clr",reads=x@reads,mc.samples=x@mc.samples,conds=x@conds,
              denom=getDenom(aldex_clr),verbose=verbose,useMC=FALSE,dirichletData=getDirichletInstances(aldex_clr),analysisData=l2p, scaleSamps = lambda[j])
    if(test == "t") {
      
      message("aldex.ttest: doing t-test")
      x.tt <- aldex.ttest(x, paired.test=FALSE, hist.plot=FALSE, verbose=verbose)
      
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
    if(test == "t" && effect && !iterate){
      
      message("aldex.effect: calculating effect sizes")
      x.effect <- aldex.effect(x, include.sample.summary=include.sample.summary, verbose=verbose)
      z <- data.frame(x.effect, x.tt, check.names=F)
      
    }else{
      
      z <- data.frame(x.tt)
    }
    sen_results[[j]] <- z
  }
  names(sen_results) = paste0("lambda_", lambda)
  return(sen_results)
}

#' Create alpha diagram for scale simulation sensitivity result
#' 
#' @param sen_results A list return by aldex.senAnalysis()
#' @param test A character string. What test was used to calculate the results
#' @param thresh A numeric between 0 and 1. What threshold should be used for significance?
#' @param taxa_to_label A positive integer. How many taxa should be labeled in the plot?
#' @param glmVar If `test = "glm"`, what variable do you want plotted?
#' @return A ggplot2 object
#' @importFrom tidyr %>%
#' @importFrom stringr str_detect
#' @export
plot_alpha <- function(sen_results, test = "t", thresh = 0.05, taxa_to_label = 10, glmVar = NULL){
  if(thresh < 0 | thresh > 1){
    stop("Please return a valid value for threshold between zero and 1.")
  }
  
  lambda <- as.numeric(sub("lambda_", "", names(sen_results)))
  B <- matrix(NA, nrow = length(sen_results), ncol = dim(sen_results[[1]])[1])
  pvals <- matrix(NA, nrow = length(sen_results), ncol = dim(sen_results[[1]])[1])
  
  if(test == "t"){
    for(i in 1:length(sen_results)){
      B[i,] <- sen_results[[i]]$effect
      pvals[i,] <- sen_results[[i]]$we.eBH
    }
  } else if(test == "glm"){
    if(is.null(glmVar)){
      stop("Please supply what variable you want to plot!")
    }
    nameEffect = names(sen_results[[1]])[base::grepl(paste0(glmVar, ".Estimate"), names(sen_results[[1]]))]
    namePval = names(sen_results[[1]])[base::grepl(paste0(glmVar, ".Pr...t...BH"), names(sen_results[[1]]))]
    for(i in 1:length(sen_results)){
      B[i,] <- sen_results[[i]][,nameEffect]
      pvals[i,] <- sen_results[[i]][,namePval]
    }
  } else{
    stop("Test not supported by plot_alpha!")
  }
  
  if(taxa_to_label > dim(sen_results[[1]])[1]){
    message("Cannot label more taxa than exist. Reverting to all taxa in the data set.")
    taxa_to_label <- dim(sen_results[[1]])
  }
  
  P <- as.data.frame(pvals) 
  P$lambda <- lambda
  P <- P[,c(ncol(P), 1:(ncol(P)-1))]
  P <- reshape(P, idvar = "lambda",
               varying = list(2:ncol(P)),
               timevar = "Sequence",
               v.names = "pval", direction = "long")
  P$Sequence = paste0("V", P$Sequence)
  
  P.toLabel <- P[(P$pval < 0.1),]
  P.toLabel <- P.toLabel[order(-P.toLabel$lambda),]
  seq_to_label <- as.numeric(sub("V", "", unique(P.toLabel$Sequence)))
  
  taxa_to_label = as.vector(na.omit(seq_to_label[1:taxa_to_label]))
  
  B.graph <- as.data.frame(B)
  B.graph$lambda <- lambda
  B.graph <- B.graph[,c(ncol(B.graph), 1:(ncol(B.graph)-1))]
  B.graph <- reshape(B.graph, idvar = "lambda",
               varying = list(2:ncol(B.graph)),
               timevar = "Sequence",
               v.names = "Effect", direction = "long")
  B.graph$Sequence <- paste0("V", B.graph$Sequence)
  B.graph <- base::merge(B.graph, P, by = c("lambda", "Sequence"))
  B.graph$Sequence <- sub("V", "", B.graph$Sequence)
  B.graph$labl = ifelse(B.graph$Sequence %in% taxa_to_label, B.graph$Sequence, NA)
  
  ##Looping graph
  seq_max = unique(B.graph$Sequence)
  top = max(B.graph$Effect) + .5
  bottom = min(B.graph$Effect) - .5
  plot(1, type="n", xlab="Lambda", ylab="Effect Size", xlim=c(min(lambda), max(lambda)), ylim=c(bottom, top), panel.first = grid())
  for(i in seq_max){
    B.tmp = B.graph[B.graph$Sequence == i,]
    points(B.tmp$lambda, B.tmp$Effect, type = "l", col = "grey")
    B.tmp = B.tmp[B.tmp$pval <= thresh, ]
    points(B.tmp$lambda, B.tmp$Effect, type = "l", col = "black")
    
    B.tmp = B.tmp[nrow(B.tmp),]
    if(nrow(B.tmp) > 0){
      if(!is.na(B.tmp$labl)){
        text(x = B.tmp$lambda + runif(1,-.05,.05) , y = B.tmp$Effect + runif(1,-.25,.25), label = B.tmp$labl)
      }
    }
  }
  abline(h = 0, type = "l", col = "red", lty = "dashed")
}

gm <- function(x, na.rm = TRUE){
  exp(mean(log(x[x > 0]), na.rm=na.rm))
}
