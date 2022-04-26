#' Sensitivity analysis using scale simulation
#' 
#' Performs scale simulation over a range of values for gamma
#' Dirichlet samples are reused for computational convenience.
#' 
#' @param aldex_clr An `aldex.clr` object
#' @param gamma A vector of positive numeric components. Used as the standard deviation of the scale simulation model.
#' @inheritParams aldex
#' @return A list of results. Each element corresponds to a single result for a given value of gamma
#' @export
aldex.senAnalysis <- function(aldex_clr, gamma, test="t", effect=TRUE,
                              include.sample.summary=FALSE, verbose=FALSE,
                              iterate=FALSE, ...){
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
    x <-  new("aldex.clr",reads=x@reads,mc.samples=x@mc.samples,conds=x@conds,
              denom=getDenom(aldex_clr),verbose=verbose,useMC=FALSE,dirichletData=getDirichletInstances(aldex_clr),analysisData=l2p, scaleSamps = gamma[j])
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
  names(sen_results) = paste0("gamma_", gamma)
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
    stop("Please return a valid value for threshold")
  }
  
  gamma <- as.numeric(sub("gamma_", "", names(sen_results)))
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
    nameEffect = names(sen_results[[1]])[stringr::str_detect(names(sen_results[[1]]), paste0(glmVar, ".Estimate"))]
    namePval = names(sen_results[[1]])[stringr::str_detect(names(sen_results[[1]]), paste0(glmVar, ".Pr...t...BH"))]
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
  
  P = pvals %>% as.data.frame %>%
    as.data.frame() %>%
    dplyr::mutate("gamma" = gamma) %>%
    dplyr::select(gamma, everything()) %>%
    tidyr::pivot_longer(cols = !gamma, names_to = "Sequence", values_to = "pval")
  
  P.toLabel = P %>% dplyr::filter(pval < 0.1) %>%
    dplyr::arrange(desc(gamma)) %>%
    dplyr::select(Sequence) %>%
    unique() %>%
    dplyr::mutate(Sequence = as.numeric(sub("V","",Sequence)))
  
  taxa_to_label = P.toLabel$Sequence[1:taxa_to_label]
  
  B %>% 
    as.data.frame() %>%
    dplyr::mutate("gamma" = gamma) %>%
    dplyr::select(gamma, dplyr::everything()) %>%
    tidyr::pivot_longer(cols = !gamma, names_to = "Sequence", values_to = "Effect") %>%
    plyr::join(P, by = c("gamma", "Sequence")) %>%
    dplyr::mutate("Sequence" = sub("V", "", Sequence)) %>%
    dplyr::mutate("labl" = sub("V", "", Sequence)) %>%
    dplyr::mutate("labl" = ifelse(labl %in% taxa_to_label, labl, NA)) %>%
    ggplot(aes(x=gamma, y = Effect, group=Sequence)) +
    geom_line() +
    gghighlight((pval <= thresh), use_direct_label  = FALSE) +
    gghighlight(!is.na(labl), unhighlighted_params = list(colour = NULL)) +
    geom_hline(yintercept=0, color="red", linetype = "dashed") +
    theme_bw() +
    ylab("Effect Size") +
    scale_y_reverse() +
    xlab("Gamma") +
    theme(text = element_text(size=18))+
    theme(legend.position = "none") 
}

gm <- function(x, na.rm = TRUE){
  exp(mean(log(x[x > 0]), na.rm=na.rm))
}
