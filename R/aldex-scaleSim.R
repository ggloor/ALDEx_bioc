#' Sensitivity analysis using scale simulation
#' 
#' Performs scale simulation over a range of values for gamma.
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
  
  ##extracting results from the aldex object
  p <- getDirichletInstances(aldex_clr)
  conds <- aldex_clr@conds
  mc.samples <- aldex_clr@mc.samples
  
  ## repeating over gamma
  for(j in 1:length(gamma)){
    l2p <- list()
    
    ##adding scale
    scale_samples <- default.scale.model(gamma[j], conds, p, mc.samples)
    for(i in 1:length(p)){
      l2p[[i]] <- sweep(log2(p[[i]]), 2,  scale_samples[i,], "-")
    }
    
    names(l2p) <- names(aldex_clr@dirichletData)
    x <-  new("aldex.clr",reads=aldex_clr@reads,mc.samples=aldex_clr@mc.samples,conds=aldex_clr@conds,
              denom=getDenom(aldex_clr),verbose=verbose,useMC=FALSE,dirichletData=getDirichletInstances(aldex_clr),analysisData=l2p, scaleSamps = scale_samples)
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

#' Create gamma diagram for scale simulation sensitivity result
#' 
#' @param sen_results A list return by aldex.senAnalysis()
#' @param test A character string. What test was used to calculate the results
#' @param thresh A numeric between 0 and 1. What threshold should be used for significance?
#' @param taxa_to_label A positive integer. How many taxa should be labeled in the plot?
#' @param glmVar If `test = "glm"`, what variable do you want plotted?
#' @param blackWhite boolean. If TRUE, returns the plot in black and white.
#' @param cex Default == 1. Controls the size of the axis and text labels in the plots.
#' @return A plot object
#' 
#' @export
plotGamma <- function(sen_results, test = "t", thresh = 0.05, taxa_to_label = 10, glmVar = NULL, blackWhite = FALSE, cex = 1){
  Sequence <- NULL ##Fixing an R CMD check note
  if(thresh < 0 | thresh > 1){
    stop("Please return a valid value for threshold between zero and 1.")
  }
  
  gamma <- as.numeric(sub("gamma_", "", names(sen_results)))
  B <- matrix(NA, nrow = length(sen_results), ncol = dim(sen_results[[1]])[1])
  pvals <- matrix(NA, nrow = length(sen_results), ncol = dim(sen_results[[1]])[1])
  
  if(test == "t"){
    for(i in 1:length(sen_results)){
      pvals[i,] <- sen_results[[i]]$we.eBH
      B[i,] <- sen_results[[i]]$effect
    }
  } else if(test == "glm"){
    if(is.null(glmVar)){
      stop("Please supply what variable you want to plot!")
    }
    nameEffect = names(sen_results[[1]])[base::grepl(paste0(glmVar, ".Est"), names(sen_results[[1]]))]
    namePval = names(sen_results[[1]])[base::grepl(paste0(glmVar, ".pval.holm"), names(sen_results[[1]]))]
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
  P$gamma <- gamma
  P <- P[,c(ncol(P), 1:(ncol(P)-1))]
  P <- reshape(P, idvar = "gamma",
               varying = list(2:ncol(P)),
               timevar = "Sequence",
               v.names = "pval", direction = "long")
  P$Sequence = paste0("V", P$Sequence)
  
  P.toLabel <- P[(P$pval < 0.1),]
  P.toLabel <- P.toLabel[order(-P.toLabel$gamma),]
  seq_to_label <- as.numeric(sub("V", "", unique(P.toLabel$Sequence)))
  
  taxa_to_label = as.vector(na.omit(seq_to_label[1:taxa_to_label]))
  
  B.graph <- as.data.frame(B)
  B.graph$gamma <- gamma
  B.graph <- B.graph[,c(ncol(B.graph), 1:(ncol(B.graph)-1))]
  B.graph <- reshape(B.graph, idvar = "gamma",
                     varying = list(2:ncol(B.graph)),
                     timevar = "Sequence",
                     v.names = "Effect", direction = "long")
  B.graph$Sequence <- paste0("V", B.graph$Sequence)
  B.graph <- base::merge(B.graph, P, by = c("gamma", "Sequence"))
  B.graph$Sequence <- sub("V", "", B.graph$Sequence)
  B.graph$labl = ifelse(B.graph$Sequence %in% taxa_to_label, B.graph$Sequence, NA)
  
  B.graphFilt <- B.graph[B.graph$pval <= thresh,]
  
  B.graphFilt$Sequence <- as.factor(B.graphFilt$Sequence)

  if(blackWhite){
    base.plot <- xyplot(Effect~ gamma,
                        groups = Sequence,
                        data = B.graph,
                        col = "grey",
                        type = "l",
                        ylab=list(label = "Effect Size", cex = cex),
                        xlab=list(label = "Gamma", cex = cex),
                        abline = c(h=0),
                        par.settings = standard.theme(color = FALSE))
    
    overlay.plot <- xyplot(Effect~gamma,
                           groups = Sequence,
                           data = B.graphFilt,
                           color = "black",
                           type = "l",
                           par.settings = standard.theme(color = FALSE))
    
    B.graphLabl <- B.graph[!is.na(B.graph$labl), ]
    B.graphLabl <- B.graphLabl[B.graphLabl$pval <= thresh,]
    
    overlay.labels <- direct.label(xyplot(Effect~gamma,
                                          groups = Sequence,
                                          data = B.graphLabl,
                                          color = "black",
                                          type = "l",
                                          par.settings = standard.theme(color = FALSE)), list("last.bumpup", cex = cex))
    
    p2 <- base.plot + as.layer(overlay.plot) + as.layer(overlay.labels)
  } else{
    base.plot <- xyplot(Effect~ gamma,
                        groups = Sequence,
                        data = B.graph,
                        col = "grey",
                        type = "l",
                        ylab=list(label = "Effect Size", cex = cex),
                        xlab=list(label = "Gamma", cex = cex),
                        abline = c(h=0))
    
    overlay.plot <- xyplot(Effect~gamma,
                           groups = Sequence,
                           data = B.graphFilt,
                           color = "black",
                           type = "l")
    
    B.graphLabl <- B.graph[!is.na(B.graph$labl), ]
    B.graphLabl <- B.graphLabl[B.graphLabl$pval <= thresh,]
    
    overlay.labels <- direct.label(xyplot(Effect~gamma,
                                          groups = Sequence,
                                          data = B.graphLabl,
                                          color = "black",
                                          type = "l"), list("last.bumpup", cex = cex))
    
    p2 <- base.plot + as.layer(overlay.plot) + as.layer(overlay.labels)
  }
  
  ## Percent significant taxa
  percSig <- 100*(rowSums(pvals <= thresh)/ncol(pvals))
  graph_df <- data.frame("gamma" = gamma, "percSig" = percSig)
  
  if(blackWhite){
    p1 <- xyplot(percSig ~ gamma,
                 data = graph_df,
                 color = "black",
                 type = "l",
                 ylim = c(0,100),
                 xlab = list(label = "Gamma", cex = 1.5),
                 ylab = list(label = "Percent of Significant Entities", cex = cex),
                 lwd = 2,
                 par.settings = standard.theme(color = FALSE))
  } else{
    p1 <- xyplot(percSig ~ gamma,
                 data = graph_df,
                 color = "black",
                 type = "l",
                 ylim = c(0,100),
                 xlab = list(label = "Gamma", cex = 1.5),
                 ylab = list(label = "Percent of Significant Entities", cex = cex),
                 lwd = 2)
  }
  return(list(p1,p2))
}

gm <- function(x, na.rm = TRUE){
  exp(mean(log(x[x > 0]), na.rm=na.rm))
}

default.scale.model <- function(gamma, conds, p, mc.samples){
  ##adding a check to remove the intercept
  ##also identifying binary variables
  conds_used <- conds
  if(is.matrix(conds)){
    inds <- unname(which(colSums(conds == 1) == nrow(conds))) ##Find if any columns are the intercept
    if(length(inds) > 0){
      conds_used <- conds[,-inds]
    }
    binary <- unname(which(apply(conds_used,2, FUN = function(vec) length(unique(vec))) == 2))
  } else{
    ## testing if the one condition has more than two levels
    if(length(unique(conds_used)) == 2){
      ##Re-code to zero one
      binary <- TRUE
      conds_used <- as.numeric(as.factor(conds_used)) - 1
      conds_used <- ifelse(conds_used == 0, -1, conds_used)
    } else{
      binary <- FALSE
    }
  }
  
  ##centering and scaling the conditions
  conds_mat <- matrix(conds_used, nrow = length(p))
  conds_mat <- apply(conds_mat, 2, FUN = function(vec){
    if(length(unique(vec)) == 2){
      vec <- as.numeric(as.factor(vec)) - 1
      return(ifelse(vec == 0, -1, vec))
    } else{
      return(vec) 
    }
  })
  conds_mat <- apply(conds_mat, 2, FUN = function(vec){
    if(length(unique(vec)) > 2){
      return(scale(vec))
    } else{
      return(vec)
    }
  }) ##Centering
  
  ## scaling gamma to the scale of the conditions
  col_var <- rep(gamma^2, ncol(conds_mat))
  scale_samples <- matrix(NA, length(p), mc.samples) ## empty container

  for(i in 1:length(p)){
    geo_means <- log(apply(p[[i]],2,gm))
    noise <- sapply(col_var, FUN = function(sd){stats::rnorm(mc.samples, 0, sqrt(sd))})
    noise <- sweep(noise, 2, conds_mat[i,], "*")
    noise_mean <- rowSums(noise)
    
    scale_samples[i,] <- geo_means + noise_mean
  }
  scale_samples <- log2(exp(scale_samples))
  return(scale_samples)
}


#' Interpret the scale model implied by a certain level of gamma
#' 
#' @param clr A aldex.clr object
#' @return A table
#' 
#' @export
interpretGamma <- function(clr){
  ## Checking the scale samples were actual
  if(is.null(clr@scaleSamps)){
    stop("No scale samples passed in aldex.clr object. This function is for interpreting gamma if scale noise is added.")
  }
  
  ## checking conds to see if a t-test or glm type test
  if(is.vector(clr@conds)){
    scale.comb <- c(clr@scaleSamps)
    print("2.5th and 97.5th quantile of the scale samples...")
    return(data.frame("Lower_Bound" = quantile(scale.comb, c(0.025)), "Upper_Bound" = quantile(scale.comb, c(0.975))))
  } else{
    ##glm type
    ## creating a new clr object without scale
    p <- length(clr@dirichletData)
    
    ## creating a new object to calculate beta^parallel
    clr.par <- clr
    
    ## transforming the dirichlet samples
    ## adding the clr-transformed dirichlet data as the analysis data
    for(i in 1:p){
      clr.par@analysisData[[i]] <- log2(clr@dirichletData[[i]])
    }
    
    ##Fitting with scale
    conditions <- clr@conds
    k <- clr@mc.samples
    scale.r <- list()
    par.r <- list()
    perp.r <- list()
    
    mc.scale <- clr@analysisData
    for(i in 1:k){
      mci_lr <- t(sapply(mc.scale, function(x) x[, i]))
      scale.r[[i]] <- lr2glm(mci_lr, conditions)
    }
    
    ## 1. Subtract of Beta^\parallel
    ## 2. This gives Beta^\perp for each sample and variable
    ## 3. Summarize over entities
    
    ##Fitting without scale
    mc.par <- clr.par@analysisData
    for(i in 1:k){
      mci_lr <- t(sapply(mc.par, function(x) x[, i]))
      par.r[[i]] <- lr2glm(mci_lr, conditions)
      perp.r[[i]] <- scale.r[[i]] - par.r[[i]] ## For each mc sample, take the difference
    }
    
    ## for each mc sample, find the average beta^perp
    perp.names <- grepl("Estimate", names(perp.r[[1]]))
    perp.est <- matrix(nrow=k,ncol=sum(perp.names))
    for(i in 1:k){
      perp.est[i,] <- colMeans(perp.r[[i]][,perp.names])
    }
    
    lower.est <- apply(perp.est, 2, FUN = function(vec) quantile(vec, 0.025))
    upper.est <- apply(perp.est, 2, FUN = function(vec) quantile(vec, 0.975))
    
    return(data.frame("Variable" = names(perp.r[[1]][,perp.names]), "Lower_Bound" = lower.est, "Upper_Bound" = upper.est))
  }  
}
