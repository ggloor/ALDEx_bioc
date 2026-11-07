# t.fast function replaces t.test
#  * runs much faster
#  * assumes unequal variance
#  * optional paired test
#  * uses multtest
t.fast <- function(data, group, paired){
  
  grp1 <- group == unique(group)[1]
  grp2 <- group == unique(group)[2]
  n1 <- sum(grp1)
  n2 <- sum(grp2)
  
  if(paired){
    # Order pairs for the mt.teststat function
    if(n1 != n2) stop("Cannot pair uneven groups.")
    i.1 <- which(grp1)
    i.2 <- which(grp2)
    paired.order <- unlist(lapply(1:length(i.1), function(i) c(i.1[i], i.2[i])))
      
    t <- multtest::mt.teststat(data[, paired.order], as.numeric(grp1)[paired.order],
                                 test = "pairt", nonpara = "n")
    df <- length(i.1) - 1
    return(list(p = pt(t, df = df, lower.tail = FALSE),
                p.less = pt(t, df = df), t = t))
    
  } else{
    t <- multtest::mt.teststat(data, as.numeric(grp1), test = "t", nonpara = "n")
    s1 <- apply(data[, grp1, drop = FALSE], 1, sd)
    s2 <- apply(data[, grp2, drop = FALSE], 1, sd)
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    return(list(p = pt(t, df = df, lower.tail = FALSE),
                p.less = pt(t, df = df), t = t))
  } 
}

# wilcox.fast function replaces wilcox.test
#  * runs much faster
#  * exact tails include the observed statistic
#  * ties use wilcox.test with its tie-adjusted normal approximation
#  * optional paired test
#    * equivalent to wilcox.test(..., correct = FALSE)
#  * uses multtest
wilcox.fast <- function(data, group, paired, alternative = c("greater", "less")){
  alternative <- match.arg(alternative)
  upper <- alternative == "greater"
  
  if(ncol(data) != length(group)) stop("Use rows for feature data.")
  grp1 <- group == unique(group)[1]
  grp2 <- group == unique(group)[2]
  n1 <- sum(grp1)
  n2 <- sum(grp2)
  
  # Check for ties in i-th Monte-Carlo instance
  data.t <- t(data)
  if(paired){
    if(n1 != n2) stop("Cannot pair uneven groups.")
    differences <- data.t[grp1, , drop = FALSE] - data.t[grp2, , drop = FALSE]
    anyTies <- any(differences == 0) ||
      any(apply(abs(differences), 2, function(x) length(unique(x))) != n1)
  }else{
    anyTies <- any(apply(data.t, 2,
                         function(x) length(unique(x))) != ncol(data))
  }
  
  # Ties trigger slower, safer wilcox.test function
  if(anyTies){
    return(apply(data.t, 2, function(i){
      wilcox.test(x = i[grp1], y = i[grp2], paired = paired, correct = FALSE, alternative = alternative)$p.value}))
  }
  
  if(paired){
    
    if(n1 != n2) stop("Cannot pair uneven groups.")
    data.diff <- data.t[grp1, , drop = FALSE] - data.t[grp2, , drop = FALSE]
    V <- apply(data.diff, 2, function(x) sum(rank(abs(x))[x > 0]))
    topscore <- (n1 * (n1+1)) / 2
    #V.lower <- ifelse(V > topscore / 2, topscore - V, V)
    if(sum(grp1) < 50){ # as per wilcox test, use exact -- ASSUMES NO TIES!!
      V.p <- psignrank(if (upper) V-1 else V, n1, lower.tail = !upper)
      return(V.p) # psignrank returns non-zero for W = mean
    }else{ # Use normal approximation
      V.std <- (V - topscore/2) / sqrt(n1*(n1 + 1) * (2*n1 + 1) / 24) # wilcox.test uses denom = 24
      return(pnorm(V.std, lower.tail = !upper))
    }
    
    
  }else{
    
    W.std <- multtest::mt.teststat(data, as.numeric(grp1), test = "wilcoxon")
    if(sum(grp1) < 50 && sum(grp2) < 50){ # as per wilcox test, use exact -- ASSUMES NO TIES!!
      W.var <- sqrt((n1*n2) * (n1+n2+1) / 12)
      W <- W.std * W.var + (n1*n2) / 2
      W.p <- pwilcox(if (upper) W-1 else W, n1, n2, lower.tail = !upper)
      return(W.p) # pwilcox returns non-zero for W = mean
    }else{ # Use normal approximation
      return(pnorm(W.std, lower.tail = !upper))
    }
  }
}
