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
    return(list(p = pt(t, df = df, lower.tail = FALSE), t= t))
    
  } else{
    t <- multtest::mt.teststat(data, as.numeric(grp1), test = "t", nonpara = "n")
    s1 <- apply(data[, grp1], 1, sd)
    s2 <- apply(data[, grp2], 1, sd)
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
    return(list(p = pt(t, df = df, lower.tail = FALSE), t= t))
  } 
}

# wilcox.fast function replaces wilcox.test
#  * runs much faster
#  * uses exact distribution for ties!
#    * this differs from ?wilcox.test
#  * optional paired test
#    * equivalent to wilcox.test(..., correct = FALSE)
#  * uses multtest
wilcox.fast <- function(data, group, paired){
  
  if(ncol(data) != length(group)) stop("Use rows for feature data.")
  grp1 <- group == unique(group)[1]
  grp2 <- group == unique(group)[2]
  n1 <- sum(grp1)
  n2 <- sum(grp2)
  
  # Check for ties in i-th Monte-Carlo instance
  data.t <- t(data)
  if(paired){
    anyTies <- any(apply(data.t[grp1, ] - data.t[grp2, ], 2,
                         function(x) length(unique(x))) != ncol(data) / 2)
  }else{
    anyTies <- any(apply(data.t, 2,
                         function(x) length(unique(x))) != ncol(data))
  }
  
  # Ties trigger slower, safer wilcox.test function
  if(anyTies){
    return(apply(data.t, 2, function(i){
      wilcox.test(x = i[grp1], y = i[grp2], paired = paired, correct = FALSE, alternative = "greater")$p.value}))
  }
  
  if(paired){
    
    if(n1 != n2) stop("Cannot pair uneven groups.")
    data.diff <- data.t[grp1, ] - data.t[grp2, ]
    V <- apply(data.diff, 2, function(x) sum(rank(abs(x))[x > 0]))
    topscore <- (n1 * (n1+1)) / 2
    V.lower <- ifelse(V > topscore / 2, topscore - V, V)
    if(sum(grp1) < 50){ # as per wilcox test, use exact -- ASSUMES NO TIES!!
      V.p <- psignrank(V.lower, n1) * 2
      return(ifelse(V.p > 1, 1, V.p)) # psignrank returns non-zero for W = mean
    }else{ # Use normal approximation
      V.std <- (topscore/2 - V.lower) / sqrt(n1*(n1 + 1) * (2*n1 + 1) / 24) # wilcox.test uses denom = 24
      return(pnorm(V.std, lower.tail = FALSE) * 2)
    }
    
    
  }else{
    
    W.std <- multtest::mt.teststat(data, as.numeric(grp1), test = "wilcoxon")
    if(sum(grp1) < 50 && sum(grp2) < 50){ # as per wilcox test, use exact -- ASSUMES NO TIES!!
      W.var <- sqrt((n1*n2) * (n1+n2+1) / 12)
      W <- abs(W.std) * W.var + (n1*n2) / 2
      W.p <- pwilcox(W - 1, n1, n2, lower.tail = FALSE) * 2
      return(ifelse(W.p > 1, 1, W.p)) # pwilcox returns non-zero for W = mean
    }else{ # Use normal approximation
      return(pnorm(abs(W.std), lower.tail = FALSE) * 2)
    }
  }
}
