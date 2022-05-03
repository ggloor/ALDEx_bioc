# returns the median clr abundances per sample, per condition
# returns the median differences in abundance between 2 conditions
# returns the median effect size and proportion of effect that overlaps 0
# data is returned in a data frame
# requires multicore
# this uses Rfast
aldex.effect <- function(clr, verbose=TRUE, include.sample.summary=FALSE, useMC=FALSE, CI=FALSE, glm.conds=NULL, paired.test=FALSE){

  # Use clr conditions slot instead of input
     if (is.vector(clr@conds)) {
       conditions <- clr@conds
     } else if (is.factor(clr@conds)) {
     	if (length(levels(clr@conds) == 2)) {
     	  conditions <- clr@conds
     	}
     } else if (is.matrix(clr@conds)){
       if(is.null(glm.conds)) stop("please provide a binary condition vector")
       conditions <- glm.conds
     } else {
       stop("please check that the conditions parameter for aldex.clr is correct.")
     }

    is.multicore = FALSE

    if ("BiocParallel" %in% rownames(installed.packages()) & useMC==TRUE){
        message("multicore environment is OK -- using the BiocParallel package")
        #require(BiocParallel)
        is.multicore = TRUE
    }
    else {
        if (verbose == TRUE) message("operating in serial mode")
    }

    nr <- numFeatures(clr) # number of features
    rn <- getFeatureNames(clr) # feature names
    # ---------------------------------------------------------------------

    # sanity check to ensure only two conditons passed to this function
    conditions <- as.factor( conditions )
    levels     <- levels( conditions )

    if ( length( conditions ) !=  numConditions(clr) ) stop("mismatch btw 'length(conditions)' and 'ncol(reads)'")

    if ( length( levels ) != 2 ) stop("only two condition levels are currently supported")

    levels <- vector( "list", length( levels ) )
    names( levels ) <- levels( conditions )

    for ( l in levels( conditions ) ) {
        levels[[l]] <- which( conditions == l )
        if ( length( levels[[l]] ) < 2 ) stop("condition level '",l,"' has less than two replicates")
    }

    # end sanity check
if (verbose == TRUE) message("sanity check complete")

    # Summarize the relative abundance (rab) win and all groups

    rab <- vector( "list", 3 )
    names(rab) <- c( "all", "win", "spl" )
    rab$win <- list()

    #this is the median value across all monte carlo replicates
    # for loops replaced with do.call
    cl2p <- NULL
    cl2p <- do.call(cbind, getMonteCarloInstances(clr))
    # this is a 2X speedup
    rab$all <- Rfast::rowMedians(cl2p)
    names(rab$all) <- rownames(cl2p)
    rm(cl2p)

 if (verbose == TRUE) message("rab.all  complete")

    for(level in levels(conditions)){
      cl2p <- NULL
      cl2p <- do.call(cbind, getMonteCarloInstances(clr)[levels[[level]]] )
      rab$win[[level]] <-  Rfast::rowMedians(cl2p)
      rm(cl2p)
 }

 if (verbose == TRUE) message("rab.win  complete")

    if (is.multicore == TRUE)  rab$spl <- bplapply( getMonteCarloInstances(clr), function(m) { t(apply( m, 1, median )) } )
    #RMV if (is.multicore == FALSE) rab$spl <- lapply( getMonteCarloInstances(clr), function(m) { t(apply( m, 1, median )) } )
    if (is.multicore == FALSE) rab$spl <- lapply( getMonteCarloInstances(clr), function(m) { Rfast::rowMedians(m) } )
if (verbose == TRUE) message("rab of samples complete")

    # ---------------------------------------------------------------------
    # Compute diffs btw and win groups

if (paired.test == FALSE ){
    l2d <- vector( "list", 2 )
    names( l2d ) <- c( "btw", "win" )
    l2d$win <- list()

    # abs( win-conditions diff ), btw smps
    #this generates a linear sample of the values rather than an exhaustive sample
    for ( level in levels(conditions) ) {
        concat <- NULL
        for ( l1 in sort( levels[[level]] ) ) {
            concat <- cbind(  getMonteCarloReplicate(clr,l1),concat )

        }

        #if the sample is huge, only sample 10000
        if ( ncol(concat) < 10000 ){
            sampl1 <- t(apply(concat, 1, function(x){sample(x, ncol(concat))}))
            sampl2 <- t(apply(concat, 1, function(x){sample(x, ncol(concat))}))
        } else {
            sampl1 <- t(apply(concat, 1, function(x){sample(x, 10000)}))
            sampl2 <- t(apply(concat, 1, function(x){sample(x, 10000)}))
        }
        l2d$win[[level]] <- cbind( l2d$win[[level]] , abs( sampl1 - sampl2 ) )
        rm(sampl1)
        rm(sampl2)
    }
if (verbose == TRUE) message("within sample difference calculated")

    # Handle the case when the groups have different spl sizes
    # get the minimum number of win spl comparisons
    ncol.wanted <- min( sapply( l2d$win, ncol ) )
    # apply multicore paradigm ML
    if (is.multicore == TRUE) l2d$win  <- bplapply( l2d$win, function(arg) { arg[,1:ncol.wanted] } )
    if (is.multicore == FALSE) l2d$win  <- lapply( l2d$win, function(arg) { arg[,1:ncol.wanted] } )

    # btw condition diff (signed)
    #get the btw condition as a random sample rather than exhaustive search
    concatl1 <- NULL
    concatl2 <- NULL
    for( l1 in levels[[1]] ) concatl1 <- cbind( getMonteCarloReplicate(clr,l1),concatl1 )
    for( l2 in levels[[2]] ) concatl2 <- cbind( getMonteCarloReplicate(clr,l2),concatl2 )

    sample.size <- min(ncol(concatl1), ncol(concatl2))

    if ( sample.size < 10000 ){
        smpl1 <- t(apply(concatl1, 1, function(x){sample(x, sample.size)}))
        smpl2 <- t(apply(concatl2, 1, function(x){sample(x, sample.size)}))
    } else {
        smpl1 <- t(apply(concatl1, 1, function(x){sample(x, 10000)}))
        smpl2 <- t(apply(concatl2, 1, function(x){sample(x, 10000)}))
    }
    l2d$btw <- smpl2 - smpl1

    rm(smpl1)
    rm(smpl2)

if (verbose == TRUE) message("between group difference calculated")

    win.max <- matrix( 0 , nrow=nr , ncol=ncol.wanted )
    l2d$effect <- matrix( 0 , nrow=nr , ncol=ncol(l2d$btw) )
    rownames(l2d$effect) <- rn

###the number of elements in l2d$btw and l2d$win may leave a remainder when
  #recycling these random vectors. Warnings are suppressed because this is not an issue
  #for this calculation. In fact, any attempt to get rid of this error would
  #decrease our power as one or both vectors would need to be truncated gg 20/06/2013

    options(warn=-1)

    for ( i in 1:nr ) {
        #mat <- rbind( l2d$win[[1]][i,] , l2d$win[[2]][i,] )
        #win.max[i,] <- apply( mat , 2 , max )
        # this is a 2x speedup
        win.max[i,] <- pmax( l2d$win[[1]][i,] , l2d$win[[2]][i,] )
        l2d$effect[i,] <- l2d$btw[i,] / win.max[i,]
        l2d$effect[i,][is.na(l2d$effect[i,])] <- 0
    }

    options(warn=0)

    rownames(win.max)   <- rn
    attr(l2d$win,"max") <- win.max
    rm(win.max)

    # ---------------------------------------------------------------------
    # Summarize diffs

    l2s <- vector( "list", 2 )
    names( l2s ) <- c( "btw", "win" )
    l2s$win <- list()

    #RMV l2s$btw <- t(apply( l2d$btw, 1, median ))
    l2s$btw <- Rfast::rowMedians(l2d$btw)
    #RMV l2s$win  <- t(apply( attr(l2d$win,"max"), 1, median ))
    l2s$win  <- Rfast::rowMedians( attr(l2d$win,"max"))
if (verbose == TRUE) message("group summaries calculated")

    if(CI == FALSE) {
    #  effect  <- t(apply( l2d$effect, 1, function(row){row[is.na(row)] <- 0 ; median(row) }))
      effect <- Rfast::rowMedians(l2d$effect)
    } else {
      effectlow <- t(apply( l2d$effect, 1, function(x) {x[is.na(x)] <- 0 ;
          quantile( x, probs=0.025, names=FALSE)} ))
      effecthigh <- t(apply( l2d$effect, 1, function(x) {x[is.na(x)] <- 0 ;
          quantile( x, probs=0.975, names=FALSE)} ))
      effect <- Rfast::rowMedians(l2d$effect)
    }
    overlap <- apply( l2d$effect, 1, function(row) { if(all(is.na(row))) warning("NAs in effect, ignore if using ALR");
                row[is.na(row)] <- 0 ;
                min( aitchison.mean( c( sum( row < 0 ) , sum( row > 0 ) ) + 0.5 ) ) } )
if (verbose == TRUE) message("unpaired effect size calculated")
} else if (paired.test == TRUE) {
  l2s <- vector( "list",2 )
  names( l2s ) <- c( "btw", "win" )
  
  
  sets <- names(levels)
  setA <- which(conditions == sets[1])
  setB <- which(conditions == sets[2])
  
  diff <- NULL
  for(i in 1:length(setA)){
    jnk1 <- getMonteCarloReplicate(clr,setA[i])
    jnk2 <- getMonteCarloReplicate(clr,setB[i])
    diff <- cbind(diff, jnk2-jnk1)
  }
  
  overlap <- apply( diff, 1, function(row) { if(all(is.na(row))) warning("NAs in effect, ignore if using ALR");
                row[is.na(row)] <- 0 ;
                min( aitchison.mean( c( sum( row < 0 ) , sum( row > 0 ) ) + 0.5 ) ) } )

  l2s$btw <- apply(diff, 1, mean)
  l2s$win <- apply(diff, 1, sd)
  
  effect <- l2s$btw/l2s$win
  if (verbose == TRUE) message("paired effect size calculated")

}

# make and fill in the data table
# i know this is inefficient, but it works and is not a bottleneck
   if(CI == FALSE | paired.test == TRUE) {
    rv <- list(
        rab = rab,
        diff = l2s,
        effect = effect,
        overlap = overlap
    )
    } else {
    rv <- list(
        rab = rab,
        diff = l2s,
        effect = effect,
        effectlow = effectlow,
        effecthigh = effecthigh,
        overlap = overlap
     )
    }

if (verbose == TRUE) message("summarizing output")

   y.rv <- data.frame(rv$rab$all)
   colnames(y.rv) <- c("rab.all")
   for(i in names(rv$rab$win)){
       nm <- paste("rab.win", i, sep=".")
       y.rv[,nm] <- data.frame(rv$rab$win[[i]])
   }
   if (include.sample.summary == TRUE){
    for(i in names(rv$rab$spl)){
       nm <- paste("rab.sample", i, sep=".")
       if (is.multicore == TRUE) y.rv[,nm] <- data.frame(t(rv$rab$spl[[i]]))
       if (is.multicore == FALSE) y.rv[,nm] <- data.frame(rv$rab$spl[[i]])
   }

   }
   for(i in names(rv$diff)){
       nm <- paste("diff", i, sep=".")
       y.rv[,nm] <- data.frame(rv$diff[[i]])
   }
   if(CI == FALSE | paired.test == TRUE) {
     y.rv[,"effect"] <- data.frame(rv$effect)
     y.rv[,"overlap"] <- data.frame(rv$overlap)
   } else {
     y.rv[,"effect"] <- data.frame(rv$effect)
     y.rv[,"effect.low"] <- data.frame(t(rv$effectlow))
     y.rv[,"effect.high"] <- data.frame(t(rv$effecthigh))
     y.rv[,"overlap"] <- data.frame(rv$overlap)
   }
    return(y.rv)

}
