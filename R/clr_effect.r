# returns the median clr abundances per sample, per condition
# returns the median differences in abundance between 2 conditions
# returns the median effect size and proportion of effect that overlaps 0
# data is returned in a data frame
# requires multicore

aldex.effect <- function(clr, verbose=TRUE, include.sample.summary=FALSE, useMC=FALSE, CI=FALSE){

  # Use clr conditions slot instead of input
    conditions <- clr@conds

    is.multicore = FALSE

    if ("BiocParallel" %in% rownames(installed.packages()) & useMC){
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


    # end sanity check
if (verbose == TRUE) message("sanity check complete")

    # Summarize the relative abundance (rab) win and all groups

    rab <- vector( "list", 3 )
    names(rab) <- c( "all", "win", "spl" )
    rab$win <- list()

    #this is the median value across all monte carlo replicates
    cl2p <- NULL
    for ( m in getMonteCarloInstances(clr) ) cl2p <- cbind( cl2p, m )
    rab$all <- t(apply( cl2p, 1, median ))
    rm(cl2p)
    gc()
 if (verbose == TRUE) message("rab.all  complete")

    #this is the median value across all monte carlo replicates per level
    for ( level in levels(conditions) ) {
        cl2p <- NULL
        for ( i in levels[[level]] ) cl2p <- cbind( cl2p, getMonteCarloSample(clr,i) )
        rab$win[[level]] <- t(apply( cl2p, 1, median ))
        rm(cl2p)
        gc()
    }
 if (verbose == TRUE) message("rab.win  complete")

    if (is.multicore == TRUE)  rab$spl <- bplapply( getMonteCarloInstances(clr), function(m) { t(apply( m, 1, median )) } )
    if (is.multicore == FALSE) rab$spl <- lapply( getMonteCarloInstances(clr), function(m) { t(apply( m, 1, median )) } )

if (verbose == TRUE) message("rab of samples complete")

    # ---------------------------------------------------------------------
    # Compute diffs btw and win groups

    l2d <- vector( "list", 2 )
    names( l2d ) <- c( "btw", "win" )
    l2d$win <- list()

    # abs( win-conditions diff ), btw smps
#this generates a linear sample of the values rather than an exhaustive sample
    for ( level in levels(conditions) ) {
        concat <- NULL
        for ( l1 in sort( levels[[level]] ) ) {
            concat <- cbind(  getMonteCarloSample(clr,l1),concat )

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
        gc()
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
    for( l1 in levels[[1]] ) concatl1 <- cbind( getMonteCarloSample(clr,l1),concatl1 )
    for( l2 in levels[[2]] ) concatl2 <- cbind( getMonteCarloSample(clr,l2),concatl2 )

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
    gc()
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
        win.max[i,] <- apply( ( rbind( l2d$win[[1]][i,] , l2d$win[[2]][i,] ) ) , 2 , max )
        l2d$effect[i,] <- l2d$btw[i,] / win.max[i,]
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

    l2s$btw <- t(apply( l2d$btw, 1, median ))
    l2s$win  <- t(apply( attr(l2d$win,"max"), 1, median ))
if (verbose == TRUE) message("group summaries calculated")

    if(CI == FALSE) {
      effect  <- t(apply( l2d$effect, 1, median ))
    } else {
      effectlow <- t(apply( l2d$effect, 1, function(x) quantile(x, probs=0.025, names=FALSE) ))
      effecthigh <- t(apply( l2d$effect, 1, function(x) quantile(x, probs=0.975, names=FALSE) ))
      effect  <- t(apply( l2d$effect, 1, median ))
    }
    overlap <- apply( l2d$effect, 1, function(row) { min( aitchison.mean( c( sum( row < 0 ) , sum( row > 0 ) ) + 0.5 ) ) } )
if (verbose == TRUE) message("effect size calculated")

# make and fill in the data table
# i know this is inefficient, but it works and is not a bottleneck
   if(CI == FALSE) {
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

   y.rv <- data.frame(t(rv$rab$all))
   colnames(y.rv) <- c("rab.all")
   for(i in names(rv$rab$win)){
       nm <- paste("rab.win", i, sep=".")
       y.rv[,nm] <- data.frame(t(rv$rab$win[[i]]))
   }
   if (include.sample.summary == TRUE){
    for(i in names(rv$rab$spl)){
       nm <- paste("rab.sample", i, sep=".")
       y.rv[,nm] <- data.frame(t(rv$rab$spl[[i]]))
   }

   }
   for(i in names(rv$diff)){
       nm <- paste("diff", i, sep=".")
       y.rv[,nm] <- data.frame(t(rv$diff[[i]]))
   }
   if(CI == FALSE) {
     y.rv[,"effect"] <- data.frame(t(rv$effect))
     y.rv[,"overlap"] <- data.frame(rv$overlap)
   } else {
     y.rv[,"effect"] <- data.frame(t(rv$effect))
     y.rv[,"effect.low"] <- data.frame(t(rv$effectlow))
     y.rv[,"effect.high"] <- data.frame(t(rv$effecthigh))
     y.rv[,"overlap"] <- data.frame(rv$overlap)
   }
    return(y.rv)

}
