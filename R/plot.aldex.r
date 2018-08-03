aldex.plot <- function( x, ..., type=c("MW","MA"),
    xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,
    all.col=rgb(0,0,0,0.2), all.pch=19, all.cex=0.4,
    called.col="red", called.pch=20, called.cex=0.6,
    thres.line.col="darkgrey", thres.lwd=1.5,
    test="welch", cutoff=0.1, rare.col="black", rare=0,
    rare.pch=20, rare.cex=0.2
){
    type <- match.arg(type)

    if (length(x$effect) == 0) stop ("Please run aldex.effect before plotting")

    if (test == "welch"){
        if (length(x$we.eBH) == 0) stop ("Welch's t test results not in dataset")
        called <- x$we.eBH <= cutoff
    }else if (test == "wilcox"){
        if (length(x$wi.eBH) == 0) stop ("Wilcoxon test results not in dataset")
        called <- x$wi.eBH <= cutoff
    }
    if (test == "glm"){
        if (length(x$glm.eBH) == 0) stop ("glm test results not in dataset")
        called <- x$glm.eBH <= cutoff
    }
    if (test == "kruskal"){
        if (length(x$kw.eBH) == 0) stop ("Kruskall-Wallace test results not in dataset")
        called <- x$kw.eBH <= cutoff
    }
    if ( type == "MW" ) {
        if ( is.null(xlab) ) xlab <- expression( "Median" ~~ Log[2] ~~ "Dispersion" )
        if ( is.null(ylab) ) ylab <- expression( "Median" ~~ Log[2] ~~ "Difference" )

        plot(x$diff.win, x$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex)
        points(x$diff.win[x$rab.all < rare], x$diff.btw[x$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
        points(x$diff.win[called], x$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
        abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
        abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
        cols <- grep("rab.win", colnames(x))
        mtext(colnames(x)[cols[1]], 2, line=2, at=min(x$diff.btw), col="grey", cex=0.8)
        mtext(colnames(x)[cols[2]], 2, line=2, at=max(x$diff.btw), col="grey", cex=0.8)
    }
    if ( type == "MA" ) {

        if ( is.null(xlab) ) xlab <- expression( "Median" ~~ Log[2] ~~ "relative abundance" )
        if ( is.null(ylab) ) ylab <- expression( "Median" ~~ Log[2] ~~ "Difference" )

        plot(x$rab.all, x$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex)
        points(x$rab.all[x$rab.all < rare], x$diff.btw[x$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
        points(x$rab.all[called], x$diff.btw[called],  col=called.col, pch=called.pch, cex=called.cex)
        cols <- grep("rab.win", colnames(x))
        mtext(colnames(x)[cols[1]], 2, line=2, at=min(x$diff.btw), col="grey", cex=0.8)
        mtext(colnames(x)[cols[2]], 2, line=2, at=max(x$diff.btw), col="grey", cex=0.8)
     }
 }
