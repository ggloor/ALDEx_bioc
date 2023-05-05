#' Plot an \code{aldex} Object
#'
#' @title Plot an \code{aldex} Object
#'
#' @description Create \code{MW}- or \code{MA}-type plots from the given \code{aldex} object.
#'
#' @param x an object of class \code{aldex} produced by the \code{aldex} function
#' @param ... optional, unused arguments included for compatibility with the S3 method signature
#' @param type which type of plot is to be produced. \code{MA} is a Bland-Altman style plot; \code{MW} is a
#' difference between to a variance within plot as described in:
#' http://dx.doi.org/10.1080/10618600.2015.1131161; \code{volcano} is a volcano plot
#' of either the difference or variance type: http://dx.doi.org/10.1186/gb-2003-4-4-210
#' @param test the method of calculating significance, one of:
#' \code{ppp} = posterior predictive p value;
#' \code{welch} = welch's t test;
#' \code{wilcox} = wilcox rank test;
#' \code{glm} = glm;
#' \code{kruskal} = Kruskal-Wallace test;
#' \code{effect} = effect size
#' @param cutoff.pval the Benjamini-Hochberg fdr cutoff, default 0.1
#' @param cutoff.effect the effect size cutoff for plotting, default 1
#' @param xlab the x-label for the plot, as per the parent \code{plot} function
#' @param ylab the y-label for the plot, as per the parent \code{plot} function
#' @param xlim the x-limits for the plot, as per the parent \code{plot} function
#' @param ylim the y-limits for the plot, as per the parent \code{plot} function
#' @param all.col the default colour of the plotted points
#' @param all.pch the default plotting symbol
#' @param all.cex the default symbol size
#' @param called.col the colour of points with false discovery rate, q <= 0.1
#' @param called.pch the symbol of points with false discovery rate, q <= 0.1
#' @param called.cex the character expansion of points with false discovery rate, q <= 0.1
#' @param thres.line.col the colour of the threshold line where within and between group variation is equivalent
#' @param thres.lwd the width of the threshold line where within and between group variation is equivalent
#' @param rare relative abundance cutoff for rare features, default 0 or the mean abundance
#' @param rare.col color for rare features, default black
#' @param rare.pch the default symbol of rare features
#' @param rare.cex the default symbol size of rare points
#'
#' @details This particular specialization of the \code{plot} function is relatively simple and provided for convenience.
#' For more advanced control of the plot is is best to use the values returned by \code{summary(x)}.
#'
#' @return None.
#'
#' @references Please use the citation given by \code{citation(package="ALDEx")}.
#'
#' @seealso \code{\link{aldex}}, \code{\link{aldex.effect}}, \code{\link{aldex.ttest}}, \code{\link{aldex.glm}}
#'
#' @examples # See the examples for 'aldex'
#' @export
aldex.plot<-function (x, ..., type = c("MW", "MA", "volcano", "volcano.var"), xlab = NULL, ylab = NULL,
    xlim = NULL, ylim = NULL, all.col = rgb(0, 0, 0, 0.2), all.pch = 19,
    all.cex = 0.4, called.col = "red", called.pch = 20, called.cex = 0.6,
    thres.line.col = "darkgrey", thres.lwd = 1.5, test = "welch",
    cutoff.pval = 0.1, cutoff.effect = 1, rare.col = "black", rare = 0, rare.pch = 20,
    rare.cex = 0.2)
{
    type <- match.arg(type)
    if (length(x$effect) == 0)
        stop("Please run aldex.effect before plotting")
    if (test == "welch" | test == 'ppp') {
        if (length(x$we.eBH) == 0 & length(x$p.val) == 0)
            stop("t test results not in dataset")
        if ( length(x$we.eBH) > 0 ){ 
        	warning('using we.eBH') 
        	called <- x$we.eBH <= cutoff.pval
       		all.p <- x$we.eBH
        } else if ( length(x$p.val) > 0 ){ 
            warning("using ppp value")
        	if(cutoff.pval > 1e-2) {
        	  warning("posterior p value cutoff changed to 1e-2")
        	  # note the called group is set here
        	  # before the 0 values are modified below
        	  called <- x$p.val <= 1e-2
        	  } else { called <- x$p.val <= cutoff.pval
        	  }
        	all.p <- x$p.val
        	# modifying 0 values for the ppp
        	adj.val <- min(all.p[all.p>0])/2
  			if(adj.val > 1e-3){adj.val = 1e-3}
  			all.p[all.p == 0] <- adj.val
        }
    }
    else if (test == "wilcox") {
        if (length(x$wi.eBH) == 0)
            stop("Wilcoxon test results not in dataset")
        called <- x$wi.eBH <= cutoff.pval
        all.p <- x$wi.eBH
    }
    else if (test == "glm") {
        if (length(x$glm.eBH) == 0)
            stop("glm test results not in dataset")
        called <- x$glm.eBH <= cutoff.pval
        all.p <- x$glm.eBH
    }
    else if (test == "kruskal") {
        if (length(x$kw.eBH) == 0)
            stop("Kruskall-Wallace test results not in dataset")
        called <- x$kw.eBH <= cutoff.pval
        all.p <- x$kw.eBH
    }
    else if (test == "effect") {
        if (cutoff.effect <= 0.49)
            stop("Please set cutoff to at least 0.5")
        called <- abs(x$effect) >= cutoff.effect
    }
    else if (test == "both") {
        if (cutoff.effect <= 0.49)
            stop("Please set cutoff to at least 0.5")
        called <- abs(x$effect) >= cutoff.effect & x$we.eBH <= cutoff.pval
    }
    if (type == "MW") {
        if (is.null(xlab))
            xlab <- expression("Median" ~ ~Log[2] ~ ~"Dispersion")
        if (is.null(ylab))
            ylab <- expression("Median" ~ ~Log[2] ~ ~"Difference")
        plot(x$diff.win, x$diff.btw, xlab = xlab, ylab = ylab,
            col = all.col, pch = all.pch, cex = all.cex, ylim=ylim, xlim=xlim)
        points(x$diff.win[x$rab.all < rare], x$diff.btw[x$rab.all <
            rare], col = rare.col, pch = rare.pch, cex = rare.cex)
        points(x$diff.win[called], x$diff.btw[called], col = called.col,
            pch = called.pch, cex = called.cex)
        abline(0, 1, col = thres.line.col, lty = 2, lwd = thres.lwd)
        abline(0, -1, col = thres.line.col, lty = 2, lwd = thres.lwd)
        cols <- grep("rab.win", colnames(x))
        mtext(colnames(x)[cols[1]], 2, line = 2, at = min(x$diff.btw),
            col = "grey", cex = 0.8)
        mtext(colnames(x)[cols[2]], 2, line = 2, at = max(x$diff.btw),
            col = "grey", cex = 0.8)
    }
    if (type == "MA") {
        if (is.null(xlab))
            xlab <- expression("Median" ~ ~Log[2] ~ ~"relative abundance")
        if (is.null(ylab))
            ylab <- expression("Median" ~ ~Log[2] ~ ~"Difference")
        plot(x$rab.all, x$diff.btw, xlab = xlab, ylab = ylab,
            col = all.col, pch = all.pch, cex = all.cex, ylim=ylim, xlim=xlim)
        points(x$rab.all[x$rab.all < rare], x$diff.btw[x$rab.all <
            rare], col = rare.col, pch = rare.pch, cex = rare.cex)
        points(x$rab.all[called], x$diff.btw[called], col = called.col,
            pch = called.pch, cex = called.cex)
        cols <- grep("rab.win", colnames(x))
        mtext(colnames(x)[cols[1]], 2, line = 2, at = min(x$diff.btw),
            col = "grey", cex = 0.8)
        mtext(colnames(x)[cols[2]], 2, line = 2, at = max(x$diff.btw),
            col = "grey", cex = 0.8)
    }
  if (type == "volcano") {
        if (is.null(ylab))
            ylab <- expression("-1 * Median Log"[10]~" q value")
        if (is.null(xlab))
            xlab <- expression("Median Log"[2]~" Difference")
        plot(x$diff.btw, -1*log10(all.p), xlab = xlab, ylab = ylab,
            col = all.col, pch = all.pch, cex = all.cex, ylim=ylim, xlim=xlim)
        points(x$diff.btw[called], -1*log10(all.p)[called], col = called.col,
            pch = called.pch, cex = called.cex)
        cols <- grep("rab.win", colnames(x))
        mtext(colnames(x)[cols[1]], 1, line = 2, at = min(x$diff.btw),
            col = "grey", cex = 0.8)
        mtext(colnames(x)[cols[2]], 1, line = 2, at = max(x$diff.btw),
            col = "grey", cex = 0.8)

    }
  if (type == "volcano.var") {
       if (is.null(ylab))
            ylab <- expression("-1 * Median Log"[10]~" q value")
        if (is.null(xlab))
            xlab <- expression("Median Log"[2]~" Dispersion")
        plot(x$diff.win, -1*log10(all.p), xlab = xlab, ylab = ylab,
            col = all.col, pch = all.pch, cex = all.cex, ylim=ylim, xlim=xlim)
        points(x$diff.win[called], -1*log10(all.p)[called], col = called.col,
            pch = called.pch, cex = called.cex)
        cols <- grep("rab.win", colnames(x))

    }
}
