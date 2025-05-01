#' Plot an \code{aldex} Object
#'
#' @title Plot an \code{aldex} Object
#'
#' @description Create \code{MW}- or \code{MA}-type plots from the given \code{aldex} object.
#'
#' @param x an object produced by the \code{aldex.glm} function
#' @param eff an object produced by the \code{aldex.glm.effect} function
#' @param ... optional, unused arguments included for compatibility with the S3 method signature
#' @param type which type of plot is to be produced. \code{MA} is a 
#' Bland-Altman style plot; \code{MW} is an effect plot showing the
#' relationship of difference between and dispersion as described in:
#' http://dx.doi.org/10.1080/10618600.2015.1131161; \code{volcano} is a volcano plot
#' http://dx.doi.org/10.1186/gb-2003-4-4-210
#' @param contrast the column name of the model matrix contrast to plot
#' @param test the method of calculating significance, one of "pval" or "fdr" or "effect".
#' @param cutoff.pval the fdr cutoff, default 0.05
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
#' @param called.cex the character expansion of points with false discovery rate, q <= 0.05
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
#' @examples # See the examples for 'aldex.glm'
#' @export
aldex.glm.plot<-function (x, ..., eff = NULL, contrast=NULL, test = c('fdr', 'pval', 'effect'), 
	type = c("MW", "MA", "volcano"), xlab = NULL, ylab = NULL,
    xlim = NULL, ylim = NULL, all.col = rgb(0, 0, 0, 0.2), all.pch = 19,
    all.cex = 0.4, called.col = "red", called.pch = 20, called.cex = 0.6,
    thres.line.col = "darkgrey", thres.lwd = 1.5, cutoff.pval = 0.05, 
    cutoff.effect = 1, rare.col = "black", rare = 0, rare.pch = 20,rare.cex = 0.2)
{
    type <- match.arg(type)
    test <- match.arg(test)
    if (length(eff) == 0){
        stop("Please run aldex.glm.effect before plotting")
    }
    if (is.null(contrast)) {
        stop("Please enter a valid contrast name")
    }
    if (!contrast %in% names(eff)){
    	stop("Please enter a valid contrast name")
    }
    if(test == "fdr"){
    	column <- paste(contrast,":pval.padj", sep="")
    	called <- x[,column] <= cutoff.pval
    	all.p <- x[,column]
    } else if(test == 'pval'){
    	column <- paste(contrast,":pval", sep="")
    	called <- x[,column] <= cutoff.pval
    	all.p <- x[,column]
    }
    else if (test == "effect") {
        if (cutoff.effect <= 0.49)
            stop("Please set cutoff to at least 0.5")
        called <- abs(eff[[contrast]]$effect) >= cutoff.effect
    }
    if (type == "MW") {
        if (is.null(xlab))
            xlab <- expression("Median" ~ ~Log[2] ~ ~"Dispersion")
        if (is.null(ylab))
            ylab <- expression("Median" ~ ~Log[2] ~ ~"Difference")
        plot(eff[[contrast]]$diff.win, eff[[contrast]]$diff.btw, xlab = xlab, ylab = ylab,
            col = all.col, pch = all.pch, cex = all.cex, ylim=ylim, xlim=xlim)
        points(eff[[contrast]]$diff.win[eff[[contrast]]$rab.all < rare],
            eff[[contrast]]$diff.btw[eff[[contrast]]$rab.all < rare], 
            col = rare.col, pch = rare.pch, cex = rare.cex)
        points(eff[[contrast]]$diff.win[called], eff[[contrast]]$diff.btw[called], 
            col = called.col, pch = called.pch, cex = called.cex)
        abline(0, 1, col = thres.line.col, lty = 2, lwd = thres.lwd)
        abline(0, -1, col = thres.line.col, lty = 2, lwd = thres.lwd)
        cols <- grep("rab.win", colnames(eff$contrast))
        mtext(colnames(eff$contrast)[cols[1]], 2, line = 2, 
            at = min(eff[[contrast]]$diff.btw), col = "grey", cex = 0.8)
        mtext(colnames(eff$contrast)[cols[2]], 2, line = 2, 
            at = max(eff[[contrast]]$diff.btw), col = "grey", cex = 0.8)
    }
    if (type == "MA") {
        if (is.null(xlab))
            xlab <- expression("Median" ~ ~Log[2] ~ ~"relative abundance")
        if (is.null(ylab))
            ylab <- expression("Median" ~ ~Log[2] ~ ~"Difference")
        plot(eff[[contrast]]$rab.all, eff[[contrast]]$diff.btw, xlab = xlab, ylab = ylab,
            col = all.col, pch = all.pch, cex = all.cex, ylim=ylim, xlim=xlim)
        points(eff[[contrast]]$rab.all[eff[[contrast]]$rab.all < rare],
            eff[[contrast]]$diff.btw[eff[[contrast]]$rab.all < rare], 
            col = rare.col, pch = rare.pch, cex = rare.cex)
        points(eff[[contrast]]$rab.all[called], eff[[contrast]]$diff.btw[called], 
            col = called.col, pch = called.pch, cex = called.cex)
        cols <- grep("rab.win", colnames(eff$contrast))
        mtext(colnames(eff$contrast)[cols[1]], 2, line = 2, at = min(eff[[contrast]]$diff.btw),
            col = "grey", cex = 0.8)
        mtext(colnames(eff$contrast)[cols[2]], 2, line = 2, at = max(eff[[contrast]]$diff.btw),
            col = "grey", cex = 0.8)
    }
  if (type == "volcano") {
        if (is.null(ylab))
            ylab <- expression("-1 * Median Log"[10]~" q value")
        if (is.null(xlab))
            xlab <- expression("Median Log"[2]~" Difference")
        plot(eff[[contrast]]$diff.btw, -1*log10(all.p), xlab = xlab, ylab = ylab,
            col = all.col, pch = all.pch, cex = all.cex, ylim=ylim, xlim=xlim)
        points(eff[[contrast]]$diff.btw[called], -1*log10(all.p)[called], 
            col = called.col, pch = called.pch, cex = called.cex)
        cols <- grep("rab.win", colnames(eff$contrast))
        abline(v=1.5, col='grey', lty=2)
        abline(h=-1*log10(cutoff.pval), col='grey', lty=2)
        mtext(colnames(eff$contrast)[cols[1]], 1, line = 2, at = min(eff[[contrast]]$diff.btw),
            col = "grey", cex = 0.8)
        mtext(colnames(eff$contrast)[cols[2]], 1, line = 2, at = max(eff[[contrast]]$diff.btw),
            col = "grey", cex = 0.8)
    }
}
