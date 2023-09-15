#' Show dispersion of the expected values returned by \code{aldex.effect}
#'
#' \code{aldex.plotFeature} generates density plots showing the dispersion
#'  of the expected values given in the output from \code{aldex.effect}.
#'  The expected values are shown in the plots. This is a diagnostic
#'  visualization to help determine if the expected values are trustworthy
#'
#' @param clrData the output object from \code{aldex.clr}
#' @param featureName the name of the feature from the input data
#' @param pooledOnly show only the pooled plots, default FALSE, shows all plots
#' @param densityOnly show only the density plots, default FALSE includes expected values
#' @param reset.par reset the plotting parameter to par(c(1,1)), default FALSE
#'
#' @author Brandon Lieng, Greg Gloor
#'
#' @seealso
#'  \code{\link{aldex.clr}},
#'  \code{\link{aldex.effect}},
#'  \code{\link{selex}}
#'
#' @references Please use the citation given by
#'  \code{citation(package="ALDEx2")}.
#'
#' @examples
#' data(selex)
#' #subset for efficiency
#' selex <- selex[1201:1600,]
#' conds <- c(rep("NS", 7), rep("S", 7))
#' x <- aldex.clr(selex, conds, mc.samples=4, denom="all")
#' aldex.plotFeature(x, "S:D:A:D")

aldex.plotFeature <- function(clrData, featureName, pooledOnly=FALSE,
                              densityOnly=FALSE, reset.par=FALSE) {
    mcInstances <- getMonteCarloInstances(clrData)
    mcInstancesByGroup <- split(mcInstances, factor(clrData@conds))

    # Generating within group pooled vectors
    withinVectors <- list()
    firstGroup <- TRUE
    for (group in mcInstancesByGroup) {
        featureClrVals <- matrix(ncol=clrData@mc.samples)
        for (sample in group) {
            featureClrVals <- rbind(featureClrVals, subset(sample, rownames(sample) %in% featureName))
        }
        if (nrow(featureClrVals) == 1) {
            stop("Feature not found in at least one group")
        }
        # Discard first row of NAs used to initialize matrix
        featureClrVals <- as.vector(featureClrVals[2:nrow(featureClrVals),])
        if (firstGroup) {
            withinVectors[[names(mcInstancesByGroup[1])]] <- featureClrVals
            firstGroup <- !firstGroup
        } else {
            withinVectors[[names(mcInstancesByGroup[2])]] <- featureClrVals
        }
    }

    # Generating difference, dispersion, and effect vectors
    differenceVector <- withinVectors[[1]] - withinVectors[[2]]
    mixtureVector <- c(withinVectors[[1]], withinVectors[[2]])
    if (!pooledOnly) {
        dispersionVector <- pmax(
                            abs(withinVectors[[1]] - sample(withinVectors[[1]])),
                            abs(withinVectors[[2]] - sample(withinVectors[[2]]))
                            )
        effectVector <- differenceVector / dispersionVector
        par(mfcol=c(2,2))
    }

    # Plotting...
    withinADensity <- density(withinVectors[[1]])
    withinBDensity <- density(withinVectors[[2]])
    maxDensityA <- max(withinADensity$y)
    maxDensityB <- max(withinBDensity$y)
    # Block here ensures axis is scaled to the largest y-value in the pool
    if (maxDensityA > maxDensityB) {
        plot(withinADensity, main="Group Distribution", xlab="clr Values", xlim=as.numeric(quantile(mixtureVector, probs=c(0.025,0.975))))
### Add in when bored
#    if(hist==T){
#      hist(withinVectors[[1]], freq=F, breaks=19, 
#         border=NULL, add=T, col=rgb(1,0,1,0.2))
#      hist(withinVectors[[2]], freq=F, breaks=19,
#         border=NULL, add=T, col=rgb(0,0,1,0.2))
#    }

    } else {
        plot(withinBDensity, main="Group Distribution", xlab="clr Values", xlim=as.numeric(quantile(mixtureVector, probs=c(0.025,0.975))))
    }
    polygon(withinADensity, col=rgb(1,0,0,0.3), border="red")
    polygon(withinBDensity, col=rgb(0,0,1,0.3), border="cyan")
    if (!densityOnly) {
        par(new=T)
        boxplotAWidths = 0.25 * withinVectors[[1]]
        boxplotBWidths = 0.25 * withinVectors[[2]]
        widths <- c(maxDensityA / 6, maxDensityB / 6)
        boxplot(withinVectors[[1]], withinVectors[[2]], boxwex=widths, at=c(maxDensityA / 2, maxDensityB / 2), horizontal=TRUE, add=TRUE, whisklty=0, staplelty=0, col=c(rgb(1,0,0,0), rgb(0,0,1,0)))
    }

    if (!pooledOnly) {
        differenceVectorDensity <- density(differenceVector)
        plot(differenceVectorDensity, main="Btw Grp Difference", xlab="Difference", ylab="Density", xlim=as.numeric(quantile(differenceVector, probs=c(0.025,0.975))))
        polygon(differenceVectorDensity, col=rgb(1,0,1,0.3))
        if (!densityOnly) {
            boxplot(differenceVector, boxwex=max(differenceVectorDensity$y) / 6, at=max(differenceVectorDensity$y) / 2, horizontal=TRUE, add=TRUE, whisklty=0, staplelty=0, col=rgb(1,0,1,0), outline=FALSE)
        }

        dispersionVectorDensity <- density(dispersionVector)
        plot(dispersionVectorDensity, main="Win Grp Dispersion", xlab="abs. Dispersion", ylab="Density", xlim=as.numeric(quantile(dispersionVector, probs=c(0.025,0.975))))
        polygon(dispersionVectorDensity, col=rgb(1,1,0,0.3))
        if (!densityOnly) {
            boxplot(dispersionVector, boxwex=max(dispersionVectorDensity$y) / 6, at=max(dispersionVectorDensity$y) / 2, horizontal=TRUE, add=TRUE, whisklty=0, staplelty=0, col=rgb(1,1,0,0), outline=FALSE)
        }

        effectVectorDensity <- density(effectVector)
        plot(effectVectorDensity, main="Effect size", xlab="Effect size", ylab="Density", xlim=as.numeric(quantile(effectVector, probs=c(0.025,0.975))))
        polygon(effectVectorDensity, col=rgb(0,1,1,0.3))
        if (!densityOnly) {
            boxplot(effectVector, boxwex=max(effectVectorDensity$y) / 6, at=max(effectVectorDensity$y) / 2, horizontal=TRUE, add=TRUE, whisklty=0, staplelty=0, col=rgb(0,1,1,0), outline=FALSE)
        }
    }
    if(reset.par == T) {par(mfrow=c(1,1))}
}
