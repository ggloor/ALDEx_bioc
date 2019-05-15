#' plotFeature.aldex
#'
#' Accepts the name of a feature and an aldex object to generate density plots
#' of within group dispersions, between group differences, and effect sizes.

aldex.plotFeature <- function(clrData, featureName, pooledOnly=FALSE,
                              densityOnly=FALSE) {
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
        plot(withinADensity, main="Within group clr values", xlab="clr Values", xlim=as.numeric(quantile(mixtureVector, probs=c(0.025,0.975))))
    } else {
        plot(withinBDensity, main="Within group clr values", xlab="clr Values", xlim=as.numeric(quantile(mixtureVector, probs=c(0.025,0.975))))
    }
    polygon(withinADensity, col=rgb(1,0,0,0.3), border="red")
    polygon(withinBDensity, col=rgb(0,0,1,0.3), border="cyan")
    if (!densityOnly) {
        par(new=T)
        boxplotAWidths = 0.25 * withinVectors[[1]]
        boxplotBWidths = 0.25 * withinVectors[[2]]
        widths <- c(maxDensityA / 6, maxDensityB / 6)
        boxplot(withinVectors[[1]], withinVectors[[2]], boxwex=widths, at=c(maxDensityA / 2, maxDensityB / 2), horizontal=TRUE, add=TRUE, col=c(rgb(1,0,0,0), rgb(0,0,1,0)))
    }

    if (!pooledOnly) {
        differenceVectorDensity <- density(differenceVector)
        plot(differenceVectorDensity, main="Diff. betw. groups", xlab="abs. difference", ylab="Density", xlim=as.numeric(quantile(differenceVector, probs=c(0.025,0.975))))
        polygon(differenceVectorDensity, col=rgb(1,0,1,0.3))
        if (!densityOnly) {
            boxplot(differenceVector, boxwex=max(differenceVectorDensity$y) / 6, at=max(differenceVectorDensity$y) / 2, horizontal=TRUE, add=TRUE, whisklty=0, staplelty=0, col=rgb(1,0,1,0), outline=FALSE)
        }

        dispersionVectorDensity <- density(dispersionVector)
        plot(dispersionVectorDensity, main="Absolute diff. betw. dispersion vectors", xlab="abs. difference", ylab="Density", xlim=as.numeric(quantile(dispersionVector, probs=c(0.025,0.975))))
        polygon(dispersionVectorDensity, col=rgb(1,1,0,0.3))
        if (!densityOnly) {
            boxplot(dispersionVector, boxwex=max(dispersionVectorDensity$y) / 6, at=max(dispersionVectorDensity$y) / 2, horizontal=TRUE, add=TRUE, whisklty=0, staplelty=0, col=rgb(1,1,0,0), outline=FALSE)
        }

        effectVectorDensity <- density(effectVector)
        plot(effectVectorDensity, main="Effect sizes", xlab="Effect size", ylab="Density", xlim=as.numeric(quantile(effectVector, probs=c(0.025,0.975))))
        polygon(effectVectorDensity, col=rgb(0,1,1,0.3))
        if (!densityOnly) {
            boxplot(effectVector, boxwex=max(effectVectorDensity$y) / 6, at=max(effectVectorDensity$y) / 2, horizontal=TRUE, add=TRUE, whisklty=0, staplelty=0, col=rgb(0,1,1,0), outline=FALSE)
        }
    }
}
