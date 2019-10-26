aldex.expectedDistance <- function(clrData) {
    if (class(clrData) != "aldex.clr" && class(clrData) != "list") {
        stop("Please supply a valid aldex.clr object.")
    }

    d.mc <- getMonteCarloInstances(clrData)
    distances <- array(0, dim=c(length(getSampleIDs(clrData)), length(getSampleIDs(clrData)), clrData@mc.samples))

    message("computing distances for each instance...")
    # Transpose call makes each row a sample and each col a feature for dist()
    for (i in 1:clrData@mc.samples) {
        distances[,,i] <- distances[,,i] + as.matrix(dist(t(sapply(d.mc, function(x) {x[,i]}))))
    }

    # Apply median() across all distance-value instances for each sample pair
    message("computing median distance across instances...")
    expectedDist <- as.matrix(apply(distances, c(1, 2), median))

    rownames(expectedDist) <- getSampleIDs(clrData)
    colnames(expectedDist) <- getSampleIDs(clrData)

    return(as.dist(expectedDist))
}
