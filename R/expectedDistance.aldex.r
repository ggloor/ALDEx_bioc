aldex.expectedDistance <- function(d.clr) {
    if (class(d.clr) != "aldex.clr" && class(d.clr) != "list") {
        stop("Please supply a valid aldex.clr object.")
    }

    d.mc <- getMonteCarloInstances(d.clr)
    distances <- array(0, dim=c(length(getSampleIDs(d.clr)), length(getSampleIDs(d.clr)), d.clr@mc.samples))

    message("computing distances for each instance...")
    # Transpose call makes each row a sample and each col a feature for dist()
    for (i in 1:d.clr@mc.samples) {
        distances[,,i] <- distances[,,i] + as.matrix(dist(t(sapply(d.mc, function(x) {x[,i]}))))
    }

    # Apply median() across all distance-value instances for each sample pair
    message("computing median distance across instances...")
    expectedDist <- as.matrix(apply(distances, c(1, 2), median))

    rownames(expectedDist) <- getSampleIDs(d.clr)
    colnames(expectedDist) <- getSampleIDs(d.clr)
    
    return(as.dist(expectedDist))
}
