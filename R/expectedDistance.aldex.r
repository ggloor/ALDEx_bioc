#' Calculate the expected values of distances between samples, given an \code{aldex} Object
#'
#' Calculates the expected value of distances between samples, given an \code{aldex} Object, using the median value of distances derived from n Monte-Carlo replicates.
#'
#' @param clrData an object of class \code{aldex} produced by the \code{aldex} function
#'
#' @return Returns a \code{dist} Object.
#'
#' @references Please use the citation given by \code{citation(package="ALDEx")}.
#'
#' @seealso \code{\link{aldex}}, \code{\link{aldex.clr}}, \code{dist}
#'
#' @examples
#' data(selex)
#' #subset for efficiency
#' selex <- selex[1201:1600,]
#' conds <- c(rep("NS", 7), rep("S", 7))
#' x <- aldex.clr(selex, conds, mc.samples = 128, denom = "all", verbose = FALSE)
#' x.dist <- aldex.expectedDistance(x)
#' plot(hclust(x.dist))
#'
#' @export
aldex.expectedDistance <- function(clrData) {
    if (class(clrData)[1] != "aldex.clr") {
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
