\name{aldex.expectedDistance}
\alias{aldex.expectedDistance}
\title{
    Calculate the expected values of distances between samples, given an \code{aldex} Object
}
\description{
    Calculates the expected value of distances between samples, given an \code{aldex} Object, using the median value of distances derived from n Monte-Carlo replicates.
}
\usage{
    \method{aldex}{expectedDistance}(clrData)
}
\arguments{
    \item{clrData}{
        an object of class \code{aldex} produced by the \code{aldex} function
    }
}
\details{
    Generates a distance matrix for each Monte-Carlo instance in an \code{aldex} Object. Calculates the median distance value across all instances.
}
\value{
    Returns a \code{dist} Object.
}
\references{
    Please use the citation given by \code{citation(package="ALDEx")}.
}
\seealso{
    \code{\link{aldex}},
    \code{\link{aldex.clr}}
    \code{dist}
}
\examples{
    data(selex)
        #subset for efficiency
        selex <- selex[1201:1600,]
    conds <- c(rep("NS", 7), rep("S", 7))
    x <- aldex.clr(selex, conds, mc.samples = 128, denom = "all", verbose = FALSE)
    x.dist <- aldex.expectedDistance(x)
}
