\name{aldex.clr}
\alias{aldex.clr}
\alias{aldex.clr.function}
\alias{aldex.clr,data.frame-method}
\alias{aldex.clr,SummarizedExperiment-method}

\title{
    Compute an \code{aldex.clr} Object 
}
\description{
    Generate Monte Carlo samples of the Dirichlet distribution for each sample. Convert each instance using the centred log-ratio transform
    This is the input for all further analyses.
}
\usage{
aldex.clr(reads, mc.samples = 128, verbose = FALSE, useMC=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reads}{
        a non-negative, integer-only containing \code{data.frame} or \code{SummarizedExperiment} that has unique names for all rows and columns, where each row is a different gene and each column represents a sequencing read-count. Rows with 0 reads in each sample are deleted prior to analysis}
  \item{mc.samples}{
        the number of Monte Carlo samples to use to estimate the underlying distributions; since we are estimating central tendencies, 128 is usually sufficient
}
  \item{verbose}{
        Print diagnostic information while running. Useful only for debugging if fails on large datasets
}  
\item{useMC}{
    use multicore by default (FALSE). Multi core processing will be attempted with the BiocParallel package, then the parallel package. If neither are installed, serial processing will be used.
}
}
\details{
An explicit description of the input format for the \code{reads} object is shown under `Examples', below.}
\value{
    The object produced by the \code{clr} function contains the clr transformed values for each Monte-Carlo Dirichlet instance, which can be accessed through \code{getMonteCarloInstances(x)}, where \code{x} is the \code{clr} function output. Each list element is named by the sample ID. \code{getFeatures(x)} returns the features, \code{getSampleIDs(x)} returns sample IDs, and \code{getFeatureNames(x)} returns the feature names.
}
\references{
    Please use the citation given by \code{citation(package="ALDEx")}.}
\author{
    Greg Gloor, Ruth Grace Wond, Andrew Fernandes and Matt Links contributed to this code
}

\seealso{
    \code{\link{aldex.ttest}},
    \code{\link{aldex.glm}},
    \code{\link{aldex.effect}},
    \code{\link{selex}}
}
\examples{

    # The 'reads' data.frame or
    # SummarizedExperiment object should have
    # row and column names that are unique,
    # and looks like the following:
    #
    #              T1a T1b  T2  T3  N1  N2  Nx
    #   Gene_00001   0   0   2   0   0   1   0
    #   Gene_00002  20   8  12   5  19  26  14
    #   Gene_00003   3   0   2   0   0   0   1
    #   Gene_00004  75  84 241 149 271 257 188
    #   Gene_00005  10  16   4   0   4  10  10
    #   Gene_00006 129 126 451 223 243 149 209
    #       ... many more rows ...
    
    data(selex)
    x <- aldex.clr(selex, mc.samples = 2, verbose = FALSE)

}