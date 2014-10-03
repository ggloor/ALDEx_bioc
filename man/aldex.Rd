\name{aldex}
\alias{aldex}
\title{
    Compute an \code{aldex} Object 
}
\description{
    Generate Monte Carlo samples of the Dirichlet distribution for each sample. Convert each instance using the centred log-ratio transform
    Return two sample test values (Welch's t, Wilcoxon) or multi-sample test values (glm or Kruskal Wallace). Returns effect size values by defaulte.
}
\usage{
aldex(reads, conditions, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{reads}{
        a non-negative, integer-only containing \code{data.frame} that has unique names for all rows and columns, where each row is a different gene and each column represents a sequencing read-count. Rows with 0 
        reads in each sample are deleted prior to analysis}
  \item{mc.samples}{
        the number of Monte Carlo samples to use to estimate the underlying distributions; since we are estimating central tendencies, 128 is usually sufficient
}
  \item{conditions}{
    a description of the data structure to be used for testing
}
  \item{test}{
    which tests to perform: t = Welch's t and Wilcoxon, glm = Kruskal Wallace and glm
}
  \item{effect}{
    calculate abundances and effect sizes
}
    \item{include.sample.summary}{
        include median clr values for each sample, defaults to FALSE
    }

  \item{verbose}{
        Print diagnostic information while running. Useful only for debugging if fails on large datasets
}
}
\details{
An explicit description of the input format for the \code{reads} object is shown under `Examples', below.
This is not intended to be the generic function. The system is intended to be used for demonstration or instructional purposes.
}

\value{returns a number of values that depends on the set of options. See the
    return values of \code{aldex.ttest}, \code{aldex.glm}, and \code{aldex.effect} for explanations and example}

\references{
    Please use the citation given by \code{citation(package="ALDEx")}.}
\author{
    Greg Gloor, Andrew Fernandes and Matt Links contributed to this code
}

\seealso{
    \code{\link{aldex.ttest}},
    \code{\link{aldex.glm}},
    \code{\link{aldex.effect}},
    \code{\link{aldex.corr}},
    \code{\link{selex}}
}
\examples{

    # The 'reads' data.frame should have row
    # and column names that are unique, and
    # looks like the following:
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
    conds <- c(rep("N", 7), rep("S",7))
    x <- aldex(selex, conds, mc.samples = 2, test="t", effect=FALSE, verbose = FALSE)

}