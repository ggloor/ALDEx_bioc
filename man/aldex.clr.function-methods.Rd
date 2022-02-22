\name{aldex.clr}
\alias{aldex.clr}
\alias{aldex.clr.function}
\alias{aldex.clr,data.frame-method}
\alias{aldex.clr,matrix-method}
\alias{aldex.clr,RangedSummarizedExperiment-method}

\title{
    Compute an \code{aldex.clr} Object
}
\description{
    Generate Monte Carlo samples of the Dirichlet distribution for each sample.
    Convert each instance using a log-ratio transform.
    This is the input for all further analyses.
}
\usage{
    aldex.clr(reads, conds, mc.samples = 128, denom="all", 
          verbose=FALSE, useMC=FALSE, scale.samples = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{reads}{
        A \code{data.frame} or \code{RangedSummarizedExperiment} object containing
        non-negative integers only and with unique names for all rows and columns,
        where each row is a different gene and each column represents a sequencing
        read-count sample. Rows with 0 reads in each sample are deleted prior to
        analysis.
    }
    \item{conds}{
        A \code{vector} containing a descriptor for the samples, allowing them to
        be grouped and compared.
      }
      \item{mc.samples}{
        The number of Monte Carlo instances to use to estimate the underlying
        distributions; since we are estimating central tendencies, 128 is usually
        sufficient, but larger numbers may be .
    }
      \item{denom}{
          An \code{any} variable (all, iqlr, zero, lvha, median, user) indicating
          features to use as the denominator for the Geometric Mean calculation
          The default "all" uses the geometric mean abundance of all features.
          Using "median" returns the median abundance of all features.
          Using "iqlr" uses the features that are between the first and third
          quartile of the variance of the clr values across all samples.
          Using "zero" uses the non-zero features in each grop
          as the denominator. This approach is an extreme case where there are
          many nonzero features in one condition but many zeros in another. Using
          "lvha" uses features that have low variance (bottom quartile) and high
          relative abundance (top quartile in every sample). It is also
          possible to supply a vector of row indices to use as the denominator.
          Here, the experimentalist is determining a-priori which rows are thought
          to be invariant. In the case of RNA-seq, this could include ribosomal
          protein genes and and other house-keeping genes. This should be used
          with caution because the offsets may be different in the original data
          and in the data used by the function because features that are 0 in all
          samples are removed by \code{aldex.clr}.
      }
     \item{verbose}{
        Print diagnostic information while running. Useful only for debugging
        if fails on large datasets.
    }
    \item{useMC}{
        Use multicore by default (FALSE). Multi core processing will be attempted
        with the BiocParallel package. Serial processing will be used if this is
        not possible. In practice serial and multicore are nearly the same speed
        because of overhead in setting up the parallel processes.
    }
    \item{scale.samples}{
        Use scale simulation if not `NULL`. If a matrix is supplied, scale simulation
        will be used assuming that matrix denotes the scale samples. If a numeric is
        supplied, scale simulation will be applied by relaxing the geometric mean 
        assumption with the numeric representing the standard deviation of the
        scale distribution.
    }
}
\details{
    An explicit description of the input format for the \code{reads} object is
    shown under `Examples', below.}
\value{
    The object produced by the \code{clr} function contains the log-ratio transformed
    values for each Monte-Carlo Dirichlet instance, which can be accessed through
    \code{getMonteCarloInstances(x)}, where \code{x} is the \code{clr} function output.
    Each list element is named by the sample ID. \code{getFeatures(x)} returns the
    features, \code{getSampleIDs(x)} returns sample IDs, and \code{getFeatureNames(x)}
    returns the feature names.
}
\references{
    Please use the citation given by \code{citation(package="ALDEx")}.}
\author{
    Greg Gloor, Thom Quinn, Ruth Grace Wong, Andrew Fernandes, Matt Links and Jia Rong Wu contributed to this code.
}

\seealso{
    \code{\link{aldex.ttest}},
    \code{\link{aldex.glm}},
    \code{\link{aldex.effect}},
    \code{\link{selex}}
}
\examples{

    # The 'reads' data.frame or
    # RangedSummarizedExperiment object should
    # have row and column names that are unique,
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
    #subset for efficiency
    selex <- selex[1201:1600,]
    conds <- c(rep("NS", 7), rep("S", 7))
    x <- aldex.clr(selex, conds, mc.samples=2, denom="all", verbose=FALSE)

}
