\name{aldex.set.mode}
\alias{aldex.set.mode}
\title{
   identify set of denominator features for log-ratio calculation
}
\description{
  calculate the features that are to be used as the denominator for
  the Geometric Mean calculation in clr_function.R}
\usage{
  aldex.set.mode(reads, conds, denom="all",verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{reads}{
    A data frame containing the samples and features per sample.
}
    \item{conds}{
    A vector describing which samples belong to what condition.
}
    \item{denom}{
    Character argument specifying which indicies to return.
    'all' returns all features in both conditons.
    'zero' returns the nonzero count features per condition.
    'iqlr' returns the features whose variance falls within the
      inter-quantile range of the CLR-transformed data.
    In cases of malformed or null queries, input defaults to 'all'.
    Additionally, the input can be a numeric vector, which contains
    a set of row indicies to center the data against. Only for advanced
    users who can pre-determine the invariant set of features within
    their data. Check that the same number of features are in the
    input and output datasets.
}
    \item{verbose}{
    Logical argument controlling messages, default=FALSE 
}
}
\details{
    An explicit example for two conditions is shown in the `Examples' below.
}
\value{
    Outputs a vector containing indicies per condition, or a single
    vector in some cases.
}
\references{
    Please use the citation given by \code{citation(package="ALDEx")}.
}
\author{
    Jia Rong Wu
}

\seealso{
    \code{\link{aldex.clr}},
    \code{\link{aldex.ttest}},
    \code{\link{aldex.effect}},
    \code{\link{selex}}
}
\examples{
    # x is the output of the \code{x <- clr(data, mc.samples)} function
    # conditions is a description of the data
    # for the selex dataset, conditions <- c(rep("N", 7), rep("S", 7))
    # input can be "all", "iqlr", "zero" or numeric for advanced users
    data(selex)
    #subset for efficiency
    selex <- selex[1201:1600,]
    conds <- c(rep("NS", 7), rep("S", 7))
    x <- aldex.clr(selex, conds, mc.samples=2, denom="all")
}
