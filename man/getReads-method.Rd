\name{getReads}

\alias{getReads}
\alias{getReads,aldex.clr-method}

\title{getReads}

\description{
  Returns the count table used as input for analysis, for an \code{aldex.clr} object.
}

\usage{
getReads(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the Monte Carlo Dirochlet instances derived from estimating the technical variance of the raw read count data, along with sample and feature information.
  }
}

\details{
  Returns the count table.
}

\value{
  A data frame representing the count table used as input for analysis.
}

\seealso{
  \code{aldex.clr}
}

\examples{
data(selex)
    #subset for efficiency
    selex <- selex[1201:1600,]
conds <- c(rep("NS", 7), rep("S", 7))
x <- aldex.clr(selex, conds, mc.samples = 2, denom = "all", verbose = FALSE)
reads <- getReads(x)
}
