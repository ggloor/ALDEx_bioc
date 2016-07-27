\name{getFeatures}

\alias{getFeatures}
\alias{getFeatures,aldex.clr-method}

\title{getFeatures}

\description{
  Returns the features as a vector, for an \code{aldex.clr} object.
}

\usage{
getFeatures(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the Monte Carlo Dirochlet instances derived from estimating the technical variance of the raw read count data, along with sample and feature information.
  }
}

\details{
  Returns the features as a vector, for an \code{aldex.clr} object.
}

\value{
  A vector of features.
}

\seealso{
  \code{aldex.clr}
}

\examples{
data(selex)
    #subset for efficiency
    selex <- selex[1201:1600,]
conds <- c(rep("NS", 7), rep("S", 7))
x <- aldex.clr(selex, conds, mc.samples = 2, denom="all", verbose = FALSE)
features <- getFeatures(x)
}
