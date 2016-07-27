\name{numFeatures}

\alias{numFeatures}
\alias{numFeatures,aldex.clr-method}

\title{numFeatures}

\description{
  Returns the number of features associated with the data, for an \code{aldex.clr} object.
}

\usage{
numFeatures(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the Monte Carlo Dirochlet instances derived from estimating the technical variance of the raw read count data, along with sample and feature information.
  }
}

\details{
  Returns the number of features associated with the data.
}

\value{
  A numeric representing the number of features associated with the data.
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
numFeatures <- numFeatures(x)
}
