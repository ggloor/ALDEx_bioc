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
x <- aldex.clr(selex, mc.samples = 2, verbose = FALSE)
features <- getFeatures(x)
}
