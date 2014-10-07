\name{getFeatureNames}

\alias{getFeatureNames}
\alias{getFeatureNames,aldex.clr-method}

\title{getFeatureNames}

\description{
  Returns the names of the features as a vector, for an \code{aldex.clr} object.
}

\usage{
getFeatureNames(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the Monte Carlo Dirochlet instances derived from estimating the technical variance of the raw read count data, along with sample and feature information.
  }
}

\details{
  Returns the names of the \code{keys} that can be used to subset the data rows. The \code{keys} values are the rsid's.
}

\value{
  A vector of feature names.
}

\seealso{
  \code{aldex.clr}
}

\examples{
data(selex)
x <- aldex.clr(selex, mc.samples = 2, verbose = FALSE)
featureNames <- getFeatureNames(x)
}
