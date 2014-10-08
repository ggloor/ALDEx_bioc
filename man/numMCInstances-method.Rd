\name{numMCInstances}

\alias{numMCInstances}
\alias{numMCInstances,aldex.clr-method}

\title{numMCInstances}

\description{
  Returns the number of Monte Carle Dirochlet instances generated for analysis, for an \code{aldex.clr} object.
}

\usage{
numMCInstances(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the Monte Carlo Dirochlet instances derived from estimating the technical variance of the raw read count data, along with sample and feature information.
  }
}

\details{
  Returns the number of Monte Carle Dirochlet instances generated for analysis.
}

\value{
  A numeric representing the number of Monte Carle Dirochlet instances generated for analysis.
}

\seealso{
  \code{aldex.clr}
}

\examples{
data(selex)
x <- aldex.clr(selex, mc.samples = 2, verbose = FALSE)
numInstances <- numMCInstances(x)
}
