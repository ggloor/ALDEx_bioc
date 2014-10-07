\name{getMonteCarloInstances}

\alias{getMonteCarloInstances}
\alias{getMonteCarloInstances,aldex.clr-method}

\title{getMonteCarloInstances}

\description{
  Returns the Monte Carlo Dirochlet instances used to create an \code{aldex.clr} object.
}

\usage{
getMonteCarloInstances(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the Monte Carlo Dirochlet instances derived from estimating the technical variance of the raw read count data, along with sample and feature information.
  }
}

\details{
  Returns the Monte Carlo Dirochlet instances used to create an \code{aldex.clr} object.
}

\value{
  A list of data frames of Monte Carlo Dirochlet instances derived from estimating the technical variance of the raw read count data.
}

\seealso{
  \code{aldex.clr}
}

\examples{
data(selex)
x <- aldex.clr(selex, mc.samples = 2, verbose = FALSE)
monteCarloInstances <- getMonteCarloInstances(x)
}
