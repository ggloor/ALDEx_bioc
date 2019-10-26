\name{getMonteCarloSample}

\alias{getMonteCarloSample}
\alias{getMonteCarloSample,aldex.clr,numeric-method}

\title{getMonteCarloSample}

\description{
  Returns a set of Monte Carlo Dirochlet replicates for one sample generated from analysis,
  for an \code{aldex.clr} object.
}

\usage{
getMonteCarloSample(.object,i)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the Monte Carlo Dirochlet
  instances derived from estimating the technical variance of the raw read count data,
  along with sample and feature information.
  }
  \item{i}{The numeric index of the desired sample.
  }
}

\details{
  Returns the designated Monte Carlo Dirochlet replicates for one sample generated from analysis.
}

\value{
  A data frame representing the designated Monte Carlo Dirochlet replicates for one sample generated
  from analysis.
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
monteCarloSample <- getMonteCarloSample(x,1)
}
