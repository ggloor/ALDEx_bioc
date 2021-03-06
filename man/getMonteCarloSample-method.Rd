\name{getMonteCarloSample}

\alias{getMonteCarloSample}
\alias{getMonteCarloSample,aldex.clr,numeric-method}

\title{getMonteCarloSample}

\description{
  Returns a single Monte Carlo Dirichlet instance for all samples
  for an \code{aldex.clr} object.
}

\usage{
getMonteCarloSample(.object,i)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the log-ratio
  transformed Monte Carlo Dirichlet instances derived from
  estimating the technical variance of the raw read count data,
  along with sample and feature information.
  }
  \item{i}{The numeric index of the desired Monte-Carlo instance.
  }
}

\details{
  Returns the designated Monte Carlo Dirichlet instance for all
  samples generated from analysis.
}

\value{
  A matrix representing the designated log-ratio transformed
  Monte Carlo Dirichlet instance for all samples.
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
