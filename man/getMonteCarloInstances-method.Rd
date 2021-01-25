\name{getMonteCarloInstances}

\alias{getMonteCarloInstances}
\alias{getMonteCarloInstances,aldex.clr-method}

\title{getMonteCarloInstances}

\description{
  Returns a list of the log-ratio transformed Monte Carlo Dirichlet
  instances created by the \code{aldex.clr} function.
}

\usage{
getMonteCarloInstances(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the
   clr-transformed Monte Carlo Dirichlet instances derived from
   estimating the technical variance of the raw read count data,
    along with sample and feature information.
  }
}

\details{
  Returns a list of the log-ratio transformed Monte Carlo Dirichlet
  instances created by the \code{aldex.clr} function.
}

\value{
  A list of data frames.
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
monteCarloInstances <- getMonteCarloInstances(x)
}
