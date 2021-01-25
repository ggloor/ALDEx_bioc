\name{getDirichletInstances}

\alias{getDirichletInstances}
\alias{getDirichletInstances,aldex.clr-method}

\title{getDirichletInstances}

\description{
  Returns a list of the Monte Carlo Dirichlet
  instances created by the \code{aldex.clr} function.
}

\usage{
getDirichletInstances(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the
   Monte Carlo Dirichlet instances derived from estimating.
  }
}

\details{
  Returns a list of the raw Monte Carlo Dirichlet
  instances created by the \code{aldex.clr} function. These
  are probability estimates.
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
monteCarloDirInstances <- getDirichletInstances(x)
}
