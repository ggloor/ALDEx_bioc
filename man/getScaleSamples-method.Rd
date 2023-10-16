\name{getScaleSamples}

\alias{getScaleSamples}
\alias{getScaleSamples,aldex.clr-method}

\title{getScaleSamples}

\description{
  Returns the log2 scale samples if scale simulation is used, for
  \code{aldex.clr} object.
}

\usage{
  getScaleSamples(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the Monte Carlo Dirichlet instances derived from estimating the technical variance of the raw read count data, along with sample and feature information.
  }
}

\details{
  Returns NULL if scale simulation was not used or a matrix of log2 scale samples if scale simuation was used.
}

\value{
  A matrix representing the log2 scale samples if scale simulation was used.
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
scale.samps <- getScaleSamples(x)
}
