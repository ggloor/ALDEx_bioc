\name{getConditions}

\alias{getConditions}
\alias{getConditions,aldex.clr-method}

\title{getConditions}

\description{
  Returns the conditions, for
  \code{aldex.clr} object.
}

\usage{
  getConditions(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the Monte Carlo Dirichlet instances derived from estimating the technical variance of the raw read count data, along with sample and feature information.
  }
}

\details{
  Returns the conditions used.
}

\value{
  A matrix representing the conditions used.
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
scale.samps <- getConditions(x)
}
