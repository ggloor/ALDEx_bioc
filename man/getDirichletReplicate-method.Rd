\name{getDirichletReplicate}

\alias{getDirichletReplicate}
\alias{getDirichletReplicate,aldex.clr,numeric-method}

\title{getDirichletReplicate}

\description{
  Returns the raw per-sample Monte Carlo Dirichlet
  replicates generated from analysis, for an \code{aldex.clr}
  object.
}

\usage{
  getDirichletReplicate(.object,i)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the Monte Carlo Dirichlet
  replicates derived from estimating the technical variance of the raw read count data,
  along with sample and feature information.
  }
  \item{i}{The numeric index of the desired sample replicate.
  }
}

\details{
  Returns the raw per-sample Monte Carlo Dirichlet replicates. These
  are estimated probabilities.
}

\value{
  A numeric matrix.
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
DirichletReplicate <- getDirichletReplicate(x,1)
}
