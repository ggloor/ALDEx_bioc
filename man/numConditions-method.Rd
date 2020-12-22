\name{numConditions}

\alias{numConditions}
\alias{numConditions,aldex.clr-method}

\title{numConditions}

\description{
  Returns the number of conditions compared for analysis,
  for an \code{aldex.clr} object.
}

\usage{
numConditions(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object.
  }
}

\details{
  Returns the number of samples compared.
}

\value{
  A numeric representing the number of samples compared.
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
conditions <- numConditions(x)
}
