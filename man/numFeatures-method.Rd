\name{numFeatures}

\alias{numFeatures}
\alias{numFeatures,aldex.clr-method}

\title{numFeatures}

\description{
  Returns the number of non-0 features associated with the data,
  for an \code{aldex.clr} object.
}

\usage{
  numFeatures(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object.
  }
}

\details{
  Returns the number of features associated with the data that are
  not 0 in all samples.
}

\value{
  A numeric representing the number of non-0 features associated
  with the data.
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
numFeatures <- numFeatures(x)
}
