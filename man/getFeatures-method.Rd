\name{getFeatures}

\alias{getFeatures}
\alias{getFeatures,aldex.clr-method}

\title{getFeatures}

\description{
  Returns the features as a vector, for an \code{aldex.clr} object.
}

\usage{
getFeatures(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object.
  }
}

\details{
  Returns the features from the first sample and first Monte-Carlo
  replicate as a vector, for an \code{aldex.clr} object. Used only
  for troubleshooting purposes.
}

\value{
  A vector of features.
}

\seealso{
  \code{aldex.clr}
}

\examples{
data(selex)
    #subset for efficiency
    selex <- selex[1201:1600,]
conds <- c(rep("NS", 7), rep("S", 7))
x <- aldex.clr(selex, conds, mc.samples = 2, denom="all", verbose = FALSE)
features <- getFeatures(x)
}
