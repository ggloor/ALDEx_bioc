\name{getFeatureNames}

\alias{getFeatureNames}
\alias{getFeatureNames,aldex.clr-method}

\title{getFeatureNames}

\description{
  Returns the names of the features as a vector, for an \code{aldex.clr} object.
}

\usage{
getFeatureNames(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object.
  }
}

\details{
  Returns the names of the \code{keys} that can be used to
  subset the data rows. The \code{keys} values are the rsid's.
}

\value{
  A vector of feature names.
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
featureNames <- getFeatureNames(x)
}
