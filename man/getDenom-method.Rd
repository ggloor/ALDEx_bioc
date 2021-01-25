\name{getDenom}

\alias{getDenom}
\alias{getDenom,aldex.clr-method}

\title{getDenom}

\description{
  Returns the offset of the features used as the denominator
  as the basis for the log-ratio, for an \code{aldex.clr} object.
}

\usage{
  getDenom(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object.
  }
}

\details{
  Returns the offset of the features used as the denominator
  as the basis for the log-ratio. A vector of numbers is the offset
  of the non-0 features used in the denominator.
}

\value{
  A vector of integer row offsets.
}

\seealso{
  \code{aldex.clr}
}

\examples{
data(selex)
    #subset for efficiency
    selex <- selex[1201:1600,]
conds <- c(rep("NS", 7), rep("S", 7))
x <- aldex.clr(selex, conds, mc.samples = 2, denom = "iqlr", verbose = FALSE)
Denom <- getDenom(x)

# to find the names of housekeeping genes used
getFeatureNames(x)[getDenom(x)]
}

