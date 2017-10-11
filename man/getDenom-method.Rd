\name{getDenom}

\alias{getDenom}
\alias{getDenom,aldex.clr-method}

\title{getDenom}

\description{
  Returns the denominator used as the basis for the log-ratio, for an \code{aldex.clr} object.
}

\usage{
getDenom(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the Monte Carlo Dirochlet instances
  derived from estimating the technical variance of the raw read count data, along with
  sample and feature information.
  }
}

\details{
  Returns the denominator used to calculate the log-ratios. "all" is the centred
  log-ratio. "iqlr" is the interquartile log-ratio. A vector of numbers is the offset
  of the variables used in the denominator
}

\value{
  A vector of values.
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

