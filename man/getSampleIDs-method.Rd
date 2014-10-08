\name{getSampleIDs}

\alias{getSampleIDs}
\alias{getSampleIDs,aldex.clr-method}

\title{getSampleIDs}

\description{
  Returns the names of the samples for an \code{aldex.clr} object. These can be used to access the original reads, as in \code{reads\$sampleID} (if the reads are a data frame).
}

\usage{
getSampleIDs(.object)
}

\arguments{
  \item{.object}{A \code{aldex.clr} object containing the Monte Carlo Dirochlet instances derived from estimating the technical variance of the raw read count data, along with sample and feature information.
  }
}

\details{
  Returns the names of the samples. These can be used to access the original reads, as in \code{reads\$sampleID} (if the reads are a data frame).
}

\value{
  A vector of sample names.
}

\seealso{
  \code{aldex.clr}
}

\examples{
data(selex)
x <- aldex.clr(selex, mc.samples = 2, verbose = FALSE)
sampleIDs <- getSampleIDs(x)
}
