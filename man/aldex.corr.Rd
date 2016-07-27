\name{aldex.corr}
\alias{aldex.corr}
\title{
    calculate Pearson's Product moment and Spearman's rank correlations
}
\description{
    calculates expected values of Pearson's Product moment and Spearman's rank
    correlations on the data returned by aldex.clr. NOTE: this function will
    be replaced by a compositionally correct method in the next release cycle.
}
\usage{
aldex.corr(clr, covar)
}
\arguments{
    \item{clr}{
      \code{clr} is the data output of the aldex.clr function
}
    \item{covar}{
    a per-sample continuous variable to be correlated with the clr values
}

}
\details{
    An explicit example for two conditions is shown in the 'Examples' below.
}
\value{
    Outputs a dataframe with the following information:
      \item{pearson.ecor}{
        a vector containing the expected Pearson's Product moment value for each feature
    }
      \item{pearson.ep}{
        a vector containing the expected P value of the Pearson Product moment value for
        each feature
    }
      \item{pearson.eBH}{
        a vector containing the expected Benjamini-Hochberg corrected P value of the
        Pearson Product moment value for each feature
    }
      \item{spearman.erho}{
        a vector containing the expected Spearman's rank correlation value for each feature
    }
      \item{spearman.ep}{
        a vector containing the expected P value of Spearman's rank correlation value for
        each feature
    }
      \item{spearman.eBH}{
        a vector containing the expected Benjamini-Hochberg corrected P value of
        Spearman's rank correlation value for each feature
    }
}
\references{
    Please use the citation given by \code{citation(package="ALDEx")}.
}
\author{
    Arianne Albert
}

\seealso{
      \code{\link{aldex.clr}},
      \code{\link{aldex.glm}},
      \code{\link{aldex.effect}},
      \code{\link{selex}}
}
\examples{
    # x is the output of the \code{x <- aldex.clr(data, conds, mc.samples, denom="all", useMC)} function
    # conditions is a description of the data
    # aldex.ttest(clr, covar)
}
