\name{aldex.ttest}
\alias{aldex.ttest}
\title{
    calculate Welch's t-test and Wilcoxon test statistics
}
\description{
    calculates expected values of the Welch's t-test and Wilcoxon rank test on the data returned by clr_function.r
}
\usage{
aldex.ttest(clr, conditions, paired.test = FALSE, hist.plot=FALSE)
}
\arguments{
  \item{clr}{
    \code{clr} is the data output of the aldex.clr function
}
  \item{conditions}{
    a description of the data structure to be used for testing
}
  \item{paired.test}{
    whether the Welch's test should be paired or not
}
  \item{hist.plot}{
    whether to plot a histogram of P values for an individual Dirichlet Monte-Carlo instance. Plot is output to the standard R plotting device.
}
}
\details{
    An explicit example for two conditions is shown in the `Examples' below.
}
\value{
    Outputs a dataframe with the following information:
    \item{we.ep}{
        a vector containing the expected P value of the Welch's t-test for each feature
    }    
    \item{we.eBH}{
        a vector containing the expected value of the Benjamini Hochberg corrected P value for each feature    
    }
    \item{wi.ep}{
        a vector containing the expected P value of the Wilcoxon test for each feature    
    }
    \item{wi.eBH}{
        a vector containing the expected value of the Benjamini Hochberg corrected P value for each feature
    }
}
\references{
    Please use the citation given by \code{citation(package="ALDEx")}.
}
\author{
    Greg Gloor
}

\seealso{
    \code{\link{aldex.clr}},
    \code{\link{aldex.glm}},
    \code{\link{aldex.effect}},
    \code{\link{selex}}
}
\examples{
    # x is the output of the \code{x <- aldex.clr(data, mc.samples)} function
    # conditions is a description of the data
    # for the selex dataset, conditions <- c(rep("N", 7), rep("S", 7))
    data(selex)
    #subset for efficiency
    selex <- selex[1201:1600,]
    conds <- c(rep("NS", 7), rep("S", 7))
    x <- aldex.clr(selex, conds, mc.samples=2, denom="all")
    ttest.test <- aldex.ttest(x, conds)

}
