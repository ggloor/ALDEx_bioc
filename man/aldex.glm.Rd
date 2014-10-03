\name{aldex.glm}
\alias{aldex.glm}
\title{
   calculate glm and Kruskal Wallis test statistics
}
\description{
    calculates expected values of the glm and Kruskal Wallis functions on the data returned by clr_function.r
}
\usage{
aldex.glm(clr, conditions, useMC=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clr}{
    \code{clr} is the data output of \code{aldex.clr}
}
  \item{conditions}{
    a description of the data structure to be used for testing
}
  \item{useMC}{
    use multicore by default (FALSE)
}
}
\details{
    An explicit example for two conditions is shown in the `Examples' below.
}
\value{
    Outputs a dataframe with the following information:
    \item{kw.ep}{
        a vector containing the expected P value of the Kruskal Wallis test for each feature    
    }
    \item{kw.eBH}{
        a vector containing the expected value of the Benjamini Hochberg corrected P value for each feature    
    }
    \item{glm.ep}{
        a vector containing the expected P value of the glm test for each feature    
    }
    \item{glm.eBH}{
        a vector containing the expected value of the Benjamini Hochberg corrected P value for each feature
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
    \code{\link{aldex.ttest}},
    \code{\link{aldex.effect}},
    \code{\link{selex}}
}
\examples{
    # x is the output of the \code{x <- aldex.clr(data, mc.samples)} function
    # conditions is a description of the data
    # for the selex dataset, conditions <- c(rep("N", 7), rep("S", 7))
    data(selex)
    x <- aldex.clr(selex, mc.samples=1)
    conditions <- c(rep("N", 7), rep("S", 7))
    glm.test <- aldex.glm(x, conditions)
}
