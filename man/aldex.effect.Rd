\name{aldex.effect}
\alias{aldex.effect}
\title{calculate effect sizes and differences between conditions}
\description{
    determines the median clr abundance of the feature in all samples and in groups
    determines the median difference between the two groups
    determines the median variation within each two group
    determines the effect size, which is the median of the ratio of the between group difference and the larger of the variance within groups
}
\usage{
aldex.effect(clr, conditions, verbose = TRUE, include.sample.summary = FALSE, useMC=FALSE)
}
\arguments{
    \item{clr}{
        \code{clr} is the data output of \code{aldex.clr}
    }
    \item{conditions}{
        a description of the data structure to be used for testing
    }
    \item{verbose}{
        Print diagnostic information while running. Useful only for debugging if fails on large datasets
    }
    \item{include.sample.summary}{
        include median clr values for each sample, defaults to FALSE
    }
      \item{useMC}{
    use multicore by default (FALSE)
}
}
\details{
    An explicit example for two conditions is shown in the `Examples' below.
}
\value{
    returns a dataframe with the following information:
    \item{rab.all}{
        a vector containing the median clr value for each feature
    }
    \item{rab.win.conditionA}{
        a vector containing the median clr value for each feature in condition A
    }
    \item{rab.win.conditionB}{
        a vector containing the median clr value for each feature in condition B
    }
    \item{diff.btw}{
        a vector containing the per-feature median difference between condition A and B
    }
    \item{diff.btw}{
        a vector containing the per-feature maximum median difference between Dirichlet instances within conditions
    }
    \item{effect}{
        a vector containing the per-feature effect size
    }
    \item{overlap}{
        a vector containing the per-feature proportion of effect size that is 0 or less
    }
}
\references{
    Please use the citation given by \code{citation(package="ALDEx")}.
}
\author{
    Greg Gloor, Andrew Fernandes, Matt Links
}
\seealso{
    \code{\link{aldex.clr}},
    \code{\link{aldex.ttest}},
    \code{\link{aldex.glm}},
    \code{\link{selex}}
}
\examples{
    # x is the output of the \code{x <- clr(data, mc.samples)} function
    # conditions is a description of the data
    # for the selex dataset, conditions <- c(rep("N", 7), rep("S", 7))
    data(selex)
    x <- aldex.clr(selex, mc.samples=2)
    conditions <- c(rep("N", 7), rep("S", 7))
    effect.test <- aldex.effect(x, conditions)
}
