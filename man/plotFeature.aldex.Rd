\name{aldex.plotFeature}
\alias{aldex.plotFeature}
\title{
    Generate dispersion, difference, and effect density plots for a feature, given an \code{aldex} Object
}
\description{
    Calculates dispersion, difference, and effect vectors from the given \code{aldex} object for the selected feature. Generates plots for each group's pooled dispersion vector, the element-wise absolute difference between groups, and the effect size of the difference between groups.
}
\usage{
    \method{aldex}{plotFeature}(clrData, featureName, pooledOnly=FALSE,
                                  densityOnly=FALSE)
}
\arguments{
    \item{clrData}{
        an object of class \code{aldex} produced by the \code{aldex} function
    }
    \item{featureName}{
        the name of the feature to generate plots for, provided as a string
    }
    \item{pooledOnly}{
        a logical value to generate only the pooled dispersion vector plot
    }
    \item{densityOnly}{
        a logical value to toggle the boxplots overlaid on each plot
    }
}
\details{
    The dispersion, difference, and effect size vectors are defined in (Fernandes et al. (2018). A reproducible effect size is more useful than an irreproducible hypothesis test to analyze high-throughput sequencing datasets. Preprint at https://arxiv.org/abs/1809.02623).
}
\value{
    None.
}
\references{
    Please use the citation given by \code{citation(package="ALDEx")}.
}
\seealso{
    \code{\link{aldex}},
    \code{\link{aldex.clr}}
}
\examples{
    data(selex)
        #subset for efficiency
        selex <- selex[1201:1600,]
    conds <- c(rep("NS", 7), rep("S", 7))
    x <- aldex.clr(selex, conds, mc.samples = 128, denom = "all", verbose = FALSE)
    aldex.plotFeature(x,"T:E:K:D")
}
