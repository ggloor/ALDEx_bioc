\name{ALDEx2m-package}
\alias{ALDEx2m-package}
\alias{ALDEx2m}
\docType{package}
\title{
    Analysis of differential abundance taking sample variation into account
}
\description{
    A differential abundance analysis for the comparison of two or more conditions. For example, single-organism and meta-RNA-seq high-throughput sequencing assays, or of selected and unselected values from in-vitro sequence selections. Uses a Dirichlet-multinomial model to infer abundance from counts, that has been optimized for three or more experimental replicates. Infers sampling variation and calculates the expected false discovery rate given the biological and sampling variation using the Wilcox rank test or Welches t-test (aldex.ttest) or the glm and Kruskal Wallis tests (aldex.glm). Reports both P and fdr values calculated by the Benjamini Hochberg correction.
}
\references{
    Please use the citation given by \code{citation(package="ALDEx")}.
}
\keyword{ package }
\seealso{
    \code{\link{aldex.clr}},
    \code{\link{aldex.ttest}},
    \code{\link{aldex.glm}},
    \code{\link{aldex.effect}},
    \code{\link{selex}}
}
\examples{
    # see examples for the aldex.clr, aldex.ttest, aldex.effect, aldex.glm functions
}

