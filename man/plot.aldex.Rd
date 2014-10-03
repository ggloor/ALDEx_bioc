\name{aldex.plot}
\alias{aldex.plot}
\title{
    Plot an \code{aldex} Object
}
\description{
    Create `MW'- or `MA'-type plots from the given \code{aldex} object.
}
\usage{
    \method{aldex}{plot}(  x, ..., type=c("MW","MA"),
    xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,
    all.col=rgb(0,0,0,0.2), all.pch=19, all.cex=0.4,
    called.col=red, called.pch=20, called.cex=0.6,
    thres.line.col=darkgrey, thres.lwd=1.5,
    test=welch, cutoff=0.1, rare.col=black, rare=0,
    rare.pch=20, rare.cex=0.2 )
}
\arguments{
    \item{x}{
        an object of class \code{aldex} produced by the \code{aldex} function
    }
    \item{...}{
        optional, unused arguments included for compatibility withe the S3 method signature
    }
    \item{type}{
        which type of plot is to be produced. MA is a Bland-Altman style plot; MW is a difference between to a variance within plot as described in the paper
    }
    \item{test}{
        the method of calculating significance, one of:
        welch = welch's t test;
        wilcox = wilcox rank test;
        glm = glm;
        kruskal = Kruskal-Wallace test
    }
    \item{cutoff}{
        the Benjamini-Hochberg fdr cutoff, default 0.1
    }
    \item{xlab}{
        the x-label for the plot, as per the parent \code{plot} function
    }
    \item{ylab}{
        the y-label for the plot, as per the parent \code{plot} function
    }
    \item{xlim}{
        the x-limits for the plot, as per the parent \code{plot} function
    }
    \item{ylim}{
        the y-limits for the plot, as per the parent \code{plot} function
    }
    \item{all.col}{
        the default colour of the plotted points
    }
    \item{all.pch}{
        the default plotting symbol
    }
    \item{all.cex}{
        the default symbol size
    }
    \item{called.col}{
        the colour of points with false discovery rate, q <= 0.1
    }
    \item{called.pch}{
        the symbol of points with false discovery rate, q <= 0.1
    }
    \item{called.cex}{
        the character expansion of points with false discovery rate, q <= 0.1
    }
    \item{thres.line.col}{
        the colour of the threshold line where within and between group variation is equivalent
    }
    \item{thres.lwd}{
        the width of the threshold line where within and between group variation is equivalent
    }
    \item{rare}{
        relative abundance cutoff for rare features, default 0 or the mean abundance 
    }
    \item{rare.col}{
        color for rare features, default black
    }
    \item{rare.pch}{
        the default symbol of rare features
    }
    \item{rare.cex}{
        the default symbol size of rare points
    }

}
\details{
    This particular specialization of the \code{plot} function is relatively simple and provided for convenience. For more advanced control of the plot is is best to use the values returned by \code{summary(x)}.
}
\value{
    None.
}
\references{
    Please use the citation given by \code{citation(package="ALDEx")}.
}
\seealso{
    \code{\link{aldex}},
    \code{\link{aldex.effect}},
    \code{\link{aldex.ttest}},
    \code{\link{aldex.glm}}
}
\examples{
    # See the examples for 'aldex'.
}
