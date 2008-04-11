\name{Legend-class}
\docType{class}
\alias{Legend-class}
\alias{drawGD,Legend-method}

\title{Class "Legend", represents a legend to add to a plot}
\description{This class represents a legend}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Legend", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{legend}:}{Object of class \code{"character"}, vector with names of the items in the legend}
    \item{\code{dp}:}{Object of class \code{"DisplayPars"} size of the legend (size), the size of the font (cex) and the colors (color) of the legend}
  }
}
\section{Methods}{
No methods defined with class "Legend" in the signature.
}
\references{http://www.stat.berkeley.edu/~steffen/ }
\author{Steffen Durinck}
\seealso{
	objects to See Also as \code{\link{gdPlot}}
}
\examples{
showClass("Legend")
}
\keyword{classes}
