\name{Title-class}
\docType{class}
\alias{Title-class}
\alias{drawGD,Title-method}

\title{Class "Title" representing the title of a plot }
\description{Represent the title of a plot}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Title", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{title}:}{Object of class \code{"character"} which will be used as title }
    \item{\code{dp}:}{Object of class \code{"DisplayPars"} specifying the size and color of the title in the final plot }
  }
}
\section{Methods}{
No methods defined with class "Title" in the signature.
}
\references{http://www.stat.berkeley.edu/~steffen/ }
\author{Steffen Durinck}
\seealso{
	objects to See Also as \code{\link{gdPlot}}
}
\examples{
showClass("Title")
}
\keyword{classes}
