\name{BaseTrack-class}
\docType{class}
\alias{BaseTrack-class}
\alias{show,BaseTrack-method}

\title{Class "BaseTrack" represents base specific data}
\description{Represents specific data, e.g. how many times was every base sequenced}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("BaseTrack", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{base}:}{Object of class \code{"numeric"}. Is a vector of base positions }
    \item{\code{value}:}{Object of class \code{"numeric"}.  Is a vector of corresponding values for every base}
    \item{\code{strand}:}{Object of class \code{"character"} represents that DNA strand}
    \item{\code{size}:}{Object of class \code{"numeric"} specifying the size of the SequenceFreqBase plot in the final plot}
    \item{\code{color}:}{Object of class \code{"character"} specifying the color to use for plotting}
    \item{\code{lty}:}{Object of class \code{"character"} specifying the line type }
    \item{\code{lwd}:}{Object of class \code{"numeric"} specifying the line width}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "BaseTrack")}: ... }
	 }
}
\references{http://www.stat.berkeley.edu/~steffen/}
\author{Steffen Durinck}

\seealso{
	objects to See Also as \code{\link{gdPlot}}
}
\examples{
showClass("BaseTrack")
}
\keyword{classes}
