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
    \item{\code{strand}:}{Object of class \code{"character"} represents
      that DNA strand}
    \item{\code{dp}:}{Object of class DisplayPars to control various
      features of how the track is displayed.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{gdObject}"}, directly.
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
if (interactive()) {
data("dummyData", package="GenomeGraphs")
ga <- new("GenomeAxis")
bt <- new("BaseTrack", base = yeastCons1[,1], value = yeastCons1[,2],
          dp = DisplayPars(color = "darkblue"))
gdPlot(list(ga, bt))
}
}
\keyword{classes}
