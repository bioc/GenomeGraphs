\name{Smoothing-class}
\Rdversion{1.1}
\docType{class}
\alias{Smoothing-class}

\title{Class "Smoothing"}
\description{Simple object to overlay line segments specified by x,y
  coordinates.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{makeSmoothing}.
}
\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"numeric"}}
    \item{\code{y}:}{Object of class \code{"numeric"}}
    \item{\code{dp}:}{Object of class \code{"DisplayPars"}}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{TrackOverlay}"}, directly.
  Class \code{"\linkS4class{gdObject}"}, by class "TrackOverlay", distance 2.
  Class \code{"\linkS4class{TrackOverlayOrNull}"}, by class "TrackOverlay", distance 2.
}
\section{Methods}{
  No methods defined with class "Smoothing" in the signature.
}
\examples{
showClass("Smoothing")
}
\keyword{classes}
