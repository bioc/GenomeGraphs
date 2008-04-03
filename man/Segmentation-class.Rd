\name{Segmentation-class}
\docType{class}
\alias{Segmentation-class}
\alias{getSegmentEnd,Segmentation-method}
\alias{getSegmentStart,Segmentation-method}
\alias{getSegments,Segmentation-method}

\title{Class "Segmentation" is used to specify segmentations to any
  class that extends Segmentable (GenericArray, BaseTrack)}
\description{A Segmentation object provides line segments to various gdObjects}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Segmentation",
  segments = list(1), segmentStart = list(1000), segmentEnd = list(1010))}.

}
\section{Slots}{
  \describe{
    \item{\code{segments}:}{Object of class \code{"list"} ~~ }
    \item{\code{segmentStart}:}{Object of class \code{"list"} ~~ }
    \item{\code{segmentEnd}:}{Object of class \code{"list"} ~~ }
    \item{\code{dp}:}{Object of class \code{"DisplayPars"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{gdObject}"}, directly.
}
\section{Methods}{
  \describe{
    \item{getSegmentEnd}{\code{signature(obj = "Segmentation")}: ... }
    \item{getSegmentStart}{\code{signature(obj = "Segmentation")}: ... }
    \item{getSegments}{\code{signature(obj = "Segmentation")}: ... }
  }
}
\author{James Bullard }
\examples{
showClass("Segmentation")
}
\keyword{classes}
