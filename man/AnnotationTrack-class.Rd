\name{AnnotationTrack-class}
\docType{class}
\alias{AnnotationTrack-class}
\alias{drawGD,AnnotationTrack-method}
\alias{getPlotId,AnnotationTrack-method}
\alias{initialize,AnnotationTrack-method}

\title{Class "AnnotationTrack"}
\description{A generic object to store annotation}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AnnotationTrack", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{chr}:}{Object of class \code{"numeric"}}
    \item{\code{strand}:}{Object of class \code{"numeric"}}
    \item{\code{regions}:}{Object of class \code{"dfOrNULL"}}
    \item{\code{dp}:}{Object of class \code{"DisplayPars"}}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{gdObject}"}, directly.
}
\section{Methods}{
  \describe{
    \item{drawGD}{\code{signature(gdObject = "AnnotationTrack")}: ... }
    \item{getPlotId}{\code{signature(obj = "AnnotationTrack")}: ... }
    \item{initialize}{\code{signature(.Object = "AnnotationTrack")}: ... }
  }
}
\author{James Bullard}
\examples{
showClass("AnnotationTrack")
}
\keyword{classes}
