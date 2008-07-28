\name{TextOverlay-class}
\docType{class}
\alias{TextOverlay-class}
\alias{drawOverlay, TextOverlay-method}

\title{Class "TextOverlay"}
\description{Textual overlay classes}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{\link{makeTextOverlay}}
}
\section{Slots}{
  \describe{
    \item{\code{text}:}{Object of class \code{"character"}}
    \item{\code{xpos}:}{Object of class \code{"numeric"}}
    \item{\code{ypos}:}{Object of class \code{"numeric"}}
    \item{\code{region}:}{Object of class \code{"numericOrNull"}}
    \item{\code{coords}:}{Object of class \code{"character"}}
    \item{\code{dp}:}{Object of class \code{"DisplayPars"}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{Overlay}"}, directly.
Class \code{"\linkS4class{gdObject}"}, by class "Overlay", distance 2.
}
\section{Methods}{
  \describe{
    \item{drawOverlay}{\code{signature(obj = "TextOverlay")}: ... }
  }
}
\examples{
showClass("TextOverlay")
}
\keyword{classes}
