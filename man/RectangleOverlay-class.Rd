\name{RectangleOverlay-class}
\docType{class}
\alias{RectangleOverlay-class}
\alias{drawOverlay,RectangleOverlay-method}

\title{Class "RectangleOverlay"}
\description{Rectangular Overlay}
\section{Objects from the Class}{
  Objects can be created by calls of the form makeRectangleOverlay (\code{\link{makeRectangleOverlay}}).
}
\section{Slots}{
  \describe{
    \item{\code{start}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{end}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{region}:}{Object of class \code{"numericOrNull"} ~~ }
    \item{\code{coords}:}{Object of class \code{"character"} ~~ }
    \item{\code{dp}:}{Object of class \code{"DisplayPars"} ~~ }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{Overlay}"}, directly.
  Class \code{"\linkS4class{gdObject}"}, by class "Overlay", distance 2.
}
\section{Methods}{
  \describe{
    \item{drawOverlay}{\code{signature(obj = "RectangleOverlay")}: ... }
  }
}
\examples{
showClass("RectangleOverlay")
}
\keyword{classes}
