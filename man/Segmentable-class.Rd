\name{Segmentable-class}
\docType{class}
\alias{Segmentable-class}
\alias{getSegmentation,Segmentable-method}

\title{Class "Segmentable"}
\description{The Segmentable class is used as an interface to determine
  whether or not a segmentation should be plotted.}
\section{Objects from the Class}{
  This object is mostly intended to be extended. 
}
\section{Slots}{
  \describe{
    \item{\code{segmentation}:}{Object of class
      \code{"SegmentationOrNULL"}}
  }
}
\section{Methods}{
  \describe{
    \item{getSegmentation}{\code{signature(obj = "Segmentable")}:
      returns the segmentation, this is essentially the method which is
      implemented by subclasses.}
  }
}
\seealso{
  \code{\linkS4class{Segmentation}}
}
\examples{
showClass("Segmentable")
}
\keyword{classes}
