\name{HighlightRegion-class}
\docType{class}
\alias{HighlightRegion-class}

\title{Class "HighlightRegion" is used to highlight vertical blocks of
  genomic regions.}
\description{HighlightRegion is used to highlight a genomic
  region of interest. The class offers the ability to highlight or block
  out regions of interest.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("HighlightRegion", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{start}:}{Object of class \code{"numeric"} genomic start position. }
    \item{\code{end}:}{Object of class \code{"numeric"} genomic end
      position. }
    \item{\code{region}:}{Object of class \code{"numericOrNull"} start
      and end number of the tracks to be covered by the region. These start
      from the first track (top = 1) to the last track:
      \code{length(gdObjects)} in the call to \code{gdObject}}
    \item{\code{coords}:}{Object of class \code{"character"} can be
      either "genomic" or "absolute", if the coordinates are "absolute"
      then one can plot things using the coordinate space defined by:
      lower-left (0,0) upper-right (1,1). In this case, start = x0, end
      = x1 and then region = (y0, y1). See the examples for more
      details. 
    }
    \item{\code{dp}:}{Object of class \code{"DisplayPars"} specifys the
      various display parameters.
    }
  }
}
\section{Extends}{
  Class \code{"\linkS4class{gdObject}"}, directly.
}
\section{Methods}{
  No methods defined with class "HighlightRegion" in the signature.
}
\author{James Bullard}
\examples{
if (interactive()) {
data("dummyData", package="GenomeGraphs")

ga <- new("GenomeAxis")
grF <- new("GeneRegion", start = 10000, end = 20000, chromosome = "I", strand = "+", biomart = yeastMart)
grR <- new("GeneRegion", start = 10000, end = 20000, chromosome = "I", strand = "-", biomart = yeastMart)
bt <- new("BaseTrack", base = yeastCons1[,1], value = yeastCons1[,2])
hr1 <- new("HighlightRegion", start = 11000, end = 13000,
           dp = DisplayPars(alpha = 1, color = "red", lty = "dashed", lwd = 3))
hr2 <- new("HighlightRegion", start = 15900, end = 16500)
          
gdPlot(list(grF, ga, grR, bt), highlightRegions = list(hr1, hr2))
}
}
\keyword{classes}
