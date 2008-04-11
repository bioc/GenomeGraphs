\name{GenericArray-class}
\docType{class}
\alias{GenericArray-class}
\alias{show,GenericArray-method}
\alias{drawGD,GenericArray-method}

\title{Class "GenericArray", representing array data }
\description{The Generic Array class is a class that can be used to
  create plots from array data such as microarrays and arrayCGH
  platforms. It can represent,  the data as line plots or dot plots and
  segments can be included as well
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GenericArray", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{intensity}:}{Object of class \code{"matrix"}, matrix containing the intensities of expression or cgh data.  Rows should be probes, columns samples}
    \item{\code{probeStart}:}{Object of class \code{"numeric"}, start position of the probes}
    \item{\code{probeEnd}:}{Object of class \code{"numeric"}, end position of the probes if available}
    \item{\code{type}:}{Object of class \code{"character"}, has two values: line and point.  If line is selected then a line will be plotted for each sample.  If point is selected. points will be plotted.}
    \item{\code{segments}:}{Object of class \code{"list"}, if segments are available they have to be represented as a list, each position in the list should correspond to each sample as defined by the columns of the intensity matrix }
    \item{\code{segmentStart}:}{Object of class \code{"list"}, containing the start positions of the segments}
    \item{\code{segmentEnd}:}{Object of class \code{"list"}, containing the end positions of the segments }
    \item{\code{color}:}{Object of class \code{"character"}, represents the color to be used to plot the intensity matrix, can be a vector of colors for multiple samples }
    \item{\code{lty}:}{Object of class \code{"character"},if line is selected as type, lty specifies which line type should be used.  Can be a vector of line types for multiple samples}
    \item{\code{pch}:}{Object of class \code{"numeric"}, if point is selected as type then pch represents which symbol should be used to plot the points}
    \item{\code{pointSize}:}{Object of class \code{"numeric"}, specifies the point size if point is selected as type}
    \item{\code{lwd}:}{Object of class \code{"numeric"}, specifies the line width if line is selected as type.  Can be a vector of line widths if multple samples are present}
    \item{\code{size}:}{Object of class \code{"numeric"} representing the size of the Array plot in the final plot}
    \item{\code{segmentColor}:}{Object of class \code{"character"}. If segments are to be plotted, this will define the color of the segment lines }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "GenericArray")}: ... }
	 }
}
\references{http://www.stat.berkeley.edu/~steffen/}
\author{Steffen Durinck}
\seealso{
	objects to See Also as \code{\link{gdPlot}}
}
\examples{
if(interactive()){
data("dummyData", package="GenomeGraphs")

minbase <- 180292097 
maxbase <- 180492096
ideog <- new("Ideogram", chromosome = "3")
expres <- new("GenericArray", intensity = intensity, probeStart = exonProbePos, 
              dp = DisplayPars(color="darkred", type="point"))
gdPlot(list(ideog, expres), minBase = minbase, maxBase =maxbase)
}
}
\keyword{classes}
