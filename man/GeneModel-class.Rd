\name{GeneModel-class}
\docType{class}
\alias{GeneModel-class}
\alias{drawGD,GeneModel-method}

\title{Class "GeneModel", represents a custom gene model}
\description{This class represents a custom gene model defined by exon boundaries.  An example of this class could be an Affymetrix gene model used to create the Affy Exon array}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GeneModel", ...)}. 
}
\section{Slots}{
	 \describe{
    \item{\code{exonStart}:}{Object of class \code{"numeric"}, vector containing the start positions of the exons that are to be drawn}
    \item{\code{exonEnd}:}{Object of class \code{"numeric"}, vector containing the end positions of the exons that are to be drawn }
    \item{\code{chromosome}:}{Object of class \code{"numeric"} , chromosome name}
    \item{\code{dp}:}{Object of class \code{"DisplayPars"}, color of the exons and size of the exon model in the final plot }
  }
}
\section{Methods}{
No methods defined with class "GeneModel" in the signature.
}
\references{http://www.stat.berkeley.edu/~steffen/ }
\author{Steffen Durinck}
\seealso{
	objects to See Also as \code{\link{gdPlot}}}
\examples{
data("unrData", package="GenomeGraphs")
affyModel = new("GeneModel", exonStart = unrPositions[,3], exonEnd =  unrPositions[,4])
gdPlot(list(affyModel), minBase = min(unrPositions[,3]), maxBase=max(unrPositions[,4]))


}
\keyword{classes}
