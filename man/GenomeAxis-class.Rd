\name{GenomeAxis-class}
\docType{class}
\alias{GenomeAxis-class}
\alias{drawGD,GenomeAxis-method}

\title{Class "GenomeAxis", representing a genomic coordinate axis}
\description{Represents a genomic coordinate axis}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GenomeAxis", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{add53}:}{Object of class \code{"logical"}, indicating if 5 to 3 prime direction needs to be plotted}
    \item{\code{add35}:}{Object of class \code{"logical"}, indicating if
      3 to 5 prime direction needs to be plotted }
    \item{\code{dp}:}{Object of class \code{"DisplayPars"}, containing the display parameters such as size of the plot and color}
    \item{\code{littleTicks}:}{Object of class \code{"logical"},
      indicating if the genome axis should be dense for improved
      locating of regions of interest.}
  }
}
\section{Methods}{
No methods defined with class "GenomeAxis" in the signature.
}
\references{http://www.stat.berkeley.edu/~steffen/ }
\author{Steffen Durinck}
\seealso{
	objects to See Also as \code{\link{gdPlot}}
}
\examples{
if(interactive()){
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")
genomeAxis = new("GenomeAxis", add53=TRUE)
plusStrand = new("GeneRegion", chromosome = "17", start = 30450000, end = 30550000, strand = "+", biomart = mart)
gdPlot(list(genomeAxis, plusStrand), minBase = 30450000, maxBase =  30550000)
}
}
\keyword{classes}
