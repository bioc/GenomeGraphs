\name{ExonArray-class}
\docType{class}
\alias{ExonArray-class}
\alias{show,ExonArray-method}

\title{Class "ExonArray" representing probe level exon array data from Affymetrix}
\description{Represents probe level exon array data from
  Affymetrix. Makes it possible to visualize alternative splicing as
  measured by the Affymetrix exon array platform and relate it to known
  transcript isoforms annotated by Ensembl
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("ExonArray", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{intensity}:}{Object of class \code{"matrix"}, array data
      matrix containing probes  as the rows and samples as the columns}
    \item{\code{probeStart}:}{Object of class \code{"numeric"} vector with the start positions of the probes}
    \item{\code{probeEnd}:}{Object of class \code{"numeric"} vector with the end positions of the probes}
    \item{\code{probeId}:}{Object of class \code{"character"} vector containing the probeset identifiers}
    \item{\code{nProbes}:}{Object of class \code{"numeric"} vector defining how many probes there are for each exon/probeset}
    \item{\code{size}:}{Object of class \code{"numeric"} specifying the size of the ExonArray plot in the final plot}
    \item{\code{displayProbesets}:}{Object of class \code{"logical"} used to indicate if probe set names should be plotted or not}
    \item{\code{probesetSize}:}{Object of class \code{"numeric"} font size of the probeset identifiers to be plotted }
    \item{\code{color}:}{Object of class \code{"character"} vector of
      colors or one color that will be used to draw  the intensity values }
    \item{\code{mapColor}:}{Object of class \code{"character"}
      specifying the color of the lines that  map the probeset positions to the Ensembl gene annotation}
    \item{\code{lwd}:}{Object of class \code{"numeric"} vector of line widths to be used to plot the intensity data.}
    \item{\code{lty}:}{Object of class \code{"character"} vector of line types to be used to plot the intensity data.}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{gdObject}"}, directly.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "ExonArray")}: ... }
	 }
}
\references{http://www.stat.berkeley.edu/~steffen/}
\author{Steffen Durinck}
\seealso{
	objects to See Also as \code{\link{gdPlot}}}
\examples{
if(interactive()){
data("unrData", package="GenomeGraphs")
library(biomaRt)
mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")

title = new("Title", title ="ENSG00000009307", dp = DisplayPars(color = "darkslategray"))
exmapcol = rep("khaki", length(unrNProbes))
exmapcol[28]="darkred"
probeSetCol = rep("grey", length(unrNProbes))
probeSetCol[27:28]="darkslategray"
probeSetLwd = rep(1, length(unrNProbes))
probeSetLwd[27:28]=3

exon = new("ExonArray", intensity = unrData, probeStart = unrPositions[,3], probeEnd=unrPositions[,4], probeId = as.character(unrPositions[,1]), nProbes = unrNProbes, dp = DisplayPars(color = "blue", mapColor = exmapcol, probeSetColor = probeSetCol, probeSetLwd = probeSetLwd), displayProbesets=FALSE)
exon2 = new("ExonArray", intensity = unrData, probeStart = unrPositions[,3], probeEnd=unrPositions[,4], probeId = as.character(unrPositions[,1]), nProbes = unrNProbes, dp = DisplayPars(color = "blue", mapColor = exmapcol, probeSetColor = probeSetCol, probeSetLwd = probeSetLwd, plotMap=FALSE), displayProbesets=TRUE)

affyModel = new("GeneModel", exonStart = unrPositions[,3], exonEnd =  unrPositions[,4])
gene = new("Gene", id = "ENSG00000009307", biomart = mart)
transcript = new("Transcript", id ="ENSG00000009307" , biomart = mart)
legend = new("Legend", legend = c("affyModel","gene"), dp = DisplayPars(color= c("darkgreen","orange")))

gdPlot(list(title,exonarray1 = exon2,exonarray2= exon, AffymetrixModel= affyModel, gene, transcript, legend), minBase = min(exon@probeStart), maxBase=max(exon@probeEnd))
}
}
\keyword{classes}
