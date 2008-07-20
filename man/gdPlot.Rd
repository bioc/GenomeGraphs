\name{gdPlot}
\alias{gdPlot}
\title{gdPlot is the main plotting function of the GenomeGraphs package}
\description{
  gdPlot is the main plotting function of the GenomeGraphs package.  A
  collection of objects are given as a list and these will then be plotten
  in the order given.
}
\usage{gdPlot(gdObjects, minBase = NA, maxBase = NA, highlightRegions = NULL, labelRot = 0)
}
\arguments{
  \item{gdObjects}{ This is either a list of \code{gdObjects} which will
    be plotted from top to bottom or a single \code{gdObjects} to be plotted.
  }
  \item{minBase}{ \code{minBase} defines the minimum base that will be
    plotted, if ommitted a minimum is determined from the objects in
    \code{gdObjects} if possible.
  }
  \item{maxBase}{\code{maxBase} defines the maximum base that will be
    plotted,if ommitted a minimum is determined from the objects in
    \code{gdObjects} if possible.
  }
  \item{highlightRegions}{\code{highlightRegions} defines a set of
    regions to overlay as rectangles on the plot. This argument is
    either a list or a single HighlightRegion object. 
  }
  \item{labelRot}{Rotate the track labels \code{labelRot} degrees. 
  }
}
\references{http://www.stat.berkeley.edu/~steffen/}
\author{Steffen Durinck and James Bullard}
\examples{

data("exampleData", package="GenomeGraphs")

minbase = min(probestart)
maxbase = max(probestart)

mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")

genesplus = new("GeneRegion", start = minbase, end = maxbase, strand = "+", chromosome = "3", biomart=mart)
genesmin = new("GeneRegion", start = minbase, end = maxbase, strand = "-", chromosome = "3", biomart=mart)

seg <- new("Segmentation", segments = segments, 
           segmentStart = segStart, segmentEnd = segEnd, dp = DisplayPars(color = "dodgerblue2", lwd=2,lty = "dashed"))

cop <- new("GenericArray", intensity  = cn, probeStart = probestart, 
           segmentation =  seg, dp = DisplayPars(size=3, color = "seagreen", type="dot"))

ideog = new("Ideogram", chromosome = "3")
expres = new("GenericArray", intensity = intensity, probeStart = exonProbePos, dp = DisplayPars(color="darkred", type="point"))
genomeAxis = new("GenomeAxis", add53 = TRUE, add35=TRUE)
gdPlot(list(ideog,expres,cop,genesplus,genomeAxis,genesmin), minBase = minbase, maxBase =maxbase)

}
\keyword{hplot}
