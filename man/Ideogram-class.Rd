\name{Ideogram-class}
\docType{class}
\alias{Ideogram-class}
\alias{drawGD,Ideogram-method}

\title{Class "Ideogram", represent an Ideogram}
\description{An ideogram is a representation of a chromosome containing the banding pattern.  Note that currently ideograms are only available for hsapiens.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Ideogram", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{chromosome}:}{Object of class \code{"character"}, representing the chromosome that needs to be drawn. E.g. 3 if chromosome 3 needs to be drawn or Y for Y chromosome.}
    \item{\code{dp}:}{Object of class \code{"DisplayPars"}, can be used to specify the size  (default 1) of the ideogram in the final plot and to specify the highlighting color}
  }
}
\section{Methods}{
No methods defined with class "Ideogram" in the signature.
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
