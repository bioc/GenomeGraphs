\name{makeSmoothing}
\alias{makeSmoothing}
\title{Create objects of class Smoothing}
\description{
  Construct objects of class Smoothing
}
\usage{
makeSmoothing(x, y, dp = NULL)
}
\arguments{
  \item{x}{x-coordinate}
  \item{y}{y-coordinate}
  \item{dp}{ The Display parameters.}
}
\value{
  An object of class Smoothing
}
\examples{
data("exampleData", package="GenomeGraphs")
seg <- makeSmoothing(probestart, lowess(probestart, cn)$y, dp = DisplayPars(color = "black", lwd=2,lty = "solid"))
cop <- makeGenericArray(intensity  = cn, probeStart = probestart, 
                        trackOverlay =  seg, dp = DisplayPars(size=3, color = "seagreen", type="dot"))
gdPlot(cop)
}
\keyword{hplot}
