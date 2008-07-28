\name{makeSegmentation}
\alias{makeSegmentation}
\title{Create objects of class segmentation}
\description{
  Construct objects of class segmentation
}
\usage{
makeSegmentation(start, end, value, dp = NULL)
}
\arguments{
  \item{start}{Either a list or a vector. If it is a list then it is a
    list of vectors of start position (this is the way it is represented
  in the segmentation class) If it is a vector it is a vector of start
  positions.}
  \item{end}{Same as start, but the corresponding end positions.}
  \item{value}{The y value of the segmentation, ie. segments(start[i],
    value[i], end[i], value[i])}
  \item{dp}{ The Display parameters.}
}
\value{
  An object of class Segmentation
}
\examples{
data("exampleData", package="GenomeGraphs")
seg <- makeSegmentation(segStart, segEnd, segments, 
                        dp = DisplayPars(color = "black", lwd=2,lty = "solid"))
cop <- makeGenericArray(intensity  = cn, probeStart = probestart, 
                        segmentation =  seg, dp = DisplayPars(size=3, color = "seagreen", type="dot"))
gdPlot(cop)
}
\keyword{hplot}
