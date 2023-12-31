\name{makeTextOverlay}
\alias{makeTextOverlay}
\title{Create objects of class TextOverlay}
\description{
Create objects of class TextOverlay
}
\usage{
makeTextOverlay(text, xpos, ypos, region = NULL, coords = c("genomic", "absolute"), dp = NULL)
}
\arguments{
  \item{text}{The text to plot}
  \item{xpos}{The xposition of the text}
  \item{ypos}{The yposition of the text}
  \item{region}{Region}
  \item{coords}{Coordinates}
  \item{dp}{The display parameters}
}
\value{
  Returns class of TextOverlay
}
\examples{
data("exampleData", package="GenomeGraphs")
seg <- makeSegmentation(segStart[[1]], segEnd[[1]], segments[[1]], 
                        dp = DisplayPars(color = "black", lwd=2,lty = "solid"))
cop <- makeGenericArray(intensity  = cn, probeStart = probestart, 
                        trackOverlay =  seg, dp = DisplayPars(size=3, color = "seagreen", type="dot"))
gdPlot(cop, overlay = makeTextOverlay("Overlay Text", xpos = .5, ypos = .5, coords = "absolute"))
}
\keyword{hplot}
