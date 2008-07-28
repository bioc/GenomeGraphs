\name{makeRectangleOverlay}
\alias{makeRectangleOverlay}
\title{Create a rectangular overlay}
\description{
  Construct ractangular overlays.
}
\usage{
makeRectangleOverlay(start, end, region = NULL, coords = c("genomic", "absolute"), dp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{Start position in coords coordinates }
  \item{end}{End position in coords coordinates }
  \item{region}{Which tracks to span, or the y (vertical range) }
  \item{coords}{ Which coordinate system to use, if absolute then the
    range is from 0,1 and region become the y coordinates }
  \item{dp}{The display parameters }
}
\details{
  The rectangular overlay can be used to plot overlays in either genomic
  or absolute coordinates. If coordinates are absolute then the region
  argument becomes the y arguments.
}
\value{
  An object of class RectangleOverlay
}
\examples{
data("exampleData", package = "GenomeGraphs")
cop <- makeGenericArray(intensity  = cn, probeStart = probestart, 
                        dp = DisplayPars(size=3, color = "seagreen", type="dot"))
gdPlot(list(makeGenomeAxis(), cop), overlays =
       makeRectangleOverlay(start = 180350000, end = 180350000 + 1e5, dp = DisplayPars(alpha = .3)))
}
\keyword{hplot}
