\name{makeGenericArray}
\alias{makeGenericArray}
\title{Creates an object of class GenericArray}
\description{
  Creates an object of class Generic Array representing microarray data.  This could be gene expression, array CGH, etc.
}
\usage{
makeGenericArray(intensity, probeStart, probeEnd, trackOverlay, dp = NULL)
}
\arguments{
  \item{intensity}{Matrix of intensities, probes in the rows, samples in the columns}
  \item{probeStart}{Vector of start positions for the probes}
  \item{probeEnd}{Vector of end positions for probes (optional)}
  \item{trackOverlay}{Object of class TrackOverlay, needs to be added if overlays should be plotted as well}
  \item{dp}{Object of class DisplayPars which handles the display parameters for plotting}
}
\value{
Object of class GenericArray
}
\references{BMC bioinformatics 2009}
\author{Jim Bullard and Steffen Durinck}
\seealso{gdPlot}
\examples{
showClass("GenericArray")
}
\keyword{hplot}
