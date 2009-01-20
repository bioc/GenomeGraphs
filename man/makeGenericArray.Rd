\name{makeGenericArray}
\alias{makeGenericArray}
\title{Creates an object of class GenericArray}
\description{
  Creates an object of class Generic Array representing microarray data.  This could be gene expression, array CGH, etc.
}
\usage{
makeGenericArray(intensity, probeStart, probeEnd, segmentation, dp = NULL)
}
\arguments{
  \item{intensity}{Matrix of intensities, probes in the rows, samples in the columns}
  \item{probeStart}{Vector of start positions for the probes}
  \item{probeEnd}{Vector of end positions for probes (optional)}
  \item{segmentation}{Object of class Segmentation, needs to be added if segments should be plotted as well}
  \item{dp}{Object of class DisplayPars which handles the display parameters for plotting}
}
\value{
Object of class GenericArray
}
\references{BMC bioinformatics 2009}
\author{Jim Bullard and Steffen Durinck}
\seealso{gdPlot}
\examples{
function (intensity, probeStart, probeEnd, segmentation, dp = NULL) 
{
    pt <- getClass("GenericArray")@prototype
    if (is.null(dp)) 
        dp <- pt@dp
    if (missing(probeEnd)) 
        probeEnd <- pt@probeEnd
    if (missing(segmentation)) 
        segmentation <- pt@segmentation
    if (missing(probeStart)) 
        stop("Need probeStart argument to know where to plot the data on the genome")
    new("GenericArray", intensity = intensity, probeStart = probeStart, 
        probeEnd = probeEnd, dp = dp, segmentation = segmentation)
  }
}
\keyword{hplot}
