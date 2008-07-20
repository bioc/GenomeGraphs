\name{makeGenericArray}
\alias{makeGenericArray}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates an object of class GenericArray}
\description{
  Creates an object of class Generic Array representing microarray data.  This could be gene expression, array CGH, etc.
}
\usage{
makeGenericArray(intensity, probeStart, probeEnd, segmentation, dp = NULL)
}
%- maybe also 'usage' for other objects documented here.
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
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
