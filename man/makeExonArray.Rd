\name{makeExonArray}
\alias{makeExonArray}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates and object of class ExonArray}
\description{
Creates an object of class ExonArray, representing exon array microarray data
}
\usage{
makeExonArray(intensity, probeStart, probeEnd, probeId, nProbes, displayProbesets = FALSE, dp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{intensity}{Matrix of intentsities, probes in the rows, samples in the columns }
  \item{probeStart}{Vector of probe start positions}
  \item{probeEnd}{Vector of probe end positions (optional)}
  \item{probeId}{Character vector containing the probe identifiers}
  \item{nProbes}{Vector indicating how many probes are in each probeset}
  \item{displayProbesets}{Logical indicating if the probeset idenifier should be displayed or not}
  \item{dp}{Object of class DisplayPars to set the display parameters}
}
\value{
 Object of ExonArray class
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\seealso{\code{\link{gdPlot}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (intensity, probeStart, probeEnd, probeId, nProbes, 
    displayProbesets = FALSE, dp = NULL) 
{
    pt <- getClass("ExonArray")@prototype
    if (is.null(dp)) 
        dp <- pt@dp
    if (missing(probeEnd)) 
        probeEnd <- pt@probeEnd
    if (missing(probeId)) 
        probeId <- pt@probeId
    if (missing(nProbes)) 
        nProbes <- pt@nProbes
    if (is.null(dp)) 
        dp <- getClass("ExonArray")@prototype@dp
    new("ExonArray", intensity = intensity, probeStart = probeStart, 
        probeEnd = probeEnd, probeId = probeId, nProbes = nProbes, 
        displayProbesets = displayProbesets, dp = dp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
