\name{makeBaseTrack}
\alias{makeBaseTrack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates an object of class BaseTrack}
\description{
Creates an object of class BaseTrack, which can represent many datasets containing a base given by a vector of positions and a corresponding vector with values for these base positions
}
\usage{
makeBaseTrack(base, value, strand, trackOverlay, dp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{base}{Numeric vector of base positions}
  \item{value}{Numeric vector with values for these base positions}
  \item{strand}{Character either + or - representing the strand}
  \item{trackOverlay}{Object of class TrackOverlay, used when overlays are needed to be drawn}
  \item{dp}{Object of class DisplayPars representing the display parameters of the plot}
}
\value{
Object of class BaseTrack
}
\references{ ~put references to the literature/web site here ~ }
\author{Jim Bullard and Steffen Durinck}
\seealso{\code{\link{DisplayPars}}, \code{\link{gdPlot}} }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (base, value, strand, segmentation, dp = NULL) 
{
    pt <- getClass("BaseTrack")@prototype
    if (is.null(dp)) 
        dp <- pt@dp
    if (missing(strand)) 
        strand <- pt@strand
    if (missing(segmentation)) 
        segmentation <- pt@segmentation
    if (missing(base)) 
        stop("Need base argument to know the base positions to plot the data on the genome")
    if (missing(value)) 
        stop("Need value argument")
    new("BaseTrack", base = base, value = value, strand = strand, 
        dp = dp, segmentation = segmentation)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
