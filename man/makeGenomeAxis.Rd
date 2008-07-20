\name{makeGenomeAxis}
\alias{makeGenomeAxis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates an object of class GenomeAxis}
\description{
Creates an object of class GenomeAxis, representing a genome coordinate axis.
}
\usage{
makeGenomeAxis(add53 = FALSE, add35 = FALSE, littleTicks = FALSE, dp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{add53}{Add a 5 to 3 prime label}
  \item{add35}{Add a 3 to 5 prime label}
  \item{littleTicks}{Add smaller ticks between larger ticks}
  \item{dp}{Set the display parameters see DisplayPars}
}
\value{
Object of class GenomeAxis
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\seealso{\code{\link{DisplayPars}},\code{\link{gdPlot}} }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (add53 = FALSE, add35 = FALSE, littleTicks = FALSE, 
    dp = NULL) 
{
    if (is.null(dp)) 
        dp <- getClass("GenomeAxis")@prototype@dp
    new("GenomeAxis", add53 = add53, add35 = add35, dp = dp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
