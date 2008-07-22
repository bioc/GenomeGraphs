\name{makeLegend}
\alias{makeLegend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates an object of class Legend}
\description{
  Creates an object of class Legend which can be used to plot a legend
}
\usage{
makeLegend(text, fill, cex)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{text}{Vector of characters representing the legend}
  \item{fill}{Vector of colors to fill the legend boxes}
  \item{cex}{Font size of the legend}
}
\value{
 Object of class Legend 
}
\references{ ~put references to the literature/web site here ~ }
\author{Jim Bullard and Steffen Durinck}
\seealso{ See Also as \code{\link{gdPlot}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (text, fill, cex) 
{
    dp <- getClass("Legend")@prototype@dp
    if (!missing(cex)) 
        setPar(dp, "cex", cex)
    if (!missing(fill)) 
        setPar(dp, "color", fill)
    new("Legend", legend = text, dp = dp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
