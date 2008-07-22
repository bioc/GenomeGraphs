\name{makeTitle}
\alias{makeTitle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates an object of class Title}
\description{
  Creates an object of class Title which can be used to add a title to the plot
}
\usage{
makeTitle(text, cex, color, size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{text}{The text that will make up the title}
  \item{cex}{Font size of the title}
  \item{color}{Font color of the title}
  \item{size}{Size of the viewport in which the title resides}
}
\value{
Object of class Title
}
\references{ ~put references to the literature/web site here ~ }
\author{Steffen Durinck and Jim Bullard}
\seealso{\code{\link{gdPlot}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (text, cex, color, size) 
{
    dp <- getClass("Title")@prototype@dp
    if (!missing(cex)) 
        setPar(dp, "cex", cex)
    if (!missing(color)) 
        setPar(dp, "color", color)
    if (!missing(size)) 
        setPar(dp, "size", size)
    new("Title", title = text, dp = dp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
