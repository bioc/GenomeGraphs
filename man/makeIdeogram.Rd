\name{makeIdeogram}
\alias{makeIdeogram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates object of class Ideogram}
\description{
Creates object of class Ideogram
}
\usage{
makeIdeogram(chromosome, dp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chromosome}{Chromosome to represent (currently human only) }
  \item{dp}{Display parameters such as color and size}
}
\details{
 
}
\value{
Object of class Ideogram
}
\references{ ~put references to the literature/web site here ~ }
\author{Jim Bullard and Steffen Durinck}
\seealso{ \code{\link{gdPlot}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (chromosome, dp = NULL) 
{
    if (missing(chromosome)) 
        stop("Need to specify chromosome for creating an Ideogram")
    if (is.numeric(chromosome)) {
        chromosome = as.character(chromosome)
    }
    if (is.null(dp)) 
        dp <- getClass("Ideogram")@prototype@dp
    new("Ideogram", chromosome = chromosome, dp = dp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
