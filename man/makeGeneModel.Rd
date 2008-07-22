\name{makeGeneModel}
\alias{makeGeneModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Creates an object of class GeneModel}
\description{
Creates an object of class GeneModel representing a custom annotation or gene model
}
\usage{
makeGeneModel(start, end, chromosome, dp = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{start}{Vector of start positions for exons}
  \item{end}{Vector of end positions for exons}
  \item{chromosome}{chromosome name}
  \item{dp}{Display parametes represented as an object of class DisplayPars}
}
\value{
Object of class GeneModel
}
\references{ ~put references to the literature/web site here ~ }
\author{Steffen Durinck and Jim Bullard}
\seealso{\code{\link{DisplayPars}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (start, end, chromosome, dp = NULL) 
{
    if (is.null(dp)) 
        dp <- getClass("GeneModel")@prototype@dp
    new("GeneModel", exonStart = start, exonEnd = end, dp = dp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
